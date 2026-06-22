mod repair;

use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Context, Result};
use clap::{Parser, Subcommand};
use serde_json::json;
use spectral3d::{Helper, Mesh, SpectralParams, N_FEATURES};

#[derive(Parser)]
#[command(name = "spectral3d-cli", about = "Spectral identity and mesh repair for 3D models")]
#[command(args_conflicts_with_subcommands = true)]
struct Cli {
    #[command(subcommand)]
    cmd: Option<Commands>,
    /// OBJ to register — the default action when no subcommand is given.
    obj: Option<PathBuf>,
}

#[derive(Subcommand)]
enum Commands {
    /// Register an OBJ, writing its identity to <name>-register.json (the default action).
    Register { obj: PathBuf },
    /// Repair a mesh into a closed solid, writing <name>-repair.obj alongside it.
    Repair { path: PathBuf },
    /// Verify an OBJ against a registration record (the <name>-register.json from `register`).
    Verify { obj: PathBuf, record: PathBuf },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.cmd {
        Some(Commands::Register { obj }) => register(obj),
        Some(Commands::Repair { path }) => repair(path),
        Some(Commands::Verify { obj, record }) => verify(obj, record),
        // No subcommand: the bare path defaults to register.
        None => match cli.obj {
            Some(obj) => register(obj),
            None => Err(anyhow!("give an OBJ to register, or `repair <path>` to fix one")),
        },
    }
}

fn register(obj: PathBuf) -> Result<()> {
    let data = fs::read(&obj).with_context(|| format!("failed to read {}", obj.display()))?;
    let params = SpectralParams::default();
    let (hash, helper) = spectral3d::register(&data, &params).map_err(|e| anyhow!("{e}"))?;

    let offsets: Vec<f64> = helper.offsets.to_vec();
    let record = json!({"hash": hash, "helper": {"offsets": offsets}});

    let stem = obj
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("can't derive an output name from {}", obj.display()))?;
    let out_path = obj.with_file_name(format!("{stem}-register.json"));
    let body = serde_json::to_string_pretty(&record).expect("json value always serializes");
    fs::write(&out_path, format!("{body}\n"))
        .with_context(|| format!("failed to write {}", out_path.display()))?;

    println!("identity {hash}");
    println!("output   {}", out_path.display());
    Ok(())
}

fn verify(obj: PathBuf, record: PathBuf) -> Result<()> {
    let data = fs::read(&obj).with_context(|| format!("failed to read {}", obj.display()))?;
    let rec_bytes =
        fs::read(&record).with_context(|| format!("failed to read {}", record.display()))?;
    let rec: serde_json::Value = serde_json::from_slice(&rec_bytes)
        .with_context(|| format!("{} isn't valid JSON", record.display()))?;

    // Pull the identity and the published offsets back out of the record.
    let expected = rec["hash"]
        .as_str()
        .ok_or_else(|| anyhow!("no \"hash\" string in {}", record.display()))?;
    let offsets_json = rec["helper"]["offsets"]
        .as_array()
        .ok_or_else(|| anyhow!("no \"helper.offsets\" array in {}", record.display()))?;
    let offsets_vec: Vec<f64> = offsets_json
        .iter()
        .map(|v| v.as_f64().ok_or_else(|| anyhow!("a helper offset isn't a number")))
        .collect::<Result<_>>()?;
    let offsets: [f64; N_FEATURES] = offsets_vec
        .try_into()
        .map_err(|v: Vec<f64>| anyhow!("expected {N_FEATURES} offsets, got {}", v.len()))?;

    let params = SpectralParams::default();
    let got = spectral3d::verify(&data, &Helper { offsets }, &params).map_err(|e| anyhow!("{e}"))?;

    // Same object up to pose/scale/noise iff the recovered hash equals the stored one.
    let matches = got == expected;
    if matches {
        println!("PASS ✓");
    } else {
        println!("FAIL ✗");
    }

    Ok(())
}

fn repair(path: PathBuf) -> Result<()> {
    let bytes = fs::read(&path).with_context(|| format!("failed to read {}", path.display()))?;
    let mesh = Mesh::parse_obj(&bytes).map_err(|e| anyhow!("{e}"))?;

    let stem = path
        .file_stem()
        .and_then(|s| s.to_str())
        .ok_or_else(|| anyhow!("can't derive an output name from {}", path.display()))?;
    let out_path = path.with_file_name(format!("{stem}-repair.obj"));

    if out_path.exists() && !confirm_overwrite(&out_path)? {
        return Ok(());
    }

    let (verts, faces, r) = repair::repair_mesh(&mesh.vertices, &mesh.faces);
    let out = repair::to_obj(&verts, &faces, &r);
    fs::write(&out_path, out).with_context(|| format!("failed to write {}", out_path.display()))?;

    println!(
        "dropped {} degenerate, welded {}, dropped {} dup faces; cut {} non-manifold, dropped {} debris; filled {} holes (+{}), flipped {} shell(s)",
        r.degenerate_dropped, r.welded, r.duplicate_faces_dropped,
        r.nonmanifold_faces_removed, r.components_dropped,
        r.holes_filled, r.cap_faces_added, r.shells_flipped
    );
    if faces.is_empty() {
        println!("nothing solid survived: the input has no closed-solid geometry (e.g. flat cards)");
    } else if r.closed {
        println!("result is now a closed manifold");
    } else {
        println!(
            "still not closed: {} open edge(s), {} non-manifold edge(s) left",
            r.open_edges_left, r.nonmanifold_edges_left
        );
    }
    println!("{}", out_path.display());
    Ok(())
}

/// Prompt before overwriting an existing output. With no terminal on stdin we
/// bail instead of quietly clobbering whatever's already sitting there.
fn confirm_overwrite(path: &Path) -> Result<bool> {
    use std::io::{IsTerminal, Write};

    if !std::io::stdin().is_terminal() {
        return Err(anyhow!(
            "{} already exists and stdin isn't a terminal — refusing to overwrite",
            path.display()
        ));
    }

    print!("{} exists — overwrite? [y/N] ", path.display());
    std::io::stdout().flush().ok();

    let mut answer = String::new();
    std::io::stdin().read_line(&mut answer)?;
    Ok(matches!(answer.trim(), "y" | "Y" | "yes" | "Yes"))
}
