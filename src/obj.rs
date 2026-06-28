//! Wavefront OBJ loading. `v` lines become vertices; `f` lines fan-triangulate
//! into faces; everything else (normals, texcoords, groups, comments) is ignored.

use std::path::Path;

use anyhow::{anyhow, Context, Result};
use spectral3d::Mesh;

/// Read an OBJ file and parse it into a triangle mesh.
pub fn load(path: &Path) -> Result<Mesh> {
    let bytes = std::fs::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    parse(&bytes).with_context(|| format!("{} isn't a valid OBJ", path.display()))
}

/// Parse Wavefront OBJ bytes into a triangle mesh.
///
/// Positive face indices are 1-based, negatives count back from the most recent
/// vertex, and a `v/vt/vn` triple keeps only the position. Polygons
/// fan-triangulate around their first vertex.
fn parse(bytes: &[u8]) -> Result<Mesh> {
    let text = String::from_utf8_lossy(bytes);
    let mut vertices: Vec<[f64; 3]> = Vec::new();
    let mut faces: Vec<[u32; 3]> = Vec::new();

    for (lineno, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split_whitespace();
        match it.next() {
            Some("v") => {
                let mut p = [0f64; 3];
                for c in p.iter_mut() {
                    let tok = it
                        .next()
                        .ok_or_else(|| anyhow!("line {}: vertex needs 3 coords", lineno + 1))?;
                    *c = tok
                        .parse::<f64>()
                        .map_err(|_| anyhow!("line {}: bad float '{tok}'", lineno + 1))?;
                }
                vertices.push(p);
            }
            Some("f") => {
                let mut idx: Vec<u32> = Vec::new();
                for tok in it {
                    let first = tok.split('/').next().unwrap_or("");
                    let i = first
                        .parse::<i64>()
                        .map_err(|_| anyhow!("line {}: bad index '{tok}'", lineno + 1))?;
                    let resolved = if i > 0 {
                        i - 1
                    } else if i < 0 {
                        vertices.len() as i64 + i
                    } else {
                        return Err(anyhow!("line {}: zero index", lineno + 1));
                    };
                    if resolved < 0 || resolved >= vertices.len() as i64 {
                        return Err(anyhow!("line {}: index {i} out of range", lineno + 1));
                    }
                    idx.push(resolved as u32);
                }
                if idx.len() < 3 {
                    return Err(anyhow!("line {}: face needs >= 3 vertices", lineno + 1));
                }
                for k in 1..idx.len() - 1 {
                    faces.push([idx[0], idx[k], idx[k + 1]]);
                }
            }
            _ => {}
        }
    }

    if vertices.is_empty() || faces.is_empty() {
        return Err(anyhow!("no geometry"));
    }
    Ok(Mesh { vertices, faces })
}
