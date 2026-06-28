//! Mesh repair lives out here, outside the spectral3d core, so the library
//! never silently alters geometry. The pipeline is deterministic and runs in a
//! fixed order:
//!
//!   1. drop degenerate faces  (zero-area, a repeated vertex) — shape-preserving
//!   2. weld coincident verts  (merge duplicates on a grid)   — shape-preserving
//!   3. drop duplicate faces   (same vertex set)              — redundant
//!   4. fill holes             (fan-cap simple boundary loops) — adds geometry
//!   5. orient shells          (flip inward shells positive)  — winding only
//!
//! The pipeline is destructive by design (a DCC-style "make solid"). On top of
//! the shape-preserving cleanup above it also:
//!
//!   a. strip non-manifold faces  (cut every face on a 3+ edge) — reshapes
//!   b. drop debris components     (sub-1e-4-area floaters)      — removes
//!   c. force-close all loops      (centroid-fan even frayed ones) — invents
//!
//! and iterates a..c to a fixed point. It trades fidelity for a guaranteed
//! closed solid. It lives in the CLI, never in the spectral3d core, and writes a
//! separate file the user reviews, so the security boundary is untouched.

use std::collections::{HashMap, HashSet};

pub struct RepairReport {
    pub degenerate_dropped: usize,
    pub welded: usize,
    pub duplicate_faces_dropped: usize,
    pub nonmanifold_faces_removed: usize,
    pub components_dropped: usize,
    pub holes_filled: usize,
    pub cap_faces_added: usize,
    pub shells_flipped: usize,
    pub open_edges_left: usize,
    pub nonmanifold_edges_left: usize,
    pub closed: bool,
}

fn norm(a: u32, b: u32) -> (u32, u32) {
    if a < b {
        (a, b)
    } else {
        (b, a)
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn nondegenerate(f: &[u32; 3]) -> bool {
    f[0] != f[1] && f[1] != f[2] && f[0] != f[2]
}

/// Map each vertex to a representative, merging anything in the same grid cell
/// (1e-6 of the bounding-box diagonal), then UN-merging any group that the weld
/// turned non-manifold. A model can split a vertex on purpose at a pinch point
/// so two sheets meet at a point without fusing; welding that blindly would
/// forge a non-manifold edge and break a mesh that was fine. The guard keeps the
/// safe seams and leaves deliberate splits alone.
fn weld_map(verts: &[[f64; 3]], faces: &[[u32; 3]]) -> Vec<u32> {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in verts {
        for k in 0..3 {
            lo[k] = lo[k].min(v[k]);
            hi[k] = hi[k].max(v[k]);
        }
    }
    let diag = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let tol = if diag > 0.0 { diag * 1e-6 } else { 1e-9 };
    let cell = |x: f64| (x / tol).round() as i64;

    let mut first: HashMap<(i64, i64, i64), u32> = HashMap::new();
    let mut rep = vec![0u32; verts.len()];
    for (i, v) in verts.iter().enumerate() {
        let key = (cell(v[0]), cell(v[1]), cell(v[2]));
        rep[i] = *first.entry(key).or_insert(i as u32);
    }

    // Un-weld any group that a merge turned non-manifold, until stable. We only
    // ever split groups back apart, so this shrinks merges monotonically and
    // terminates.
    loop {
        let mut count: HashMap<(u32, u32), i32> = HashMap::new();
        for &[i, j, k] in faces {
            let (a, b, c) = (rep[i as usize], rep[j as usize], rep[k as usize]);
            if a == b || b == c || a == c {
                continue;
            }
            for (x, y) in [(a, b), (b, c), (c, a)] {
                *count.entry(norm(x, y)).or_insert(0) += 1;
            }
        }
        let mut bad: HashSet<u32> = HashSet::new();
        for (&(x, y), &c) in &count {
            if c > 2 {
                bad.insert(x);
                bad.insert(y);
            }
        }
        if bad.is_empty() {
            break;
        }
        let mut changed = false;
        for (i, r) in rep.iter_mut().enumerate() {
            if *r != i as u32 && bad.contains(r) {
                *r = i as u32;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }
    rep
}

/// Union-find over faces linked by a shared edge, returning each face's
/// component root.
fn components(faces: &[[u32; 3]]) -> Vec<usize> {
    let n = faces.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }
    let mut seen: HashMap<(u32, u32), usize> = HashMap::new();
    for (fi, &[i, j, k]) in faces.iter().enumerate() {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            let key = norm(a, b);
            if let Some(&fj) = seen.get(&key) {
                let (ra, rb) = (find(&mut parent, fi), find(&mut parent, fj));
                if ra != rb {
                    parent[ra] = rb;
                }
            } else {
                seen.insert(key, fi);
            }
        }
    }
    (0..n).map(|f| find(&mut parent, f)).collect()
}

/// Fan-triangulate every simple boundary loop. Returns the cap faces and the
/// loop count; frayed (non-simple) boundaries are left for the audit to report.
fn cap_faces(faces: &[[u32; 3]]) -> (Vec<[u32; 3]>, usize) {
    let mut count: HashMap<(u32, u32), i32> = HashMap::new();
    for &[i, j, k] in faces {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            *count.entry(norm(a, b)).or_insert(0) += 1;
        }
    }
    let mut succ: HashMap<u32, u32> = HashMap::new();
    for &[i, j, k] in faces {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            if count[&norm(a, b)] == 1 {
                succ.insert(a, b);
            }
        }
    }

    let mut used: HashSet<u32> = HashSet::new();
    let mut fill = Vec::new();
    let mut holes = 0;
    let mut starts: Vec<u32> = succ.keys().copied().collect();
    starts.sort_unstable();
    for s in starts {
        if used.contains(&s) {
            continue;
        }
        let mut chain = Vec::new();
        let mut cur = s;
        while !used.contains(&cur) {
            used.insert(cur);
            chain.push(cur);
            match succ.get(&cur) {
                Some(&n) => cur = n,
                None => break,
            }
        }
        if cur == s && chain.len() >= 3 {
            // Fan wound opposite to the boundary so the cap supplies the missing
            // half-edges and closes the shell instead of doubling it.
            let v0 = chain[0];
            for w in 1..chain.len() - 1 {
                fill.push([v0, chain[w + 1], chain[w]]);
            }
            holes += 1;
        }
    }
    (fill, holes)
}

/// Flip every shell whose signed volume came out negative, so all shells share
/// one (outward) winding. A closed shell's signed volume is reference-free, so
/// the origin works as the apex.
fn orient_shells(verts: &[[f64; 3]], faces: &mut [[u32; 3]]) -> usize {
    if faces.is_empty() {
        return 0;
    }
    let root = components(faces);
    let mut vol: HashMap<usize, f64> = HashMap::new();
    for (fi, f) in faces.iter().enumerate() {
        let (a, b, c) = (
            verts[f[0] as usize],
            verts[f[1] as usize],
            verts[f[2] as usize],
        );
        *vol.entry(root[fi]).or_insert(0.0) += dot(a, cross(b, c)) / 6.0;
    }
    let inward: HashSet<usize> = vol
        .iter()
        .filter(|&(_, &v)| v < 0.0)
        .map(|(&r, _)| r)
        .collect();
    if inward.is_empty() {
        return 0;
    }
    for (fi, f) in faces.iter_mut().enumerate() {
        if inward.contains(&root[fi]) {
            f.swap(1, 2);
        }
    }
    inward.len()
}

fn edge_audit(faces: &[[u32; 3]]) -> (usize, usize) {
    let mut count: HashMap<(u32, u32), i32> = HashMap::new();
    for &[i, j, k] in faces {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            *count.entry(norm(a, b)).or_insert(0) += 1;
        }
    }
    (
        count.values().filter(|&&c| c == 1).count(),
        count.values().filter(|&&c| c > 2).count(),
    )
}

/// Cut every face that touches a non-manifold edge (3+ faces). One sweep does
/// it: dropping faces can only lower an edge's count, never raise it, so after
/// removing all faces on the offending edges none remain above two. Aggressive
/// on purpose — it leaves clean holes for the cap stage to close.
fn strip_nonmanifold(faces: &mut Vec<[u32; 3]>) -> usize {
    let mut count: HashMap<(u32, u32), i32> = HashMap::new();
    for &[i, j, k] in faces.iter() {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            *count.entry(norm(a, b)).or_insert(0) += 1;
        }
    }
    let bad: HashSet<(u32, u32)> = count
        .iter()
        .filter(|&(_, &c)| c > 2)
        .map(|(&e, _)| e)
        .collect();
    if bad.is_empty() {
        return 0;
    }
    let before = faces.len();
    faces.retain(|&[i, j, k]| {
        !(bad.contains(&norm(i, j)) || bad.contains(&norm(j, k)) || bad.contains(&norm(k, i)))
    });
    before - faces.len()
}

fn tri_area(verts: &[[f64; 3]], f: &[u32; 3]) -> f64 {
    let (a, b, c) = (
        verts[f[0] as usize],
        verts[f[1] as usize],
        verts[f[2] as usize],
    );
    let ab = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let ac = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let cr = cross(ab, ac);
    0.5 * (cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]).sqrt()
}

/// Drop floating fragments under 1e-4 of total area — the slivers that stripping
/// non-manifold faces tends to shed. The threshold sits far below any genuine
/// part (a wheel is whole percents of the model), so multi-part solids survive.
fn drop_small_components(verts: &[[f64; 3]], faces: &mut Vec<[u32; 3]>) -> usize {
    if faces.is_empty() {
        return 0;
    }
    let root = components(faces);
    let mut comp_area: HashMap<usize, f64> = HashMap::new();
    let mut total = 0.0;
    for (fi, f) in faces.iter().enumerate() {
        let s = tri_area(verts, f);
        *comp_area.entry(root[fi]).or_insert(0.0) += s;
        total += s;
    }
    if total <= 0.0 {
        return 0;
    }
    // Never drop the biggest component, so a model made entirely of small parts
    // (a pile of leaf cards) survives as geometry instead of vanishing.
    let biggest = comp_area
        .iter()
        .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
        .map(|(&r, _)| r);
    let floor = total * 1e-4;
    let drop: HashSet<usize> = comp_area
        .iter()
        .filter(|&(&r, &s)| s < floor && Some(r) != biggest)
        .map(|(&r, _)| r)
        .collect();
    if drop.is_empty() {
        return 0;
    }
    let keep: Vec<[u32; 3]> = faces
        .iter()
        .enumerate()
        .filter(|(fi, _)| !drop.contains(&root[*fi]))
        .map(|(_, f)| *f)
        .collect();
    *faces = keep;
    drop.len()
}

/// Close every remaining boundary loop, frayed ones included, by tracing the
/// directed boundary half-edges into cycles and fanning each off a fresh
/// centroid hub. Branching vertices keep all their boundary edges (a Vec), and
/// the walk consumes the smallest-indexed one each step for determinism. Adds
/// one vertex per loop, so it genuinely reshapes — hence force-only. Returns
/// (hub verts to append from `base`, cap faces, loops closed).
fn cap_faces_force(
    verts: &[[f64; 3]],
    faces: &[[u32; 3]],
    base: u32,
) -> (Vec<[f64; 3]>, Vec<[u32; 3]>, usize) {
    let mut count: HashMap<(u32, u32), i32> = HashMap::new();
    for &[i, j, k] in faces {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            *count.entry(norm(a, b)).or_insert(0) += 1;
        }
    }
    let mut succ: HashMap<u32, Vec<u32>> = HashMap::new();
    for &[i, j, k] in faces {
        for (a, b) in [(i, j), (j, k), (k, i)] {
            if count[&norm(a, b)] == 1 {
                succ.entry(a).or_default().push(b);
            }
        }
    }
    // pop() yields the smallest target, so traversal is deterministic
    for v in succ.values_mut() {
        v.sort_unstable_by(|a, b| b.cmp(a));
    }
    let mut sources: Vec<u32> = succ.keys().copied().collect();
    sources.sort_unstable();

    let mut extra: Vec<[f64; 3]> = Vec::new();
    let mut fill: Vec<[u32; 3]> = Vec::new();
    let mut holes = 0;
    for &s in &sources {
        while succ.get(&s).is_some_and(|v| !v.is_empty()) {
            // Walk a closed trail, splitting it into simple cycles as we go: the
            // moment the path revisits a vertex, the segment since that vertex
            // is a self-contained loop — cap it and pop it off, then keep
            // walking. This keeps every fan over distinct vertices, so a hub
            // spoke never lands on three faces (the bowtie that forged
            // non-manifold edges before).
            let mut path: Vec<u32> = Vec::new();
            let mut pos: HashMap<u32, usize> = HashMap::new();
            let mut cur = s;
            while let Some(nxt) = succ.get_mut(&cur).and_then(|v| v.pop()) {
                pos.insert(cur, path.len());
                path.push(cur);
                cur = nxt;
                if let Some(&p) = pos.get(&cur) {
                    if path.len() - p >= 3 {
                        let cycle = &path[p..];
                        let mut c = [0.0; 3];
                        for &v in cycle {
                            for (t, ct) in c.iter_mut().enumerate() {
                                *ct += verts[v as usize][t];
                            }
                        }
                        for ct in &mut c {
                            *ct /= cycle.len() as f64;
                        }
                        let hub = base + extra.len() as u32;
                        extra.push(c);
                        let n = cycle.len();
                        for w in 0..n {
                            let a = cycle[w];
                            let b = cycle[(w + 1) % n];
                            // boundary half-edge a->b; the cap supplies b->a
                            fill.push([hub, b, a]);
                        }
                        holes += 1;
                    }
                    for &v in &path[p..] {
                        pos.remove(&v);
                    }
                    path.truncate(p);
                }
            }
        }
    }
    (extra, fill, holes)
}

/// Run the full pipeline. Faces come back referencing a fresh, compacted vertex
/// list (returned alongside), so the output OBJ carries no orphan vertices.
pub fn repair_mesh(
    verts: &[[f64; 3]],
    faces0: &[[u32; 3]],
) -> (Vec<[f64; 3]>, Vec<[u32; 3]>, RepairReport) {
    // 1. drop degenerate faces
    let mut faces: Vec<[u32; 3]> = faces0.iter().copied().filter(nondegenerate).collect();
    let mut degenerate_dropped = faces0.len() - faces.len();

    // 2. weld coincident vertices (remap face indices to representatives)
    let rep = weld_map(verts, &faces);
    let welded = (0..verts.len()).filter(|&i| rep[i] != i as u32).count();
    for f in &mut faces {
        for k in 0..3 {
            f[k] = rep[f[k] as usize];
        }
    }
    let before = faces.len();
    faces.retain(nondegenerate);
    degenerate_dropped += before - faces.len();

    // 3. drop duplicate faces (same vertex set, any winding)
    let before = faces.len();
    let mut seen: HashSet<[u32; 3]> = HashSet::new();
    faces.retain(|f| {
        let mut key = *f;
        key.sort_unstable();
        seen.insert(key)
    });
    let duplicate_faces_dropped = before - faces.len();

    // Force-cap may invent centroid hubs, so carry a working vertex list that
    // can grow past the input.
    let mut work_verts: Vec<[f64; 3]> = verts.to_vec();

    let mut nonmanifold_faces_removed = 0;
    let mut components_dropped = 0;
    let mut holes_filled = 0;
    let mut cap_faces_added = 0;

    // Iterate to a fixed point: cut non-manifold faces, shed debris, then cap
    // every loop. Capping a frayed junction can expose a few fresh non-manifold
    // edges, so the next round's cut clears them. Bounded, and it breaks the
    // moment the mesh is a clean closed manifold. On an already-clean mesh every
    // stage is a no-op and the first round breaks out at once.
    for _ in 0..6 {
        nonmanifold_faces_removed += strip_nonmanifold(&mut faces);
        components_dropped += drop_small_components(&work_verts, &mut faces);

        let (caps, h) = cap_faces(&faces);
        cap_faces_added += caps.len();
        holes_filled += h;
        faces.extend(caps);

        let (extra, fcaps, fh) = cap_faces_force(&work_verts, &faces, work_verts.len() as u32);
        work_verts.extend(extra);
        cap_faces_added += fcaps.len();
        holes_filled += fh;
        faces.extend(fcaps);

        let (oe, nm) = edge_audit(&faces);
        if oe == 0 && nm == 0 {
            break;
        }
    }

    // 5. orient shells consistently
    let shells_flipped = orient_shells(&work_verts, &mut faces);

    // audit + compact onto a fresh vertex list
    let (open_edges_left, nonmanifold_edges_left) = edge_audit(&faces);
    let mut remap: HashMap<u32, u32> = HashMap::new();
    let mut out_verts: Vec<[f64; 3]> = Vec::new();
    for f in &mut faces {
        for slot in f.iter_mut() {
            let old = *slot;
            *slot = *remap.entry(old).or_insert_with(|| {
                out_verts.push(work_verts[old as usize]);
                (out_verts.len() - 1) as u32
            });
        }
    }

    let report = RepairReport {
        degenerate_dropped,
        welded,
        duplicate_faces_dropped,
        nonmanifold_faces_removed,
        components_dropped,
        holes_filled,
        cap_faces_added,
        shells_flipped,
        open_edges_left,
        nonmanifold_edges_left,
        // An empty mesh is vacuously closed — don't call that a win.
        closed: !faces.is_empty() && open_edges_left == 0 && nonmanifold_edges_left == 0,
    };
    (out_verts, faces, report)
}

/// Serialize the repaired mesh as a clean OBJ (vertices + faces only), with a
/// header recording exactly what the pipeline did.
pub fn to_obj(verts: &[[f64; 3]], faces: &[[u32; 3]], r: &RepairReport) -> String {
    let mut out = String::new();
    out.push_str(&format!(
        "# repaired by spectral3d-cli v{}\n",
        env!("CARGO_PKG_VERSION")
    ));
    out.push_str(&format!(
        "# dropped {} degenerate face(s), welded {} vertex/vertices, dropped {} duplicate face(s)\n",
        r.degenerate_dropped, r.welded, r.duplicate_faces_dropped
    ));
    if r.nonmanifold_faces_removed > 0 || r.components_dropped > 0 {
        out.push_str(&format!(
            "# destructive: cut {} non-manifold face(s), dropped {} debris component(s)\n",
            r.nonmanifold_faces_removed, r.components_dropped
        ));
    }
    out.push_str(&format!(
        "# filled {} hole(s) (+{} faces), flipped {} shell(s)\n",
        r.holes_filled, r.cap_faces_added, r.shells_flipped
    ));
    out.push_str(&format!(
        "# result: closed={}, {} open edge(s), {} non-manifold edge(s)\n",
        r.closed, r.open_edges_left, r.nonmanifold_edges_left
    ));
    out.push('\n');
    for v in verts {
        out.push_str(&format!("v {} {} {}\n", v[0], v[1], v[2]));
    }
    for f in faces {
        out.push_str(&format!("f {} {} {}\n", f[0] + 1, f[1] + 1, f[2] + 1));
    }
    out
}
