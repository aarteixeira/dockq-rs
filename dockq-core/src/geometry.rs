//! Geometry kernels. CONTRACT for task #5 (geometry agent). All math in f32 to match
//! the reference (Biopython float32 coords + Cython float32 kernel).
//!
//! Validated against the Python oracle: `operations.residue_distances`,
//! `operations.get_fnat_stats`, and `Bio.SVDSuperimposer.SVDSuperimposer`.

use crate::model::DistMatrix;
use nalgebra::{Matrix3, Vector3};
use rayon::prelude::*;

/// Above this many residue rows, parallelize `residue_distances` over rows.
/// Below it, the serial path wins (rayon fan-out overhead dominates small matrices).
/// Chosen from the microbench in this crate's tests/report (serial faster for the
/// typical few-hundred-residue interfaces; parallel pays off for large all-vs-all).
const PAR_ROW_THRESHOLD: usize = 128;

/// Minimum squared inter-atom distance per residue pair.
///
/// `coords1` / `coords2` are flat atom-coordinate lists (concatenated over residues in
/// order); `atoms_per_res*` say how consecutive atoms belong to each residue.
/// Result `[i][j]` = min over atoms (a in res i, b in res j) of |a-b|^2.
///
/// This replicates the reference Cython kernel **exactly**, including its early-break
/// quirk: `min_d` starts at 100000.0, and the moment an update sets `min_d > 1000.0`
/// the inner atom loops break. Consequence — if the *first* atom pair of a residue pair
/// is farther than sqrt(1000) Å, the stored value is that first pair's squared distance,
/// not the true minimum. Real interface data never hits this (distances are small), but
/// matching it keeps us bit-equivalent to the oracle on arbitrary inputs.
pub fn residue_distances(
    coords1: &[[f32; 3]],
    coords2: &[[f32; 3]],
    atoms_per_res1: &[usize],
    atoms_per_res2: &[usize],
) -> DistMatrix {
    let n_res_i = atoms_per_res1.len();
    let n_res_j = atoms_per_res2.len();

    // Prefix offsets into the flat coord arrays, so each row can be computed
    // independently (needed for the parallel path).
    let mut row_starts = Vec::with_capacity(n_res_i + 1);
    let mut acc = 0usize;
    row_starts.push(0usize);
    for &a in atoms_per_res1 {
        acc += a;
        row_starts.push(acc);
    }
    let mut col_starts = Vec::with_capacity(n_res_j + 1);
    let mut acc = 0usize;
    col_starts.push(0usize);
    for &a in atoms_per_res2 {
        acc += a;
        col_starts.push(acc);
    }

    let mut out = DistMatrix::zeros(n_res_i, n_res_j);

    let fill_row = |i: usize, row: &mut [f32]| {
        let x0 = row_starts[i];
        let x1 = row_starts[i + 1];
        for j in 0..n_res_j {
            let y0 = col_starts[j];
            let y1 = col_starts[j + 1];
            row[j] = min_sq_dist(&coords1[x0..x1], &coords2[y0..y1]);
        }
    };

    if n_res_i >= PAR_ROW_THRESHOLD {
        out.data
            .par_chunks_mut(n_res_j)
            .enumerate()
            .for_each(|(i, row)| fill_row(i, row));
    } else {
        for (i, row) in out.data.chunks_mut(n_res_j).enumerate() {
            fill_row(i, row);
        }
    }

    out
}

/// Minimum squared distance between two atom blocks, replicating the Cython inner
/// loops verbatim (init 100000.0, early break once `min_d > 1000.0` after an update).
#[inline]
fn min_sq_dist(block_a: &[[f32; 3]], block_b: &[[f32; 3]]) -> f32 {
    let mut min_d: f32 = 100000.0;
    'outer: for a in block_a {
        for b in block_b {
            let dx = a[0] - b[0];
            let dy = a[1] - b[1];
            let dz = a[2] - b[2];
            let this_d = dx * dx + dy * dy + dz * dz;
            if this_d < min_d {
                min_d = this_d;
                if min_d > 1000.0 {
                    break 'outer;
                }
            }
        }
    }
    min_d
}

/// Fnat statistics from two squared-distance matrices, thresholded at `threshold_sq`.
/// Returns `(n_shared, n_nonnative, n_native, n_model)`:
/// native = native_d < t; model = model_d < t; shared = model & native;
/// nonnative = model & !native. Replicates the Cython `get_fnat_stats` loop exactly,
/// except `threshold_sq` is already squared by the caller (the Cython version squares
/// its `threshold` argument internally).
pub fn fnat_stats(
    model_d: &DistMatrix,
    native_d: &DistMatrix,
    threshold_sq: f32,
) -> (u32, u32, u32, u32) {
    let mut n_native = 0u32;
    let mut n_model = 0u32;
    let mut n_shared = 0u32;
    let mut n_nonnat = 0u32;

    // Iterate the native matrix shape, exactly as the reference does.
    let rows = native_d.rows;
    let cols = native_d.cols;
    for i in 0..rows {
        for j in 0..cols {
            let nd = native_d.get(i, j);
            let md = model_d.get(i, j);
            if nd < threshold_sq {
                n_native += 1;
                if md < threshold_sq {
                    n_shared += 1;
                }
            }
            if md < threshold_sq {
                n_model += 1;
                if nd >= threshold_sq {
                    n_nonnat += 1;
                }
            }
        }
    }

    (n_shared, n_nonnat, n_native, n_model)
}

/// Kabsch / SVD superposition, replicating Biopython's `SVDSuperimposer` exactly.
///
/// Finds the rotation+translation putting `coords` onto `reference`:
/// ```text
/// av1 = mean(coords);  av2 = mean(reference)
/// a   = (coords-av1)^T @ (reference-av2)            # 3x3 correlation matrix
/// u, d, vt = svd(a)
/// rot = u @ vt
/// if det(rot) < 0: vt[2] = -vt[2]; rot = u @ vt     # fix reflection
/// tran = av2 - av1 @ rot
/// ```
/// Returns `(rot, tran)` with `rot` row-major: applied as
/// `transformed[k] = sum_j p[j]*rot[j][k] + tran[k]`.
pub fn kabsch(reference: &[[f32; 3]], coords: &[[f32; 3]]) -> ([[f32; 3]; 3], [f32; 3]) {
    let n = coords.len();
    debug_assert_eq!(n, reference.len(), "kabsch: coord count mismatch");
    let inv_n = 1.0f32 / n as f32;

    // Centroids (sum then divide, matching numpy `sum(coords)/n`).
    let mut av1 = Vector3::<f32>::zeros();
    let mut av2 = Vector3::<f32>::zeros();
    for k in 0..n {
        av1 += Vector3::new(coords[k][0], coords[k][1], coords[k][2]);
        av2 += Vector3::new(reference[k][0], reference[k][1], reference[k][2]);
    }
    av1 *= inv_n;
    av2 *= inv_n;

    // Correlation matrix a = (coords-av1)^T @ (reference-av2), a 3x3.
    // a[r][c] = sum_k (coords[k][r]-av1[r]) * (reference[k][c]-av2[c]).
    let mut a = Matrix3::<f32>::zeros();
    for k in 0..n {
        let p = Vector3::new(
            coords[k][0] - av1[0],
            coords[k][1] - av1[1],
            coords[k][2] - av1[2],
        );
        let q = Vector3::new(
            reference[k][0] - av2[0],
            reference[k][1] - av2[1],
            reference[k][2] - av2[2],
        );
        // outer product p (3x1) * q^T (1x3) accumulated.
        a += p * q.transpose();
    }

    // SVD: nalgebra gives a = U * Σ * Vᵀ, with `u` ≙ numpy u and `v_t` ≙ numpy vt.
    let svd = a.svd(true, true);
    let u = svd.u.expect("kabsch: SVD U not computed");
    let mut vt = svd.v_t.expect("kabsch: SVD Vᵀ not computed");

    // rot = u @ vt.
    let mut rot = u * vt;

    // Reflection fix: if det < 0, negate the 3rd ROW of vt and recompute.
    if rot.determinant() < 0.0 {
        // numpy: vt[2] = -vt[2]  (row index 2)
        vt[(2, 0)] = -vt[(2, 0)];
        vt[(2, 1)] = -vt[(2, 1)];
        vt[(2, 2)] = -vt[(2, 2)];
        rot = u * vt;
    }

    // tran = av2 - av1 @ rot.  (av1 @ rot)[c] = sum_r av1[r] * rot[r][c].
    let av1_rot = Vector3::new(
        av1[0] * rot[(0, 0)] + av1[1] * rot[(1, 0)] + av1[2] * rot[(2, 0)],
        av1[0] * rot[(0, 1)] + av1[1] * rot[(1, 1)] + av1[2] * rot[(2, 1)],
        av1[0] * rot[(0, 2)] + av1[1] * rot[(1, 2)] + av1[2] * rot[(2, 2)],
    );
    let tran = [av2[0] - av1_rot[0], av2[1] - av1_rot[1], av2[2] - av1_rot[2]];

    // rot row-major: rot_out[r][c] = rot[(r, c)].
    let rot_out = [
        [rot[(0, 0)], rot[(0, 1)], rot[(0, 2)]],
        [rot[(1, 0)], rot[(1, 1)], rot[(1, 2)]],
        [rot[(2, 0)], rot[(2, 1)], rot[(2, 2)]],
    ];

    (rot_out, tran)
}

/// Apply a Kabsch transform: `transformed[k][c] = sum_j coords[k][j]*rot[j][c] + tran[c]`.
/// (Matches numpy `coords @ rot + tran`.)
pub fn apply_transform(
    coords: &[[f32; 3]],
    rot: &[[f32; 3]; 3],
    tran: &[f32; 3],
) -> Vec<[f32; 3]> {
    coords
        .iter()
        .map(|p| {
            [
                p[0] * rot[0][0] + p[1] * rot[1][0] + p[2] * rot[2][0] + tran[0],
                p[0] * rot[0][1] + p[1] * rot[1][1] + p[2] * rot[2][1] + tran[1],
                p[0] * rot[0][2] + p[1] * rot[1][2] + p[2] * rot[2][2] + tran[2],
            ]
        })
        .collect()
}

/// RMSD = sqrt(sum|a-b|^2 / N) (Biopython's private `_rms`, no superposition).
pub fn rmsd(a: &[[f32; 3]], b: &[[f32; 3]]) -> f32 {
    debug_assert_eq!(a.len(), b.len(), "rmsd: coord count mismatch");
    let n = a.len();
    let mut sum = 0.0f32;
    for k in 0..n {
        let dx = a[k][0] - b[k][0];
        let dy = a[k][1] - b[k][1];
        let dz = a[k][2] - b[k][2];
        sum += dx * dx + dy * dy + dz * dz;
    }
    (sum / n as f32).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn dm(rows: usize, cols: usize, data: &[f32]) -> DistMatrix {
        DistMatrix {
            rows,
            cols,
            data: data.to_vec(),
        }
    }

    // --- residue_distances --------------------------------------------------
    // Expected values computed with the reference Cython kernel in the
    // .venv-baseline (DockQ.operations.residue_distances).

    #[test]
    fn residue_distances_small() {
        // 2 residues in chain1 (atoms_per_res = [2, 1]), 2 in chain2 ([1, 2]).
        // Flat coords (3 atoms in c1, 3 atoms in c2).
        let c1 = [
            [0.0f32, 0.0, 0.0],
            [1.0, 0.0, 0.0], // res0: 2 atoms
            [0.0, 5.0, 0.0], // res1: 1 atom
        ];
        let c2 = [
            [0.0f32, 0.5, 0.0], // res0: 1 atom
            [2.0, 0.0, 0.0],
            [0.0, 4.5, 0.0], // res1: 2 atoms
        ];
        let apr1 = [2usize, 1];
        let apr2 = [1usize, 2];
        let d = residue_distances(&c1, &c2, &apr1, &apr2);
        // Oracle output:
        // [[0.25, 1.0], [20.25, 0.25]]
        assert_relative_eq!(d.get(0, 0), 0.25, max_relative = 1e-6);
        assert_relative_eq!(d.get(0, 1), 1.0, max_relative = 1e-6);
        assert_relative_eq!(d.get(1, 0), 20.25, max_relative = 1e-6);
        assert_relative_eq!(d.get(1, 1), 0.25, max_relative = 1e-6);
    }

    #[test]
    fn residue_distances_early_break_quirk() {
        // First atom pair is far (>sqrt(1000) Å): the Cython early-break stores
        // that first pair's distance, NOT the true minimum.
        let c1 = [[0.0f32, 0.0, 0.0], [0.0, 0.0, 0.0]];
        let c2 = [[100.0f32, 0.0, 0.0], [0.0, 0.0, 0.5]];
        let apr = [2usize];
        let d = residue_distances(&c1, &c2, &apr, &apr);
        // Oracle: 10000.0 (NOT 0.25)
        assert_relative_eq!(d.get(0, 0), 10000.0, max_relative = 1e-6);

        // Same atoms, close pair first: true minimum is found.
        let c2b = [[0.0f32, 0.0, 0.5], [100.0, 0.0, 0.0]];
        let db = residue_distances(&c1, &c2b, &apr, &apr);
        assert_relative_eq!(db.get(0, 0), 0.25, max_relative = 1e-6);
    }

    // --- fnat_stats ---------------------------------------------------------

    #[test]
    fn fnat_stats_basic() {
        // 2x2 squared-distance matrices, threshold_sq = 25.0 (i.e. 5.0 Å).
        // native < 25: positions (0,0)=1, (0,1)=100, (1,0)=4, (1,1)=49
        //   -> native contacts at (0,0),(1,0)
        // model:        (0,0)=1, (0,1)=9,  (1,0)=400, (1,1)=16
        //   -> model contacts at (0,0),(0,1),(1,1)
        let native = dm(2, 2, &[1.0, 100.0, 4.0, 49.0]);
        let model = dm(2, 2, &[1.0, 9.0, 400.0, 16.0]);
        let (shared, nonnat, n_native, n_model) = fnat_stats(&model, &native, 25.0);
        // shared = model & native = {(0,0)} -> 1
        // n_native = 2 ; n_model = 3
        // nonnat = model & !native = {(0,1),(1,1)} -> 2
        assert_eq!((shared, nonnat, n_native, n_model), (1, 2, 2, 3));
    }

    // --- kabsch / apply_transform / rmsd ------------------------------------

    #[test]
    fn kabsch_90deg_about_z() {
        // reference = coords rotated +90° about z about the origin, then translated.
        // A +90° rotation about z maps (x,y,z) -> (-y, x, z).
        let coords = [
            [1.0f32, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [2.0, 3.0, 1.0],
            [-1.0, 2.0, 0.5],
        ];
        // reference[k] = R90z @ coords[k] + t, t = (5, -2, 1)
        let t = [5.0f32, -2.0, 1.0];
        let reference: Vec<[f32; 3]> = coords
            .iter()
            .map(|p| [-p[1] + t[0], p[0] + t[1], p[2] + t[2]])
            .collect();

        let (rot, tran) = kabsch(&reference, &coords);
        let transformed = apply_transform(&coords, &rot, &tran);

        // Superposition must be (near) exact.
        let r = rmsd(&transformed, &reference);
        assert!(r < 1e-4, "rms {r} too large");

        // rot should be the +90°-about-z rotation in the row-major,
        // right-multiply convention: p_out = p @ rot. For p=(1,0,0) -> (0,1,0):
        // rot row0 = (0,1,0); p=(0,1,0) -> (-1,0,0): rot row1 = (-1,0,0); rot[2]=(0,0,1).
        assert_relative_eq!(rot[0][0], 0.0, epsilon = 1e-4);
        assert_relative_eq!(rot[0][1], 1.0, epsilon = 1e-4);
        assert_relative_eq!(rot[1][0], -1.0, epsilon = 1e-4);
        assert_relative_eq!(rot[1][1], 0.0, epsilon = 1e-4);
        assert_relative_eq!(rot[2][2], 1.0, epsilon = 1e-4);
        assert_relative_eq!(tran[0], 5.0, epsilon = 1e-3);
        assert_relative_eq!(tran[1], -2.0, epsilon = 1e-3);
        assert_relative_eq!(tran[2], 1.0, epsilon = 1e-3);
    }

    #[test]
    fn kabsch_self_superposition_rms_zero() {
        let coords = [
            [1.5f32, -2.3, 0.7],
            [3.1, 0.0, -1.2],
            [-0.5, 4.4, 2.0],
            [2.2, 2.2, 2.2],
        ];
        let (rot, tran) = kabsch(&coords, &coords);
        let transformed = apply_transform(&coords, &rot, &tran);
        let r = rmsd(&transformed, &coords);
        assert!(r < 1e-5, "self-superposition rms {r} should be ~0");
    }

    #[test]
    fn rmsd_matches_known() {
        // Two point sets; sum|a-b|^2 = (1+1+4) over 3 points.
        let a = [[0.0f32, 0.0, 0.0], [1.0, 1.0, 1.0], [2.0, 0.0, 0.0]];
        let b = [[1.0f32, 0.0, 0.0], [1.0, 1.0, 0.0], [2.0, 0.0, 2.0]];
        // diffs^2: (1)+(1)+(4) = 6 ; /3 = 2 ; sqrt = 1.4142135
        let r = rmsd(&a, &b);
        assert_relative_eq!(r, (2.0f32).sqrt(), max_relative = 1e-6);
    }

    #[test]
    fn kabsch_real_data_sanity() {
        // Real backbone-like coords (model) vs a known rigid-body image (reference),
        // verifying the full kabsch->apply->rmsd pipeline recovers ~0 rms.
        let model = [
            [12.501f32, 39.048, 28.539],
            [15.552, 39.41, 30.835],
            [17.984, 37.36, 28.738],
            [16.412, 34.13, 27.97],
            [16.06, 33.69, 31.74],
        ];
        // R = rotation by 120° about (1,1,1)/sqrt(3) maps (x,y,z)->(z,x,y).
        let tr = [10.0f32, -5.0, 3.0];
        let reference: Vec<[f32; 3]> = model
            .iter()
            .map(|p| [p[2] + tr[0], p[0] + tr[1], p[1] + tr[2]])
            .collect();
        let (rot, tran) = kabsch(&reference, &model);
        let transformed = apply_transform(&model, &rot, &tran);
        let r = rmsd(&transformed, &reference);
        assert!(r < 1e-3, "real-data sanity rms {r} should be ~0");
    }
}
