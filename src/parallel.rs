//! Parallel computation of molecular integrals using Rayon
//!
//! This module provides multithreaded versions of integral computations,
//! particularly useful for computing batches of two-electron integrals.

use crate::common::CGBF;
use crate::{hgp, rys, two_electron};
use rayon::prelude::*;

/// Method to use for computing two-electron integrals
#[derive(Debug, Clone, Copy)]
pub enum ERIMethod {
    /// Standard Taketa-Huzinaga-O-ohata method
    Standard,
    /// Rys quadrature method
    Rys,
    /// Head-Gordon-Pople method
    HeadGordonPople,
}

/// Compute a batch of two-electron integrals in parallel
///
/// # Arguments
/// * `basis` - Slice of contracted Gaussian basis functions
/// * `indices` - List of (i,j,k,l) quartet indices to compute
/// * `method` - Which ERI method to use
///
/// # Returns
/// Vector of integral values in the same order as `indices`
///
/// # Example
/// ```
/// use rmolints::parallel::{compute_eris_parallel, ERIMethod};
/// use rmolints::common::*;
///
/// let s = CGBF {
///     origin: Vec3::new(0.0, 0.0, 0.0),
///     shell: (0, 0, 0),
///     primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
/// };
/// let basis = vec![s.clone(), s.clone()];
/// let indices = vec![(0, 0, 0, 0), (0, 1, 0, 1)];
///
/// let results = compute_eris_parallel(&basis, &indices, ERIMethod::Standard);
/// assert_eq!(results.len(), 2);
/// ```
pub fn compute_eris_parallel(
    basis: &[CGBF],
    indices: &[(usize, usize, usize, usize)],
    method: ERIMethod,
) -> Vec<f64> {
    indices
        .par_iter()
        .map(|&(i, j, k, l)| {
            let bra1 = &basis[i];
            let bra2 = &basis[j];
            let ket1 = &basis[k];
            let ket2 = &basis[l];

            match method {
                ERIMethod::Standard => two_electron::electron_repulsion(bra1, bra2, ket1, ket2),
                ERIMethod::Rys => rys::electron_repulsion_rys(bra1, bra2, ket1, ket2),
                ERIMethod::HeadGordonPople => hgp::electron_repulsion_hgp(bra1, bra2, ket1, ket2),
            }
        })
        .collect()
}

/// Compute all unique two-electron integrals for a basis set in parallel
///
/// This computes the full ERI tensor exploiting 8-fold symmetry:
/// (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk) = (kl|ij) = (lk|ij) = (kl|ji) = (lk|ji)
///
/// # Arguments
/// * `basis` - Slice of contracted Gaussian basis functions
/// * `method` - Which ERI method to use
///
/// # Returns
/// Vector of (i, j, k, l, value) tuples for all unique integrals where i >= j, k >= l, ij >= kl
///
/// # Example
/// ```
/// use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
/// use rmolints::common::*;
///
/// let s1 = CGBF {
///     origin: Vec3::new(0.0, 0.0, 0.0),
///     shell: (0, 0, 0),
///     primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
/// };
/// let s2 = CGBF {
///     origin: Vec3::new(0.0, 0.0, 1.0),
///     shell: (0, 0, 0),
///     primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
/// };
/// let basis = vec![s1, s2];
///
/// let eris = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);
/// // For 2 basis functions with 8-fold symmetry, we get 6 unique integrals
/// assert_eq!(eris.len(), 6);
/// ```
pub fn compute_eri_tensor_parallel(
    basis: &[CGBF],
    method: ERIMethod,
) -> Vec<(usize, usize, usize, usize, f64)> {
    let n = basis.len();

    // Generate all unique index quartets exploiting 8-fold symmetry
    let mut indices = Vec::new();
    for i in 0..n {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..=i {
                let l_max = if k == i { j } else { k };
                for l in 0..=l_max {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        indices.push((i, j, k, l));
                    }
                }
            }
        }
    }

    // Compute integrals in parallel
    indices
        .par_iter()
        .map(|&(i, j, k, l)| {
            let bra1 = &basis[i];
            let bra2 = &basis[j];
            let ket1 = &basis[k];
            let ket2 = &basis[l];

            let value = match method {
                ERIMethod::Standard => two_electron::electron_repulsion(bra1, bra2, ket1, ket2),
                ERIMethod::Rys => rys::electron_repulsion_rys(bra1, bra2, ket1, ket2),
                ERIMethod::HeadGordonPople => hgp::electron_repulsion_hgp(bra1, bra2, ket1, ket2),
            };

            (i, j, k, l, value)
        })
        .collect()
}

/// Compute a full (non-symmetric) ERI tensor in parallel
///
/// This computes all N^4 integrals without exploiting symmetry.
/// Useful when you need direct indexing [i][j][k][l].
///
/// # Arguments
/// * `basis` - Slice of contracted Gaussian basis functions
/// * `method` - Which ERI method to use
///
/// # Returns
/// 4D vector indexed as [i][j][k][l]
///
/// # Warning
/// This computes N^4 integrals which can be very large!
/// For N=10, this is 10,000 integrals. For N=100, it's 100 million integrals.
pub fn compute_eri_tensor_full_parallel(
    basis: &[CGBF],
    method: ERIMethod,
) -> Vec<Vec<Vec<Vec<f64>>>> {
    let n = basis.len();

    // Generate all index quartets
    let indices: Vec<(usize, usize, usize, usize)> = (0..n)
        .flat_map(|i| {
            (0..n).flat_map(move |j| (0..n).flat_map(move |k| (0..n).map(move |l| (i, j, k, l))))
        })
        .collect();

    // Compute all integrals in parallel
    let values = compute_eris_parallel(basis, &indices, method);

    // Reshape into 4D array
    let mut tensor = vec![vec![vec![vec![0.0; n]; n]; n]; n];
    for (idx, &(i, j, k, l)) in indices.iter().enumerate() {
        tensor[i][j][k][l] = values[idx];
    }

    tensor
}

/// Compute two-electron integrals for specified basis function pairs in parallel
///
/// More flexible than `compute_eri_tensor_parallel` - you specify exactly which
/// integral shells you want computed.
///
/// # Arguments
/// * `bra1_indices` - Indices of first basis functions
/// * `bra2_indices` - Indices of second basis functions
/// * `ket1_indices` - Indices of third basis functions
/// * `ket2_indices` - Indices of fourth basis functions
/// * `basis` - Slice of contracted Gaussian basis functions
/// * `method` - Which ERI method to use
///
/// # Returns
/// Vector of integral values, one per set of indices
pub fn compute_eri_shells_parallel(
    quartets: &[(usize, usize, usize, usize)],
    basis: &[CGBF],
    method: ERIMethod,
) -> Vec<f64> {
    compute_eris_parallel(basis, quartets, method)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::*;
    use approx::assert_relative_eq;

    fn make_s_orbital(alpha: f64, pos: Vec3) -> CGBF {
        CGBF {
            origin: pos,
            shell: (0, 0, 0),
            primitives: vec![Primitive {
                exponent: alpha,
                coefficient: 1.0,
            }],
        }
    }

    #[test]
    fn test_parallel_single_integral() {
        let s = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let basis = vec![s];
        let indices = vec![(0, 0, 0, 0)];

        let results = compute_eris_parallel(&basis, &indices, ERIMethod::Standard);
        assert_eq!(results.len(), 1);
        assert_relative_eq!(results[0], 1.128379, epsilon = 1e-5);
    }

    #[test]
    fn test_parallel_multiple_integrals() {
        let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 1.0));
        let basis = vec![s1, s2];
        let indices = vec![(0, 0, 0, 0), (0, 1, 0, 1), (1, 1, 1, 1)];

        let results = compute_eris_parallel(&basis, &indices, ERIMethod::Standard);
        assert_eq!(results.len(), 3);
        assert_relative_eq!(results[0], 1.128379, epsilon = 1e-5);
        assert_relative_eq!(results[2], 1.128379, epsilon = 1e-5);
    }

    #[test]
    fn test_parallel_vs_serial_standard() {
        let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 1.0));
        let basis = vec![s1.clone(), s2.clone()];

        // Serial computation
        let serial = two_electron::electron_repulsion(&s1, &s2, &s1, &s2);

        // Parallel computation
        let indices = vec![(0, 1, 0, 1)];
        let parallel = compute_eris_parallel(&basis, &indices, ERIMethod::Standard);

        assert_relative_eq!(serial, parallel[0], epsilon = 1e-10);
    }

    #[test]
    fn test_parallel_vs_serial_rys() {
        let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 1.0));
        let basis = vec![s1.clone(), s2.clone()];

        let serial = rys::electron_repulsion_rys(&s1, &s2, &s1, &s2);
        let indices = vec![(0, 1, 0, 1)];
        let parallel = compute_eris_parallel(&basis, &indices, ERIMethod::Rys);

        assert_relative_eq!(serial, parallel[0], epsilon = 1e-10);
    }

    #[test]
    fn test_parallel_vs_serial_hgp() {
        let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 1.0));
        let basis = vec![s1.clone(), s2.clone()];

        let serial = hgp::electron_repulsion_hgp(&s1, &s2, &s1, &s2);
        let indices = vec![(0, 1, 0, 1)];
        let parallel = compute_eris_parallel(&basis, &indices, ERIMethod::HeadGordonPople);

        assert_relative_eq!(serial, parallel[0], epsilon = 1e-10);
    }

    #[test]
    fn test_eri_tensor_symmetry() {
        let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 1.0));
        let basis = vec![s1, s2];

        let eris = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);

        // For 2 basis functions, should have 6 unique integrals with 8-fold symmetry
        // (0,0,0,0), (1,0,0,0), (1,0,1,0), (1,1,0,0), (1,1,1,0), (1,1,1,1)
        assert_eq!(eris.len(), 6);

        // Check specific values
        for (i, j, k, l, value) in &eris {
            if (*i, *j, *k, *l) == (0, 0, 0, 0) {
                assert_relative_eq!(*value, 1.128379, epsilon = 1e-5);
            } else if (*i, *j, *k, *l) == (1, 1, 1, 1) {
                assert_relative_eq!(*value, 1.128379, epsilon = 1e-5);
            }
        }
    }

    #[test]
    fn test_all_methods_agree() {
        let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
        let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.5));
        let basis = vec![s1, s2];
        let indices = vec![(0, 1, 0, 1)];

        let std_result = compute_eris_parallel(&basis, &indices, ERIMethod::Standard);
        let rys_result = compute_eris_parallel(&basis, &indices, ERIMethod::Rys);
        let hgp_result = compute_eris_parallel(&basis, &indices, ERIMethod::HeadGordonPople);

        assert_relative_eq!(std_result[0], rys_result[0], epsilon = 1e-5);
        assert_relative_eq!(std_result[0], hgp_result[0], epsilon = 1e-5);
    }
}
