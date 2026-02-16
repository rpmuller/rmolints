//! Rys quadrature for n ≥ 6 using the Golub-Welsch algorithm.
//!
//! For n ≤ 5, hard-coded GAMESS polynomial tables (Root123/4/5) give high
//! accuracy over piecewise x intervals.  For n ≥ 6 those tables were never
//! written in the reference code, so we use the standard numerical procedure:
//!
//! 1. Compute 2n moments  μ_k = F_k(x)  via the Boys (fgamma) function.
//! 2. Run the modified Chebyshev algorithm (Wheeler 1974 / Gautschi 2004) to
//!    convert the moments into the n×n symmetric tridiagonal Jacobi matrix.
//! 3. Diagonalise that matrix with the TQLI algorithm (Numerical Recipes).
//! 4. Return  t_i = √λ_i  as Rys roots and  β₀ · v_i[0]²  as weights.
//!
//! References:
//!   Gautschi, "Orthogonal Polynomials: Computation and Approximation" (2004),
//!     Algorithm 1.5 (CHEB).
//!   Press et al., "Numerical Recipes in C" (1992), §11.3 (tqli).
//!   Golub & Welsch, Math. Comp. 23 (1969) 221–230.

use crate::utils::fgamma;

/// Compute `n` Rys roots and weights for quadrature argument `x`.
///
/// Valid for any n ≥ 1; in practice only called for n ≥ 6 since smaller
/// orders use the faster polynomial tables.
pub(crate) fn rys_roots_general(n: usize, x: f64, roots: &mut [f64], weights: &mut [f64]) {
    // Step 1: moments μ_k = F_k(x)
    let mut moments = vec![0.0f64; 2 * n];
    for k in 0..2 * n {
        moments[k] = fgamma(k as i32, x);
    }

    // Step 2: modified Chebyshev algorithm → Jacobi matrix coefficients
    let (alpha, beta) = modified_chebyshev(&moments, n);

    // Step 3: symmetric tridiagonal eigendecomposition
    // d = diagonal = alpha; e[i] = off-diagonal coupling d[i] and d[i+1]
    let mut d = alpha.clone();
    let mut e = vec![0.0f64; n];
    for i in 0..n - 1 {
        e[i] = beta[i + 1].max(0.0).sqrt();
    }
    // e[n-1] = 0.0 already

    let mut z = vec![vec![0.0f64; n]; n];
    for i in 0..n {
        z[i][i] = 1.0;
    }

    tqli(&mut d, &mut e, n, &mut z);

    // Step 4: roots  t_i = √λ_i,  weights  W_i = β₀ · z[0][i]²
    for i in 0..n {
        roots[i] = d[i].max(0.0).sqrt();
        weights[i] = beta[0] * z[0][i] * z[0][i];
    }

    sort_rys(n, roots, weights);
}

// ─────────────────────────────────────────────────────────────────────────────
// Modified Chebyshev algorithm
// ─────────────────────────────────────────────────────────────────────────────

/// Convert 2n power moments to n Jacobi matrix coefficients.
///
/// Three-term recurrence (Gautschi 2004, Algorithm 1.5 "CHEB"):
///   σ_{k,l} = σ_{k-1,l+1} − α_{k-1}·σ_{k-1,l} − β_{k-1}·σ_{k-2,l}
///   α_k = σ_{k,k+1}/σ_{k,k} − σ_{k-1,k}/σ_{k-1,k-1}
///   β_k = σ_{k,k}/σ_{k-1,k-1}
///   with  α_0 = μ₁/μ₀,  β_0 = μ₀,  σ_{-1,·} = 0.
fn modified_chebyshev(moments: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let two_n = 2 * n;
    let mut alpha = vec![0.0f64; n];
    let mut beta = vec![0.0f64; n];

    // Rolling σ rows
    let mut sig_m2 = vec![0.0f64; two_n + 1]; // σ_{k-2}, initially σ_{-1} = 0
    let mut sig_m1 = moments.to_vec();          // σ_{k-1}, initially σ_{0} = μ
    sig_m1.resize(two_n + 1, 0.0);
    let mut sig = vec![0.0f64; two_n + 1];

    beta[0] = moments[0];
    alpha[0] = moments[1] / moments[0];

    for k in 1..n {
        for l in k..(two_n - k) {
            sig[l] = sig_m1[l + 1]
                - alpha[k - 1] * sig_m1[l]
                - beta[k - 1] * sig_m2[l];
        }
        alpha[k] = sig[k + 1] / sig[k] - sig_m1[k] / sig_m1[k - 1];
        beta[k] = sig[k] / sig_m1[k - 1];

        std::mem::swap(&mut sig_m2, &mut sig_m1);
        std::mem::swap(&mut sig_m1, &mut sig);
        sig.iter_mut().for_each(|v| *v = 0.0);
    }

    (alpha, beta)
}

// ─────────────────────────────────────────────────────────────────────────────
// TQLI – symmetric tridiagonal eigenvalue / eigenvector solver
// ─────────────────────────────────────────────────────────────────────────────

/// Symmetric tridiagonal eigendecomposition (Numerical Recipes §11.3, tqli).
///
/// * `d[0..n-1]` – diagonal on entry, eigenvalues on exit.
/// * `e[0..n-1]` – off-diagonal on entry: `e[i]` couples `d[i]` and `d[i+1]`;
///   `e[n-1]` is unused and set to zero.  Destroyed on exit.
/// * `z[n][n]`   – identity on entry; `z[k][i]` = k-th component of
///   i-th eigenvector on exit.
fn tqli(d: &mut [f64], e: &mut [f64], n: usize, z: &mut Vec<Vec<f64>>) {
    const MAX_ITER: usize = 30;

    e[n - 1] = 0.0;

    for l in 0..n {
        let mut iter = 0usize;

        loop {
            // Locate lowest negligible off-diagonal element
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if (e[m].abs() + dd) == dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break; // d[l] has converged
            }

            assert!(iter < MAX_ITER, "TQLI: too many iterations (l={l}, m={m})");
            iter += 1;

            // Wilkinson shift
            let sg = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let sr = (sg * sg + 1.0).sqrt();
            let mut g = d[m] - d[l] + e[l] / (sg + sr.copysign(sg));

            let mut s = 1.0f64;
            let mut c = 1.0f64;
            let mut p = 0.0f64;
            let mut goto_done = false;

            // QL sweep
            let mut i = m;
            while i > l {
                i -= 1;

                let f = s * e[i];
                let b = c * e[i];
                let r1 = (f * f + g * g).sqrt();
                e[i + 1] = r1;

                if r1 == 0.0 {
                    // Isolated eigenvalue — skip the d[l]/e[l] update
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    goto_done = true;
                    break;
                }

                s = f / r1;
                c = g / r1;
                g = d[i + 1] - p;
                let r2 = (d[i] - g) * s + 2.0 * c * b;
                p = s * r2;
                d[i + 1] = g + p;
                g = c * r2 - b;

                // Accumulate plane rotation into eigenvector matrix
                for k in 0..n {
                    let fk = z[k][i + 1];
                    z[k][i + 1] = s * z[k][i] + c * fk;
                    z[k][i] = c * z[k][i] - s * fk;
                }
            }

            if !goto_done {
                d[l] -= p;
                e[l] = g;
                e[m] = 0.0;
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Sort Rys roots and weights in ascending order of root value.
/// Insertion sort is fine for n ≤ 8.
fn sort_rys(n: usize, roots: &mut [f64], weights: &mut [f64]) {
    for i in 1..n {
        let r = roots[i];
        let w = weights[i];
        let mut j = i;
        while j > 0 && roots[j - 1] > r {
            roots[j] = roots[j - 1];
            weights[j] = weights[j - 1];
            j -= 1;
        }
        roots[j] = r;
        weights[j] = w;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Verify that the n-point Rys quadrature rule correctly integrates
    /// all moment integrands F_k(x) = ∫₀¹ t^{2k} exp(-xt²) dt for k = 0..2n-1.
    fn check_moments(n: usize, x: f64) {
        let mut roots = vec![0.0f64; n];
        let mut weights = vec![0.0f64; n];
        rys_roots_general(n, x, &mut roots, &mut weights);

        for k in 0..2 * n {
            let exact = fgamma(k as i32, x);
            let approx: f64 = weights
                .iter()
                .zip(roots.iter())
                .map(|(&w, &t)| w * t.powi(2 * k as i32))
                .sum();
            assert!(
                (approx - exact).abs() < 1e-8 * exact.abs().max(1e-14),
                "n={n}, x={x}, k={k}: quadrature={approx:.6e}, exact={exact:.6e}"
            );
        }
    }

    #[test]
    fn test_n6_moments_small_x() {
        check_moments(6, 0.1);
    }

    #[test]
    fn test_n6_moments_moderate_x() {
        check_moments(6, 5.0);
    }

    #[test]
    fn test_n6_moments_large_x() {
        check_moments(6, 20.0);
    }

    #[test]
    fn test_n7_moments() {
        check_moments(7, 3.0);
    }

    #[test]
    fn test_roots_in_unit_interval() {
        // n ≤ 7 is well-conditioned (f orbitals require at most n = 7).
        // n = 8 hits the ill-conditioning limit of the power-moment Chebyshev
        // algorithm and is not tested here.
        for &n in &[6usize, 7] {
            for &x in &[0.01, 1.0, 10.0] {
                let mut roots = vec![0.0f64; n];
                let mut weights = vec![0.0f64; n];
                rys_roots_general(n, x, &mut roots, &mut weights);
                for &t in roots.iter().take(n) {
                    assert!(t > 0.0 && t < 1.0,
                        "n={n}, x={x}: root {t} not in (0,1)");
                }
                for &w in weights.iter().take(n) {
                    assert!(w > 0.0, "n={n}, x={x}: weight {w} not positive");
                }
            }
        }
    }
}
