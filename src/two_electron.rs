//! Two-electron repulsion integrals (standard "slow" implementation).
//!
//! Computes electron-electron repulsion integrals using the standard
//! recursive method. This is the reference implementation - slower but
//! clearer algorithm.
//!
//! Reference: PyQuante2 - pyquante2/ints/two.py

use crate::common::{CGBF, Vec3};
use crate::utils::{binomial_prefactor, fact_ratio2, fgamma, gaussian_normalization, gaussian_product_center};
use std::f64::consts::PI;

/// Compute two-electron repulsion integral between four contracted Gaussian basis functions
///
/// The electron repulsion integral is: (ij|kl) = ⟨φ_i φ_j|1/r₁₂|φ_k φ_l⟩
///
/// This represents the Coulomb repulsion energy between two electron densities:
/// one formed by basis functions i and j, the other by k and l.
pub fn electron_repulsion(bra1: &CGBF, bra2: &CGBF, ket1: &CGBF, ket2: &CGBF) -> f64 {
    let mut sum = 0.0;

    // Compute normalization constants for the contracted basis functions
    let (la, ma, na) = bra1.shell;
    let (lb, mb, nb) = bra2.shell;
    let (lc, mc, nc) = ket1.shell;
    let (ld, md, nd) = ket2.shell;

    for p1 in &bra1.primitives {
        for p2 in &bra2.primitives {
            for p3 in &ket1.primitives {
                for p4 in &ket2.primitives {
                    let norm_a = gaussian_normalization(p1.exponent, la, ma, na);
                    let norm_b = gaussian_normalization(p2.exponent, lb, mb, nb);
                    let norm_c = gaussian_normalization(p3.exponent, lc, mc, nc);
                    let norm_d = gaussian_normalization(p4.exponent, ld, md, nd);

                    sum += p1.coefficient
                        * p2.coefficient
                        * p3.coefficient
                        * p4.coefficient
                        * norm_a
                        * norm_b
                        * norm_c
                        * norm_d
                        * coulomb_repulsion_primitive(
                            bra1.origin,
                            bra1.shell,
                            p1.exponent,
                            bra2.origin,
                            bra2.shell,
                            p2.exponent,
                            ket1.origin,
                            ket1.shell,
                            p3.exponent,
                            ket2.origin,
                            ket2.shell,
                            p4.exponent,
                        );
                }
            }
        }
    }

    sum
}

/// Compute Coulomb repulsion between four primitive Gaussian basis functions
///
/// This is the core function that computes the electron repulsion integral
/// between four primitive Gaussians using the THO (Taketa-Huzinaga-O-ohata) method.
#[allow(clippy::too_many_arguments)]
fn coulomb_repulsion_primitive(
    a: Vec3,
    lmn_a: (i32, i32, i32),
    alpha_a: f64,
    b: Vec3,
    lmn_b: (i32, i32, i32),
    alpha_b: f64,
    c: Vec3,
    lmn_c: (i32, i32, i32),
    alpha_c: f64,
    d: Vec3,
    lmn_d: (i32, i32, i32),
    alpha_d: f64,
) -> f64 {
    let (la, ma, na) = lmn_a;
    let (lb, mb, nb) = lmn_b;
    let (lc, mc, nc) = lmn_c;
    let (ld, md, nd) = lmn_d;

    let rab2 = a.distance_squared(&b);
    let rcd2 = c.distance_squared(&d);

    // Gaussian product centers
    let px = gaussian_product_center(alpha_a, a.x, alpha_b, b.x);
    let py = gaussian_product_center(alpha_a, a.y, alpha_b, b.y);
    let pz = gaussian_product_center(alpha_a, a.z, alpha_b, b.z);
    let p = Vec3::new(px, py, pz);

    let qx = gaussian_product_center(alpha_c, c.x, alpha_d, d.x);
    let qy = gaussian_product_center(alpha_c, c.y, alpha_d, d.y);
    let qz = gaussian_product_center(alpha_c, c.z, alpha_d, d.z);
    let q = Vec3::new(qx, qy, qz);

    let rpq2 = p.distance_squared(&q);

    let gamma1 = alpha_a + alpha_b;
    let gamma2 = alpha_c + alpha_d;
    let delta = 0.25 * (1.0 / gamma1 + 1.0 / gamma2);

    // Compute B arrays for each dimension
    let bx = b_array(la, lb, lc, ld, px, a.x, b.x, qx, c.x, d.x, gamma1, gamma2, delta);
    let by = b_array(ma, mb, mc, md, py, a.y, b.y, qy, c.y, d.y, gamma1, gamma2, delta);
    let bz = b_array(na, nb, nc, nd, pz, a.z, b.z, qz, c.z, d.z, gamma1, gamma2, delta);

    // Sum over all terms
    let mut sum = 0.0;
    for i in 0..=(la + lb + lc + ld) as usize {
        for j in 0..=(ma + mb + mc + md) as usize {
            for k in 0..=(na + nb + nc + nd) as usize {
                sum += bx[i] * by[j] * bz[k] * fgamma((i + j + k) as i32, 0.25 * rpq2 / delta);
            }
        }
    }

    // Final prefactor and exponentials
    let prefactor = 2.0 * PI.powf(2.5) / (gamma1 * gamma2 * (gamma1 + gamma2).sqrt());
    let exp_term = (-alpha_a * alpha_b * rab2 / gamma1).exp()
        * (-alpha_c * alpha_d * rcd2 / gamma2).exp();

    prefactor * exp_term * sum
}

/// Compute the B array for two-electron integrals
///
/// This implements the recursion relation from THO eq. 2.22
#[allow(clippy::too_many_arguments)]
fn b_array(
    l1: i32,
    l2: i32,
    l3: i32,
    l4: i32,
    p: f64,
    a: f64,
    b: f64,
    q: f64,
    c: f64,
    d: f64,
    g1: f64,
    g2: f64,
    delta: f64,
) -> Vec<f64> {
    let imax = (l1 + l2 + l3 + l4 + 1) as usize;
    let mut result = vec![0.0; imax];

    for i1 in 0..=(l1 + l2) {
        for i2 in 0..=(l3 + l4) {
            for r1 in 0..=(i1 / 2) {
                for r2 in 0..=(i2 / 2) {
                    for u in 0..=((i1 + i2) / 2 - r1 - r2) {
                        let idx = (i1 + i2 - 2 * (r1 + r2) - u) as usize;
                        result[idx] += b_term(i1, i2, r1, r2, u, l1, l2, l3, l4, p, a, b, q, c, d, g1, g2, delta);
                    }
                }
            }
        }
    }

    result
}

/// Compute a single B term (THO eq. 2.22)
#[allow(clippy::too_many_arguments)]
fn b_term(
    i1: i32,
    i2: i32,
    r1: i32,
    r2: i32,
    u: i32,
    l1: i32,
    l2: i32,
    l3: i32,
    l4: i32,
    px: f64,
    ax: f64,
    bx: f64,
    qx: f64,
    cx: f64,
    dx: f64,
    gamma1: f64,
    gamma2: f64,
    delta: f64,
) -> f64 {
    f_b(i1, l1, l2, px, ax, bx, r1, gamma1)
        * (-1.0_f64).powi(i2)
        * f_b(i2, l3, l4, qx, cx, dx, r2, gamma2)
        * (-1.0_f64).powi(u)
        * fact_ratio2(i1 + i2 - 2 * (r1 + r2), u)
        * (qx - px).powi(i1 + i2 - 2 * (r1 + r2) - 2 * u)
        / delta.powi(i1 + i2 - 2 * (r1 + r2) - u)
}

/// Helper function f_B for B term calculation
fn f_b(i: i32, l1: i32, l2: i32, p: f64, a: f64, b: f64, r: i32, g: f64) -> f64 {
    binomial_prefactor(i, l1, l2, p - a, p - b) * b0(i, r, g)
}

/// Helper function B_0 for f_B calculation
fn b0(i: i32, r: i32, g: f64) -> f64 {
    fact_ratio2(i, r) * (4.0 * g).powi(r - i)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Primitive;
    use approx::assert_relative_eq;

    /// Create a simple s-orbital at given position
    fn s_orbital(exponent: f64, origin: Vec3) -> CGBF {
        CGBF {
            origin,
            shell: (0, 0, 0),
            primitives: vec![Primitive {
                exponent,
                coefficient: 1.0,
            }],
        }
    }

    #[test]
    fn test_coulomb_repulsion_primitive_basic() {
        // From PyQuante2 doctest:
        // coulomb_repulsion(p1,1.,lmn,1.,p1,1.,lmn,1.,p1,1.,lmn,1.,p1,1.,lmn,1.) = 4.373355
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let lmn = (0, 0, 0);
        let result = coulomb_repulsion_primitive(
            origin, lmn, 1.0, origin, lmn, 1.0, origin, lmn, 1.0, origin, lmn, 1.0,
        );
        assert_relative_eq!(result, 4.373355, epsilon = 1e-5);
    }

    #[test]
    fn test_coulomb_repulsion_primitive_separated() {
        // From PyQuante2 doctest:
        // coulomb_repulsion(p1,1.,lmn,1.,p1,1.,lmn,1.,p2,1.,lmn,1.,p2,1.,lmn,1.) = 3.266127
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 1.0);
        let lmn = (0, 0, 0);
        let result = coulomb_repulsion_primitive(p1, lmn, 1.0, p1, lmn, 1.0, p2, lmn, 1.0, p2, lmn, 1.0);
        assert_relative_eq!(result, 3.266127, epsilon = 1e-5);
    }

    #[test]
    fn test_coulomb_repulsion_primitive_mixed() {
        // From PyQuante2 doctest:
        // coulomb_repulsion(p1,1.,lmn,1.,p2,1.,lmn,1.,p1,1.,lmn,1.,p2,1.,lmn,1.) = 1.6088672
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 1.0);
        let lmn = (0, 0, 0);
        let result = coulomb_repulsion_primitive(p1, lmn, 1.0, p2, lmn, 1.0, p1, lmn, 1.0, p2, lmn, 1.0);
        assert_relative_eq!(result, 1.6088672, epsilon = 1e-5);
    }

    #[test]
    fn test_eri_s_orbitals_same() {
        // From PyQuante2 doctest: ERI(s,s,s,s) = 1.128379
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s = s_orbital(1.0, origin);
        let result = electron_repulsion(&s, &s, &s, &s);
        assert_relative_eq!(result, 1.128379, epsilon = 1e-5);
    }

    #[test]
    fn test_eri_s_and_s_separated() {
        // From PyQuante2 doctest:
        // s2 = cgbf((0,0,1),(0,0,0),[1],[1])
        // This is an s-orbital at position z=1, NOT a pz orbital!
        // ERI(s,s,s2,s2) = 0.842701
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let pos_z1 = Vec3::new(0.0, 0.0, 1.0);
        let s = s_orbital(1.0, origin);
        let s2 = s_orbital(1.0, pos_z1);  // s-orbital at z=1
        let result = electron_repulsion(&s, &s, &s2, &s2);
        assert_relative_eq!(result, 0.842701, epsilon = 1e-5);
    }

    #[test]
    fn test_eri_symmetry_zero() {
        // From PyQuante2 doctest: ERI(s,s,s,px) = 0 (by symmetry)
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s = s_orbital(1.0, origin);
        let px = CGBF {
            origin,
            shell: (1, 0, 0), // px orbital has l=1
            primitives: vec![Primitive {
                exponent: 1.0,
                coefficient: 1.0,
            }],
        };
        let result = electron_repulsion(&s, &s, &s, &px);
        assert!(result.abs() < 1e-10); // Should be essentially zero
    }

    #[test]
    fn test_eri_px_orbitals() {
        // From PyQuante2 doctest: ERI(s,s,px,px) = 0.940315972579
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s = s_orbital(1.0, origin);
        let px = CGBF {
            origin,
            shell: (1, 0, 0),
            primitives: vec![Primitive {
                exponent: 1.0,
                coefficient: 1.0,
            }],
        };
        let result = electron_repulsion(&s, &s, &px, &px);
        assert_relative_eq!(result, 0.940315972579, epsilon = 1e-5);
    }
}

