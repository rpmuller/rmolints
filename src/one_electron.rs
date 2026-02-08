//! One-electron integrals: overlap, kinetic energy, and nuclear attraction.
//!
//! These integrals are computed between Gaussian basis functions and are
//! fundamental building blocks for quantum chemistry calculations.
//!
//! Reference: PyQuante2 - pyquante2/ints/one.py

use crate::common::{CGBF, Vec3};
use crate::utils::{binomial_prefactor, fact2, fgamma, gaussian_product_center};
use std::f64::consts::PI;

/// Compute overlap integral between two contracted Gaussian basis functions
///
/// The overlap integral is: S_ij = ⟨φ_i|φ_j⟩
///
/// This is the high-level interface that handles contracted basis functions.
/// The normalization constants are NOT included here - they should be
/// part of the basis function definition.
pub fn overlap(bra: &CGBF, ket: &CGBF) -> f64 {
    let mut sum = 0.0;
    for p1 in &bra.primitives {
        for p2 in &ket.primitives {
            sum += p1.coefficient
                * p2.coefficient
                * overlap_primitive(
                    p1.exponent,
                    bra.shell,
                    bra.origin,
                    p2.exponent,
                    ket.shell,
                    ket.origin,
                );
        }
    }
    sum
}

/// Compute overlap integral between two primitive Gaussian basis functions
///
/// Full form of the overlap integral from THO eq. 2.12
fn overlap_primitive(
    alpha1: f64,
    lmn1: (i32, i32, i32),
    a: Vec3,
    alpha2: f64,
    lmn2: (i32, i32, i32),
    b: Vec3,
) -> f64 {
    let (l1, m1, n1) = lmn1;
    let (l2, m2, n2) = lmn2;

    let rab2 = a.distance_squared(&b);
    let gamma = alpha1 + alpha2;

    let px = gaussian_product_center(alpha1, a.x, alpha2, b.x);
    let py = gaussian_product_center(alpha1, a.y, alpha2, b.y);
    let pz = gaussian_product_center(alpha1, a.z, alpha2, b.z);

    let pre = (PI / gamma).powf(1.5) * (-alpha1 * alpha2 * rab2 / gamma).exp();

    let wx = overlap_1d(l1, l2, px - a.x, px - b.x, gamma);
    let wy = overlap_1d(m1, m2, py - a.y, py - b.y, gamma);
    let wz = overlap_1d(n1, n2, pz - a.z, pz - b.z, gamma);

    pre * wx * wy * wz
}

/// One-dimensional component of the overlap integral (THO eq. 2.12)
fn overlap_1d(l1: i32, l2: i32, pax: f64, pbx: f64, gamma: f64) -> f64 {
    let mut total = 0.0;
    let imax = 1 + ((l1 + l2) as f64 / 2.0).floor() as i32;

    for i in 0..imax {
        total += binomial_prefactor(2 * i, l1, l2, pax, pbx) * fact2(2 * i - 1)
            / (2.0 * gamma).powi(i);
    }
    total
}

/// Compute kinetic energy integral between two contracted Gaussian basis functions
///
/// The kinetic energy integral is: T_ij = ⟨φ_i|-½∇²|φ_j⟩
pub fn kinetic(bra: &CGBF, ket: &CGBF) -> f64 {
    let mut sum = 0.0;
    for p1 in &bra.primitives {
        for p2 in &ket.primitives {
            sum += p1.coefficient
                * p2.coefficient
                * kinetic_primitive(
                    p1.exponent,
                    bra.shell,
                    bra.origin,
                    p2.exponent,
                    ket.shell,
                    ket.origin,
                );
        }
    }
    sum
}

/// Compute kinetic energy integral between two primitive Gaussian basis functions
///
/// The full form of the kinetic energy integral
fn kinetic_primitive(
    alpha1: f64,
    lmn1: (i32, i32, i32),
    a: Vec3,
    alpha2: f64,
    lmn2: (i32, i32, i32),
    b: Vec3,
) -> f64 {
    let (l1, m1, n1) = lmn1;
    let (l2, m2, n2) = lmn2;

    let term0 = alpha2 * (2.0 * (l2 + m2 + n2) as f64 + 3.0)
        * overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2, m2, n2), b);

    let term1 = -2.0 * alpha2.powi(2)
        * (overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2 + 2, m2, n2), b)
            + overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2, m2 + 2, n2), b)
            + overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2, m2, n2 + 2), b));

    let term2 = -0.5
        * (l2 as f64
            * (l2 - 1) as f64
            * overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2 - 2, m2, n2), b)
            + m2 as f64
                * (m2 - 1) as f64
                * overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2, m2 - 2, n2), b)
            + n2 as f64
                * (n2 - 1) as f64
                * overlap_primitive(alpha1, (l1, m1, n1), a, alpha2, (l2, m2, n2 - 2), b));

    term0 + term1 + term2
}

/// Compute nuclear attraction integral between two basis functions and a point charge
///
/// The nuclear attraction integral is: V_ij = ⟨φ_i|-Z/r|φ_j⟩
pub fn nuclear_attraction(bra: &CGBF, ket: &CGBF, center: Vec3) -> f64 {
    let mut sum = 0.0;
    for p1 in &bra.primitives {
        for p2 in &ket.primitives {
            sum += p1.coefficient
                * p2.coefficient
                * nuclear_attraction_primitive(
                    p1.exponent,
                    bra.shell,
                    bra.origin,
                    p2.exponent,
                    ket.shell,
                    ket.origin,
                    center,
                );
        }
    }
    sum
}

/// Compute nuclear attraction integral between two primitive Gaussians
///
/// Full form of the nuclear attraction integral
fn nuclear_attraction_primitive(
    alpha1: f64,
    lmn1: (i32, i32, i32),
    a: Vec3,
    alpha2: f64,
    lmn2: (i32, i32, i32),
    b: Vec3,
    c: Vec3,
) -> f64 {
    let (l1, m1, n1) = lmn1;
    let (l2, m2, n2) = lmn2;

    let gamma = alpha1 + alpha2;

    let px = gaussian_product_center(alpha1, a.x, alpha2, b.x);
    let py = gaussian_product_center(alpha1, a.y, alpha2, b.y);
    let pz = gaussian_product_center(alpha1, a.z, alpha2, b.z);
    let p = Vec3::new(px, py, pz);

    let rab2 = a.distance_squared(&b);
    let rcp2 = c.distance_squared(&p);

    let ax = a_array(l1, l2, px - a.x, px - b.x, px - c.x, gamma);
    let ay = a_array(m1, m2, py - a.y, py - b.y, py - c.y, gamma);
    let az = a_array(n1, n2, pz - a.z, pz - b.z, pz - c.z, gamma);

    let mut total = 0.0;
    for i in 0..=(l1 + l2) as usize {
        for j in 0..=(m1 + m2) as usize {
            for k in 0..=(n1 + n2) as usize {
                total += ax[i] * ay[j] * az[k] * fgamma((i + j + k) as i32, rcp2 * gamma);
            }
        }
    }

    -2.0 * PI / gamma * (-alpha1 * alpha2 * rab2 / gamma).exp() * total
}

/// Compute the A array used in nuclear attraction integrals (THO eq. 2.18, 3.1)
fn a_array(l1: i32, l2: i32, pa: f64, pb: f64, cp: f64, gamma: f64) -> Vec<f64> {
    let imax = (l1 + l2 + 1) as usize;
    let mut a = vec![0.0; imax];

    for i in 0..imax {
        let i_i32 = i as i32;
        for r in 0..=(i_i32 / 2) {
            for u in 0..=((i_i32 - 2 * r) / 2) {
                let idx = (i_i32 - 2 * r - u) as usize;
                a[idx] += a_term(i_i32, r, u, l1, l2, pa, pb, cp, gamma);
            }
        }
    }

    a
}

/// Compute a single A term (THO eq. 2.18)
fn a_term(i: i32, r: i32, u: i32, l1: i32, l2: i32, pa: f64, pb: f64, cp: f64, gamma: f64) -> f64 {
    use crate::utils::factorial;

    (-1.0_f64).powi(i)
        * binomial_prefactor(i, l1, l2, pa, pb)
        * (-1.0_f64).powi(u)
        * factorial(i)
        * cp.powi(i - 2 * r - 2 * u)
        * (0.25 / gamma).powi(r + u)
        / factorial(r)
        / factorial(u)
        / factorial(i - 2 * r - 2 * u)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Primitive;
    use approx::assert_relative_eq;

    /// Create a simple s-orbital primitive (l=m=n=0)
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
    fn test_overlap_1d() {
        // From PyQuante2 doctest: overlap1d(0,0,0,0,1) = 1.0
        let result = overlap_1d(0, 0, 0.0, 0.0, 1.0);
        assert_relative_eq!(result, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_overlap_primitive_s_orbitals() {
        // From PyQuante2 doctest:
        // overlap(1,(0,0,0),array((0,0,0),'d'),1,(0,0,0),array((0,0,0),'d')) = 1.968701
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let result = overlap_primitive(1.0, (0, 0, 0), origin, 1.0, (0, 0, 0), origin);
        assert_relative_eq!(result, 1.968701, epsilon = 1e-5);
    }

    #[test]
    fn test_overlap_contracted() {
        // Test with contracted Gaussian (same as primitive case with single primitive)
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s1 = s_orbital(1.0, origin);
        let s2 = s_orbital(1.0, origin);
        let result = overlap(&s1, &s2);
        assert_relative_eq!(result, 1.968701, epsilon = 1e-5);
    }

    #[test]
    fn test_kinetic_primitive_s_orbitals() {
        // From PyQuante2 doctest:
        // kinetic(1,(0,0,0),array((0,0,0),'d'),1,(0,0,0),array((0,0,0),'d')) = 2.953052
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let result = kinetic_primitive(1.0, (0, 0, 0), origin, 1.0, (0, 0, 0), origin);
        assert_relative_eq!(result, 2.953052, epsilon = 1e-5);
    }

    #[test]
    fn test_kinetic_contracted() {
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s1 = s_orbital(1.0, origin);
        let s2 = s_orbital(1.0, origin);
        let result = kinetic(&s1, &s2);
        assert_relative_eq!(result, 2.953052, epsilon = 1e-5);
    }

    #[test]
    fn test_nuclear_attraction_primitive() {
        // From PyQuante2 doctest:
        // nuclear_attraction(1,(0,0,0),array((0,0,0),'d'),1,(0,0,0),array((0,0,0),'d'),array((0,0,0),'d')) = -3.141593
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let result = nuclear_attraction_primitive(1.0, (0, 0, 0), origin, 1.0, (0, 0, 0), origin, origin);
        assert_relative_eq!(result, -3.141593, epsilon = 1e-5);
    }

    #[test]
    fn test_nuclear_attraction_contracted() {
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s1 = s_orbital(1.0, origin);
        let s2 = s_orbital(1.0, origin);
        let result = nuclear_attraction(&s1, &s2, origin);
        assert_relative_eq!(result, -3.141593, epsilon = 1e-5);
    }

    #[test]
    fn test_a_array() {
        // From PyQuante2 doctest: A_array(0,0,0,0,0,1) = [1.0]
        let result = a_array(0, 0, 0.0, 0.0, 0.0, 1.0);
        assert_eq!(result.len(), 1);
        assert_relative_eq!(result[0], 1.0, epsilon = 1e-10);
    }
}

