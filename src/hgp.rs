//! Two-electron repulsion integrals using the Head-Gordon and Pople method.
//!
//! This implementation uses the algorithm developed by Head-Gordon and Pople,
//! which is particularly efficient for medium to high angular momentum cases.
//! It relies on the Boys function for computing auxiliary integrals.
//!
//! Reference: PyQuante2 - pyquante2/ints/hgp.py
//! Based on Saika and Obara scheme with Head-Gordon & Pople optimizations

use crate::common::{CGBF, Vec3};
use crate::utils::{fgamma, gaussian_normalization, gaussian_product_center};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Compute two-electron repulsion integral using Head-Gordon-Pople method
///
/// The electron repulsion integral is: (ij|kl) = ⟨φ_i φ_j|1/r₁₂|φ_k φ_l⟩
pub fn electron_repulsion_hgp(bra1: &CGBF, bra2: &CGBF, ket1: &CGBF, ket2: &CGBF) -> f64 {
    let mut sum = 0.0;

    // Get angular momentum
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
                        * hrr(
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

/// Horizontal Recursion Relation (HRR)
///
/// Transfers angular momentum from center b to center a, and from center d to center c
/// Recursively reduces to VRR base cases
#[allow(clippy::too_many_arguments)]
fn hrr(
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

    // Transfer angular momentum from b to a (x-direction)
    if lb > 0 {
        return hrr(a, (la + 1, ma, na), alpha_a, b, (lb - 1, mb, nb), alpha_b, c, lmn_c, alpha_c, d, lmn_d, alpha_d)
            + (a.x - b.x) * hrr(a, lmn_a, alpha_a, b, (lb - 1, mb, nb), alpha_b, c, lmn_c, alpha_c, d, lmn_d, alpha_d);
    }
    // Transfer angular momentum from b to a (y-direction)
    else if mb > 0 {
        return hrr(a, (la, ma + 1, na), alpha_a, b, (lb, mb - 1, nb), alpha_b, c, lmn_c, alpha_c, d, lmn_d, alpha_d)
            + (a.y - b.y) * hrr(a, lmn_a, alpha_a, b, (lb, mb - 1, nb), alpha_b, c, lmn_c, alpha_c, d, lmn_d, alpha_d);
    }
    // Transfer angular momentum from b to a (z-direction)
    else if nb > 0 {
        return hrr(a, (la, ma, na + 1), alpha_a, b, (lb, mb, nb - 1), alpha_b, c, lmn_c, alpha_c, d, lmn_d, alpha_d)
            + (a.z - b.z) * hrr(a, lmn_a, alpha_a, b, (lb, mb, nb - 1), alpha_b, c, lmn_c, alpha_c, d, lmn_d, alpha_d);
    }
    // Transfer angular momentum from d to c (x-direction)
    else if ld > 0 {
        return hrr(a, lmn_a, alpha_a, b, lmn_b, alpha_b, c, (lc + 1, mc, nc), alpha_c, d, (ld - 1, md, nd), alpha_d)
            + (c.x - d.x) * hrr(a, lmn_a, alpha_a, b, lmn_b, alpha_b, c, lmn_c, alpha_c, d, (ld - 1, md, nd), alpha_d);
    }
    // Transfer angular momentum from d to c (y-direction)
    else if md > 0 {
        return hrr(a, lmn_a, alpha_a, b, lmn_b, alpha_b, c, (lc, mc + 1, nc), alpha_c, d, (ld, md - 1, nd), alpha_d)
            + (c.y - d.y) * hrr(a, lmn_a, alpha_a, b, lmn_b, alpha_b, c, lmn_c, alpha_c, d, (ld, md - 1, nd), alpha_d);
    }
    // Transfer angular momentum from d to c (z-direction)
    else if nd > 0 {
        return hrr(a, lmn_a, alpha_a, b, lmn_b, alpha_b, c, (lc, mc, nc + 1), alpha_c, d, (ld, md, nd - 1), alpha_d)
            + (c.z - d.z) * hrr(a, lmn_a, alpha_a, b, lmn_b, alpha_b, c, lmn_c, alpha_c, d, (ld, md, nd - 1), alpha_d);
    }

    // Base case: all angular momentum on centers a and c
    // Use VRR to compute
    vrr(a, lmn_a, alpha_a, b, alpha_b, c, lmn_c, alpha_c, d, alpha_d, 0)
}

/// Vertical Recursion Relation (VRR)
///
/// Builds up integrals with angular momentum using Boys function
/// Computes integrals of the form (la,ma,na,0|lc,mc,nc,0)_m
#[allow(clippy::too_many_arguments)]
fn vrr(
    a: Vec3,
    lmn_a: (i32, i32, i32),
    alpha_a: f64,
    b: Vec3,
    alpha_b: f64,
    c: Vec3,
    lmn_c: (i32, i32, i32),
    alpha_c: f64,
    d: Vec3,
    alpha_d: f64,
    m: i32,
) -> f64 {
    let (la, ma, na) = lmn_a;
    let (lc, mc, nc) = lmn_c;

    // Gaussian product centers
    let px = gaussian_product_center(alpha_a, a.x, alpha_b, b.x);
    let py = gaussian_product_center(alpha_a, a.y, alpha_b, b.y);
    let pz = gaussian_product_center(alpha_a, a.z, alpha_b, b.z);
    let p = Vec3::new(px, py, pz);

    let qx = gaussian_product_center(alpha_c, c.x, alpha_d, d.x);
    let qy = gaussian_product_center(alpha_c, c.y, alpha_d, d.y);
    let qz = gaussian_product_center(alpha_c, c.z, alpha_d, d.z);
    let q = Vec3::new(qx, qy, qz);

    let zeta = alpha_a + alpha_b;
    let eta = alpha_c + alpha_d;

    let wx = gaussian_product_center(zeta, px, eta, qx);
    let wy = gaussian_product_center(zeta, py, eta, qy);
    let wz = gaussian_product_center(zeta, pz, eta, qz);
    let _w = Vec3::new(wx, wy, wz);

    let mtot = la + ma + na + lc + mc + nc + m;

    // Initialize VRR cache
    let mut cache: HashMap<(i32, i32, i32, i32, i32, i32, i32), f64> = HashMap::new();

    // Compute base case (0,0,0,0,0,0,m) for all needed m values
    let rab2 = a.distance_squared(&b);
    let kab = (2.0_f64).sqrt() * PI.powf(1.25) / zeta
        * (-alpha_a * alpha_b / zeta * rab2).exp();

    let rcd2 = c.distance_squared(&d);
    let kcd = (2.0_f64).sqrt() * PI.powf(1.25) / eta
        * (-alpha_c * alpha_d / eta * rcd2).exp();

    let rpq2 = p.distance_squared(&q);
    let t = zeta * eta / (zeta + eta) * rpq2;

    // Compute Boys function values
    let mut fg_terms = vec![0.0; (mtot + 1) as usize];
    fg_terms[mtot as usize] = fgamma(mtot, t);
    for im in (0..mtot).rev() {
        fg_terms[im as usize] = (2.0 * t * fg_terms[(im + 1) as usize] + (-t).exp()) / (2.0 * im as f64 + 1.0);
    }

    // Base cases
    for im in 0..=mtot {
        cache.insert((0, 0, 0, 0, 0, 0, im), kab * kcd / (zeta + eta).sqrt() * fg_terms[im as usize]);
    }

    // Build up x-direction (la)
    for i in 0..la {
        for im in 0..=(mtot - i - 1) {
            let mut val = (px - a.x) * cache[&(i, 0, 0, 0, 0, 0, im)]
                + (wx - px) * cache[&(i, 0, 0, 0, 0, 0, im + 1)];

            if i > 0 {
                val += i as f64 / (2.0 * zeta)
                    * (cache[&(i - 1, 0, 0, 0, 0, 0, im)]
                        - eta / (zeta + eta) * cache[&(i - 1, 0, 0, 0, 0, 0, im + 1)]);
            }

            cache.insert((i + 1, 0, 0, 0, 0, 0, im), val);
        }
    }

    // Build up y-direction (ma)
    for j in 0..ma {
        for i in 0..=la {
            for im in 0..=(mtot - i - j - 1) {
                let mut val = (py - a.y) * cache[&(i, j, 0, 0, 0, 0, im)]
                    + (wy - py) * cache[&(i, j, 0, 0, 0, 0, im + 1)];

                if j > 0 {
                    val += j as f64 / (2.0 * zeta)
                        * (cache[&(i, j - 1, 0, 0, 0, 0, im)]
                            - eta / (zeta + eta) * cache[&(i, j - 1, 0, 0, 0, 0, im + 1)]);
                }

                cache.insert((i, j + 1, 0, 0, 0, 0, im), val);
            }
        }
    }

    // Build up z-direction (na)
    for k in 0..na {
        for j in 0..=ma {
            for i in 0..=la {
                for im in 0..=(mtot - i - j - k - 1) {
                    let mut val = (pz - a.z) * cache[&(i, j, k, 0, 0, 0, im)]
                        + (wz - pz) * cache[&(i, j, k, 0, 0, 0, im + 1)];

                    if k > 0 {
                        val += k as f64 / (2.0 * zeta)
                            * (cache[&(i, j, k - 1, 0, 0, 0, im)]
                                - eta / (zeta + eta) * cache[&(i, j, k - 1, 0, 0, 0, im + 1)]);
                    }

                    cache.insert((i, j, k + 1, 0, 0, 0, im), val);
                }
            }
        }
    }

    // Build up c-center x-direction (lc)
    for ic in 0..lc {
        for k in 0..=na {
            for j in 0..=ma {
                for i in 0..=la {
                    for im in 0..=(mtot - i - j - k - ic - 1) {
                        let mut val = (qx - c.x) * cache[&(i, j, k, ic, 0, 0, im)]
                            + (wx - qx) * cache[&(i, j, k, ic, 0, 0, im + 1)];

                        if ic > 0 {
                            val += ic as f64 / (2.0 * eta)
                                * (cache[&(i, j, k, ic - 1, 0, 0, im)]
                                    - zeta / (zeta + eta) * cache[&(i, j, k, ic - 1, 0, 0, im + 1)]);
                        }

                        if i > 0 {
                            val += i as f64 / (2.0 * (zeta + eta)) * cache[&(i - 1, j, k, ic, 0, 0, im + 1)];
                        }

                        cache.insert((i, j, k, ic + 1, 0, 0, im), val);
                    }
                }
            }
        }
    }

    // Build up c-center y-direction (mc)
    for jc in 0..mc {
        for ic in 0..=lc {
            for k in 0..=na {
                for j in 0..=ma {
                    for i in 0..=la {
                        for im in 0..=(mtot - i - j - k - ic - jc - 1) {
                            let mut val = (qy - c.y) * cache[&(i, j, k, ic, jc, 0, im)]
                                + (wy - qy) * cache[&(i, j, k, ic, jc, 0, im + 1)];

                            if jc > 0 {
                                val += jc as f64 / (2.0 * eta)
                                    * (cache[&(i, j, k, ic, jc - 1, 0, im)]
                                        - zeta / (zeta + eta) * cache[&(i, j, k, ic, jc - 1, 0, im + 1)]);
                            }

                            if j > 0 {
                                val += j as f64 / (2.0 * (zeta + eta)) * cache[&(i, j - 1, k, ic, jc, 0, im + 1)];
                            }

                            cache.insert((i, j, k, ic, jc + 1, 0, im), val);
                        }
                    }
                }
            }
        }
    }

    // Build up c-center z-direction (nc)
    for kc in 0..nc {
        for jc in 0..=mc {
            for ic in 0..=lc {
                for k in 0..=na {
                    for j in 0..=ma {
                        for i in 0..=la {
                            for im in 0..=(mtot - i - j - k - ic - jc - kc - 1) {
                                let mut val = (qz - c.z) * cache[&(i, j, k, ic, jc, kc, im)]
                                    + (wz - qz) * cache[&(i, j, k, ic, jc, kc, im + 1)];

                                if kc > 0 {
                                    val += kc as f64 / (2.0 * eta)
                                        * (cache[&(i, j, k, ic, jc, kc - 1, im)]
                                            - zeta / (zeta + eta) * cache[&(i, j, k, ic, jc, kc - 1, im + 1)]);
                                }

                                // Fixed: check j > 0 before accessing j - 1
                                if j > 0 {
                                    val += j as f64 / (2.0 * (zeta + eta)) * cache[&(i, j - 1, k, ic, jc, kc, im + 1)];
                                }

                                cache.insert((i, j, k, ic, jc, kc + 1, im), val);
                            }
                        }
                    }
                }
            }
        }
    }

    // Return the final value
    *cache.get(&(la, ma, na, lc, mc, nc, m)).unwrap_or(&0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::Primitive;
    use approx::assert_relative_eq;

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
    fn test_hgp_s_orbitals_same() {
        let origin = Vec3::new(0.0, 0.0, 0.0);
        let s = s_orbital(1.0, origin);
        let result = electron_repulsion_hgp(&s, &s, &s, &s);
        assert_relative_eq!(result, 1.128379, epsilon = 1e-5);
    }

    #[test]
    fn test_hgp_separated_orbitals() {
        let p1 = Vec3::new(0.0, 0.0, 0.0);
        let p2 = Vec3::new(0.0, 0.0, 1.0);
        let s1 = s_orbital(1.0, p1);
        let s2 = s_orbital(1.0, p2);
        let result = electron_repulsion_hgp(&s1, &s1, &s2, &s2);
        assert_relative_eq!(result, 0.842701, epsilon = 1e-5);
    }

    #[test]
    fn test_hgp_px_orbitals() {
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
        let result = electron_repulsion_hgp(&s, &s, &px, &px);
        assert_relative_eq!(result, 0.940315972579, epsilon = 1e-5);
    }
}

