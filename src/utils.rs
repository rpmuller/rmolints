//! Utility functions for integral calculations.
//!
//! Contains mathematical helper functions like factorials, binomials,
//! and special functions needed for computing molecular integrals.

use std::f64::consts::PI;

/// Compute factorial of n
#[inline]
pub fn factorial(n: i32) -> f64 {
    if n <= 1 {
        1.0
    } else {
        (2..=n).fold(1.0, |acc, i| acc * i as f64)
    }
}

/// Compute double factorial of n: n!! = n * (n-2) * (n-4) * ...
#[inline]
pub fn fact2(n: i32) -> f64 {
    if n <= 1 {
        1.0
    } else {
        let mut result = 1.0;
        let mut i = n;
        while i > 1 {
            result *= i as f64;
            i -= 2;
        }
        result
    }
}

/// Compute binomial coefficient C(n, k) = n! / (k! * (n-k)!)
#[inline]
pub fn binomial(n: i32, k: i32) -> f64 {
    if k < 0 || k > n {
        0.0
    } else if k == 0 || k == n {
        1.0
    } else {
        factorial(n) / (factorial(k) * factorial(n - k))
    }
}

/// Compute the binomial prefactor used in integral calculations
///
/// This is from Augspurger and Dykstra, used in overlap and nuclear attraction integrals.
pub fn binomial_prefactor(s: i32, ia: i32, ib: i32, xpa: f64, xpb: f64) -> f64 {
    let mut total = 0.0;
    for t in 0..=s {
        if s - ia <= t && t <= ib {
            total += binomial(ia, s - t)
                * binomial(ib, t)
                * xpa.powi(ia - s + t)
                * xpb.powi(ib - t);
        }
    }
    total
}

/// Compute the normalization constant for a primitive Gaussian
///
/// For a Gaussian with exponent alpha and angular momentum (l, m, n):
/// N = (2α/π)^(3/4) * (4α)^((l+m+n)/2) / sqrt((2l-1)!! * (2m-1)!! * (2n-1)!!)
#[inline]
pub fn gaussian_normalization(alpha: f64, l: i32, m: i32, n: i32) -> f64 {
    use std::f64::consts::PI;

    let l_plus_m_plus_n = l + m + n;
    let numerator = (2.0 * alpha / PI).powf(0.75) * (4.0 * alpha).powi(l_plus_m_plus_n) / 2.0_f64.powi(l_plus_m_plus_n);
    let denominator = (fact2(2 * l - 1) * fact2(2 * m - 1) * fact2(2 * n - 1)).sqrt();

    numerator / denominator
}

/// Compute the center of the Gaussian product
///
/// When two Gaussians centered at A and B with exponents alpha1 and alpha2
/// are multiplied, the resulting Gaussian is centered at P.
#[inline]
pub fn gaussian_product_center(alpha1: f64, a: f64, alpha2: f64, b: f64) -> f64 {
    (alpha1 * a + alpha2 * b) / (alpha1 + alpha2)
}

/// Compute factorial ratio: a! / (b! * (a-2b)!)
///
/// Used in two-electron integral calculations
#[inline]
pub fn fact_ratio2(a: i32, b: i32) -> f64 {
    if a < 0 || b < 0 || a < 2 * b {
        0.0
    } else {
        factorial(a) / (factorial(b) * factorial(a - 2 * b))
    }
}

/// Incomplete gamma function used in nuclear attraction integrals
///
/// F_m(t) = ∫₀¹ u^(2m) exp(-t u²) du
///
/// This is related to the standard incomplete gamma function by:
/// F_m(t) = 0.5 * t^(-m-0.5) * P(m+0.5, t) * Γ(m+0.5)
/// where P is the regularized incomplete gamma function.
pub fn fgamma(m: i32, t: f64) -> f64 {
    const SMALL: f64 = 1e-12;

    let t = t.max(SMALL);  // Avoid numerical issues with very small t
    let a = m as f64 + 0.5;

    // Compute incomplete gamma function
    let gamma_inc = gamm_inc(a, t);

    // Return F_m(t) = 0.5 * t^(-m-0.5) * gamma_inc
    0.5 * t.powf(-a) * gamma_inc
}

/// Regularized incomplete gamma function P(a,x) * Γ(a)
///
/// Uses series representation for x < a+1, continued fractions otherwise.
/// Based on Numerical Recipes and PyQuante2 implementation.
fn gamm_inc(a: f64, x: f64) -> f64 {
    if x < a + 1.0 {
        // Use series representation
        let (gam, gln) = gser(a, x);
        gln.exp() * gam
    } else {
        // Use continued fractions
        let (gamc, gln) = gcf(a, x);
        gln.exp() * (1.0 - gamc)
    }
}

/// Series representation of incomplete gamma function
fn gser(a: f64, x: f64) -> (f64, f64) {
    const ITMAX: usize = 100;
    const EPS: f64 = 3e-7;

    let gln = lgamma(a);

    if x == 0.0 {
        return (0.0, gln);
    }

    let mut ap = a;
    let mut del = 1.0 / a;
    let mut sum = del;

    for _ in 0..ITMAX {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * EPS {
            break;
        }
    }

    let gamser = sum * (-x + a * x.ln() - gln).exp();
    (gamser, gln)
}

/// Continued fraction representation of incomplete gamma function
fn gcf(a: f64, x: f64) -> (f64, f64) {
    const ITMAX: usize = 100;
    const EPS: f64 = 3e-7;
    const FPMIN: f64 = 1e-30;

    let gln = lgamma(a);
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / FPMIN;
    let mut d = 1.0 / b;
    let mut h = d;

    for i in 1..=ITMAX {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < FPMIN {
            d = FPMIN;
        }
        c = b + an / c;
        if c.abs() < FPMIN {
            c = FPMIN;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < EPS {
            break;
        }
    }

    let gammcf = (-x + a * x.ln() - gln).exp() * h;
    (gammcf, gln)
}

/// Natural logarithm of the gamma function
fn lgamma(x: f64) -> f64 {
    // Use Lanczos approximation
    const G: f64 = 7.0;
    const COEF: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];

    if x < 0.5 {
        // Use reflection formula: ln Γ(1-z) = ln(π/sin(πz)) - ln Γ(z)
        let s = PI / ((PI * x).sin());
        s.ln() - lgamma(1.0 - x)
    } else {
        let z = x - 1.0;
        let mut sum = COEF[0];
        for (i, &c) in COEF.iter().enumerate().skip(1) {
            sum += c / (z + i as f64);
        }
        let t = z + G + 0.5;
        (2.0 * PI).sqrt().ln() + (z + 0.5) * t.ln() - t + sum.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_factorial() {
        assert_eq!(factorial(0), 1.0);
        assert_eq!(factorial(1), 1.0);
        assert_eq!(factorial(5), 120.0);
    }

    #[test]
    fn test_fact2() {
        assert_eq!(fact2(-1), 1.0);
        assert_eq!(fact2(0), 1.0);
        assert_eq!(fact2(1), 1.0);
        assert_eq!(fact2(3), 3.0);
        assert_eq!(fact2(5), 15.0); // 5 * 3 * 1
        assert_eq!(fact2(6), 48.0); // 6 * 4 * 2
    }

    #[test]
    fn test_binomial() {
        assert_eq!(binomial(5, 0), 1.0);
        assert_eq!(binomial(5, 5), 1.0);
        assert_eq!(binomial(5, 2), 10.0);
        assert_eq!(binomial(4, 2), 6.0);
    }

    #[test]
    fn test_binomial_prefactor() {
        // From PyQuante2 doctest
        assert_eq!(binomial_prefactor(0, 0, 0, 0.0, 0.0), 1.0);
    }

    #[test]
    fn test_gaussian_product_center() {
        // Equal weights should give midpoint
        let center = gaussian_product_center(1.0, 0.0, 1.0, 2.0);
        assert!((center - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fgamma_at_zero() {
        // F_0(0) = 1
        assert!((fgamma(0, 0.0) - 1.0).abs() < 1e-10);
        // F_1(0) = 1/3
        assert!((fgamma(1, 0.0) - 1.0 / 3.0).abs() < 1e-10);
        // F_2(0) = 1/5
        assert!((fgamma(2, 0.0) - 1.0 / 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_fgamma_small_t() {
        // For small t, should be close to 1/(2m+1)
        let result = fgamma(0, 1e-15);
        assert!((result - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_fact_ratio2() {
        // fact_ratio2(a, b) = a! / (b! * (a-2b)!)
        assert_eq!(fact_ratio2(0, 0), 1.0); // 0! / (0! * 0!) = 1
        assert_eq!(fact_ratio2(2, 1), 2.0); // 2! / (1! * 0!) = 2
        assert_eq!(fact_ratio2(4, 2), 12.0); // 4! / (2! * 0!) = 24/2 = 12
        assert_eq!(fact_ratio2(4, 1), 12.0); // 4! / (1! * 2!) = 24/(1*2) = 12
    }
}
