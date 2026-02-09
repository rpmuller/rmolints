//! Benchmark comparing original HGP vs VRR-contracted HGP
//!
//! This benchmark measures the performance difference between:
//! - Original: HRR called O(n_primitives^4) times
//! - Contracted: VRR tensors accumulated, HRR called once
//!
//! Expected speedup increases with number of primitives per CGBF.

use rmolints::common::{CGBF, Primitive, Vec3};
use rmolints::hgp::{electron_repulsion_hgp, electron_repulsion_hgp_contracted};
use std::time::Instant;

fn create_cgbf(exponents: &[f64], coefficients: &[f64], shell: (i32, i32, i32), origin: Vec3) -> CGBF {
    let primitives = exponents
        .iter()
        .zip(coefficients.iter())
        .map(|(&exp, &coef)| Primitive {
            exponent: exp,
            coefficient: coef,
        })
        .collect();
    CGBF {
        origin,
        shell,
        primitives,
    }
}

fn benchmark_quartet(name: &str, bra1: &CGBF, bra2: &CGBF, ket1: &CGBF, ket2: &CGBF, n_iterations: usize) {
    let n_primitives = bra1.primitives.len() * bra2.primitives.len() * ket1.primitives.len() * ket2.primitives.len();

    println!("\n{}", name);
    println!("  Primitive quartets: {}", n_primitives);
    println!("  Shell types: {:?} {:?} {:?} {:?}", bra1.shell, bra2.shell, ket1.shell, ket2.shell);

    // Warmup
    for _ in 0..10 {
        electron_repulsion_hgp(bra1, bra2, ket1, ket2);
        electron_repulsion_hgp_contracted(bra1, bra2, ket1, ket2);
    }

    // Benchmark original
    let start = Instant::now();
    let mut sum = 0.0;
    for _ in 0..n_iterations {
        sum += electron_repulsion_hgp(bra1, bra2, ket1, ket2);
    }
    let time_original = start.elapsed();

    // Benchmark contracted
    let start = Instant::now();
    let mut sum_contracted = 0.0;
    for _ in 0..n_iterations {
        sum_contracted += electron_repulsion_hgp_contracted(bra1, bra2, ket1, ket2);
    }
    let time_contracted = start.elapsed();

    // Verify results match
    let diff = (sum - sum_contracted).abs() / n_iterations as f64;
    println!("  Result difference: {:.2e}", diff);

    let time_original_us = time_original.as_micros() as f64 / n_iterations as f64;
    let time_contracted_us = time_contracted.as_micros() as f64 / n_iterations as f64;
    let speedup = time_original_us / time_contracted_us;

    println!("  Original:   {:8.2} µs/integral", time_original_us);
    println!("  Contracted: {:8.2} µs/integral", time_contracted_us);
    println!("  Speedup:    {:8.2}x ({}% faster)", speedup, ((speedup - 1.0) * 100.0) as i32);

    if diff > 1e-10 {
        println!("  ⚠️  WARNING: Results differ by more than 1e-10!");
    }
}

fn main() {
    println!("VRR-Level Contraction Optimization Benchmark");
    println!("============================================");

    let origin = Vec3::new(0.0, 0.0, 0.0);
    let n_iterations = 10000;

    // Test 1: 2 primitives (4 quartets) - minimal contraction
    let exps_2 = vec![1.0, 0.5];
    let coefs_2 = vec![0.6, 0.4];
    let s_2 = create_cgbf(&exps_2, &coefs_2, (0, 0, 0), origin);
    benchmark_quartet("2 primitives (s-s-s-s)", &s_2, &s_2, &s_2, &s_2, n_iterations);

    // Test 2: 3 primitives (81 quartets) - typical STO-3G
    let exps_3 = vec![3.42525091, 0.62391373, 0.16885540];
    let coefs_3 = vec![0.15432897, 0.53532814, 0.44463454];
    let s_3 = create_cgbf(&exps_3, &coefs_3, (0, 0, 0), origin);
    benchmark_quartet("3 primitives (s-s-s-s, STO-3G-like)", &s_3, &s_3, &s_3, &s_3, n_iterations);

    // Test 3: 6 primitives (1296 quartets) - typical 6-31G
    let exps_6 = vec![5484.6717, 825.23495, 188.04696, 52.9645, 16.89757, 5.799635];
    let coefs_6 = vec![0.0018311, 0.0139501, 0.0684451, 0.2327143, 0.4701930, 0.3585209];
    let s_6 = create_cgbf(&exps_6, &coefs_6, (0, 0, 0), origin);
    benchmark_quartet("6 primitives (s-s-s-s, 6-31G-like)", &s_6, &s_6, &s_6, &s_6, n_iterations / 2);

    // Test 4: p-orbitals with 3 primitives
    let px_3 = create_cgbf(&exps_3, &coefs_3, (1, 0, 0), origin);
    benchmark_quartet("3 primitives (s-s-px-px)", &s_3, &s_3, &px_3, &px_3, n_iterations);

    // Test 5: Mixed p-orbitals with 3 primitives
    let py_3 = create_cgbf(&exps_3, &coefs_3, (0, 1, 0), origin);
    benchmark_quartet("3 primitives (px-px-py-py)", &px_3, &px_3, &py_3, &py_3, n_iterations);

    // Test 6: Higher angular momentum
    let dxx_3 = create_cgbf(&exps_3, &coefs_3, (2, 0, 0), origin);
    benchmark_quartet("3 primitives (s-s-dxx-dxx)", &s_3, &s_3, &dxx_3, &dxx_3, n_iterations / 2);

    // Test 7: Separated orbitals
    let origin2 = Vec3::new(0.0, 0.0, 1.5);
    let s_3_sep = create_cgbf(&exps_3, &coefs_3, (0, 0, 0), origin2);
    benchmark_quartet("3 primitives separated (s-s-s-s)", &s_3, &s_3, &s_3_sep, &s_3_sep, n_iterations);

    println!("\n============================================");
    println!("Summary:");
    println!("- Speedup increases with number of primitives");
    println!("- Best gains with 6+ primitives (typical 6-31G)");
    println!("- Higher angular momentum benefits more (larger VRR)");
}
