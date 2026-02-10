use rmolints::basis::{build_basis, build_sto3g_basis, BasisSet};
use rmolints::molecule::Molecule;
use rmolints::parallel::{
    compute_eri_tensor_parallel, compute_eri_tensor_screened_parallel, ERIMethod,
    SCHWARZ_THRESHOLD,
};
use std::time::Instant;

fn benchmark_pair(name: &str, basis: &[rmolints::common::CGBF], method: ERIMethod) {
    let n = basis.len();

    // Warmup both code paths
    let _ = compute_eri_tensor_parallel(basis, method);
    let _ = compute_eri_tensor_screened_parallel(basis, method, SCHWARZ_THRESHOLD);

    // Run each 3 times, take best
    let mut times_unscreened = Vec::new();
    let mut times_screened = Vec::new();
    let mut unscreened_count = 0;
    let mut screened_count = 0;

    for _ in 0..3 {
        let start = Instant::now();
        let u = compute_eri_tensor_parallel(basis, method);
        times_unscreened.push(start.elapsed().as_secs_f64() * 1000.0);
        unscreened_count = u.len();

        let start = Instant::now();
        let s = compute_eri_tensor_screened_parallel(basis, method, SCHWARZ_THRESHOLD);
        times_screened.push(start.elapsed().as_secs_f64() * 1000.0);
        screened_count = s.len();

        // Correctness: every screened integral should match unscreened
        let umap: std::collections::HashMap<(usize, usize, usize, usize), f64> =
            u.iter().map(|&(i, j, k, l, v)| ((i, j, k, l), v)).collect();
        let max_err = s.iter().fold(0.0_f64, |acc, &(i, j, k, l, v)| {
            acc.max((v - umap[&(i, j, k, l)]).abs())
        });
        assert!(max_err < 1e-10, "screened result differs: {max_err:.2e}");
    }

    let t_u = times_unscreened.iter().cloned().fold(f64::INFINITY, f64::min);
    let t_s = times_screened.iter().cloned().fold(f64::INFINITY, f64::min);
    let skipped = unscreened_count - screened_count;
    let skip_pct = 100.0 * skipped as f64 / unscreened_count as f64;
    let speedup = t_u / t_s;

    println!(
        "  {:<20} N={:3}  {:7} total  {:6} skipped ({:5.1}%)  \
         unscreened: {:8.2} ms  screened: {:8.2} ms  speedup: {:.2}x",
        name, n, unscreened_count, skipped, skip_pct, t_u, t_s, speedup
    );
}

fn main() {
    println!("\n=== Schwarz Screening Benchmark ===");
    println!("Threshold: {:.0e}  Method: HGP-Contracted\n", SCHWARZ_THRESHOLD);

    let method = ERIMethod::HeadGordonPopleContracted;

    println!("--- STO-3G basis ---");
    for (name, mol) in &[
        ("H2", Molecule::h2(1.4)),
        ("H2O", Molecule::h2o()),
        ("Benzene", Molecule::benzene()),
    ] {
        let basis = build_sto3g_basis(mol);
        benchmark_pair(name, &basis, method);
    }

    println!("\n--- 6-31G basis ---");
    for (name, mol) in &[("H2O", Molecule::h2o()), ("Benzene", Molecule::benzene())] {
        let basis = build_basis(mol, BasisSet::_631G);
        benchmark_pair(name, &basis, method);
    }

    println!("\n--- 6-31G(d) basis ---");
    for (name, mol) in &[("H2O", Molecule::h2o()), ("Benzene", Molecule::benzene())] {
        let basis = build_basis(mol, BasisSet::_631GStar);
        benchmark_pair(name, &basis, method);
    }
}
