use rmolints::basis::build_sto3g_basis;
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn main() {
    println!("\n=== Rys Optimization Comparison ===\n");
    println!("Comparing Rys vs HGP on benzene\n");

    let molecule = Molecule::benzene();
    let basis = build_sto3g_basis(&molecule);
    let n = basis.len();

    println!("Molecule: Benzene (C6H6)");
    println!("Basis functions: {}\n", n);

    // Run each method 5 times to get stable measurements
    let methods = vec![
        ("HGP", ERIMethod::HeadGordonPople),
        ("Rys", ERIMethod::Rys),
        ("Standard", ERIMethod::Standard),
    ];

    for (name, method) in &methods {
        println!("{}", name);
        println!("{}", "=".repeat(name.len()));

        // Warmup
        let _ = compute_eri_tensor_parallel(&basis, *method);

        let mut times = Vec::new();
        for i in 0..5 {
            let start = Instant::now();
            let eris = compute_eri_tensor_parallel(&basis, *method);
            let elapsed = start.elapsed().as_secs_f64() * 1000.0;
            times.push(elapsed);
            println!("  Run {}: {:.2} ms ({} ERIs)", i + 1, elapsed, eris.len());
        }

        let avg = times.iter().sum::<f64>() / times.len() as f64;
        let min = times.iter().cloned().fold(f64::INFINITY, f64::min);
        let max = times.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let stddev = (times.iter().map(|&t| (t - avg).powi(2)).sum::<f64>() / times.len() as f64).sqrt();

        println!("  Average: {:.2} ms", avg);
        println!("  Min: {:.2} ms", min);
        println!("  Max: {:.2} ms", max);
        println!("  Std Dev: {:.2} ms\n", stddev);
    }

    println!("\n=== Summary ===\n");
    println!("Run this benchmark multiple times to verify consistency.");
    println!("Look for the minimum times (best performance).");
}
