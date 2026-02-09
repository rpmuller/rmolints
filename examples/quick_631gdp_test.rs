use rmolints::basis::{build_basis, BasisSet};
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn main() {
    println!("\n=== Quick 6-31G(d,p) Performance Test ===\n");

    let h2o = Molecule::h2o();
    let basis = build_basis(&h2o, BasisSet::_631GStarStar);

    println!("Molecule: H2O");
    println!("Basis set: 6-31G(d,p)");
    println!("Basis functions: {}", basis.len());

    // Warmup
    println!("\nWarmup...");
    let _ = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPople);

    // Benchmark each method
    let methods = vec![
        ("HGP", ERIMethod::HeadGordonPople),
        ("Rys", ERIMethod::Rys),
        ("Standard", ERIMethod::Standard),
    ];

    println!("\nBenchmarking (3 runs each, reporting minimum):\n");

    for (name, method) in methods {
        let mut times = Vec::new();
        for _ in 0..3 {
            let start = Instant::now();
            let eris = compute_eri_tensor_parallel(&basis, method);
            let elapsed = start.elapsed().as_secs_f64() * 1000.0;
            times.push(elapsed);

            if times.len() == 1 {
                println!("{}: {} unique ERIs", name, eris.len());
            }
        }

        let min_time = times.iter().cloned().fold(f64::INFINITY, f64::min);
        let avg_time = times.iter().sum::<f64>() / times.len() as f64;

        println!("  Min: {:.2} ms, Avg: {:.2} ms", min_time, avg_time);
    }
}
