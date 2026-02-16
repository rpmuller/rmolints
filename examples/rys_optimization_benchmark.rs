use rmolints::basis::build_sto3g_basis;
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn benchmark_method(name: &str, molecule: &Molecule, method: ERIMethod, n_runs: usize) -> f64 {
    let basis = build_sto3g_basis(molecule);
    let n = basis.len();

    // Warmup
    let _ = compute_eri_tensor_parallel(&basis, method);

    let mut times = Vec::with_capacity(n_runs);
    for _ in 0..n_runs {
        let start = Instant::now();
        let eris = compute_eri_tensor_parallel(&basis, method);
        let elapsed = start.elapsed().as_secs_f64() * 1000.0;
        std::hint::black_box(&eris);
        times.push(elapsed);
    }

    times.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = times[n_runs / 2];
    let min = times[0];
    let max = times[n_runs - 1];

    println!(
        "  {:30} {:3} basis fns  median: {:8.3} ms  (min: {:8.3}, max: {:8.3})",
        name, n, median, min, max
    );
    median
}

fn main() {
    println!("\n=== Rys Optimization Benchmark ===\n");

    let n_runs = 7;

    let molecules: Vec<(&str, Molecule)> = vec![
        ("H2", Molecule::h2(1.4)),
        ("H2O", Molecule::h2o()),
        ("NH3", Molecule::nh3()),
        ("Benzene", Molecule::benzene()),
    ];

    println!("Rys Quadrature:");
    println!("{}", "-".repeat(90));
    let mut rys_times = Vec::new();
    for (name, mol) in &molecules {
        let t = benchmark_method(name, mol, ERIMethod::Rys, n_runs);
        rys_times.push(t);
    }

    println!("\nHGP Contracted (reference):");
    println!("{}", "-".repeat(90));
    let mut hgp_times = Vec::new();
    for (name, mol) in &molecules {
        let t = benchmark_method(name, mol, ERIMethod::HeadGordonPopleContracted, n_runs);
        hgp_times.push(t);
    }

    println!("\n\n=== Rys vs HGP Ratio ===\n");
    println!("| {:12} | {:>10} | {:>10} | {:>8} |", "Molecule", "Rys (ms)", "HGP (ms)", "Ratio");
    println!("|{}|{}|{}|{}|", "-".repeat(14), "-".repeat(12), "-".repeat(12), "-".repeat(10));
    for (i, (name, _)) in molecules.iter().enumerate() {
        let ratio = rys_times[i] / hgp_times[i];
        println!(
            "| {:12} | {:10.3} | {:10.3} | {:8.2}x |",
            name, rys_times[i], hgp_times[i], ratio
        );
    }

    println!("\nBenchmark complete!");
}
