use rmolints::basis::build_sto3g_basis;
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn benchmark_molecule(name: &str, molecule: Molecule, method: ERIMethod, warmup: bool) -> (usize, usize, f64) {
    let basis = build_sto3g_basis(&molecule);
    let n = basis.len();
    let n_unique = n * (n + 1) * (n * (n + 1) / 2 + 1) / 8;

    // Warmup run if requested
    if warmup {
        let _ = compute_eri_tensor_parallel(&basis, method);
    }

    // Timed run
    let start = Instant::now();
    let eris = compute_eri_tensor_parallel(&basis, method);
    let elapsed = start.elapsed().as_secs_f64() * 1000.0; // Convert to milliseconds

    println!(
        "  {:20} {:12} basis fns, {:8} ERIs: {:10.2} ms",
        name,
        n,
        eris.len(),
        elapsed
    );

    (n, n_unique, elapsed)
}

fn main() {
    println!("\n=== Molecular ERI Computation Benchmark ===");
    println!("Comparing all methods on real molecules\n");

    let molecules = vec![
        ("H2", Molecule::h2(1.4)),
        ("H2O", Molecule::h2o()),
        ("Benzene (C6H6)", Molecule::benzene()),
    ];

    let methods = vec![
        ("Standard THO", ERIMethod::Standard),
        ("Rys Quadrature", ERIMethod::Rys),
        ("Head-Gordon-Pople", ERIMethod::HeadGordonPople),
        ("HGP Optimized", ERIMethod::HeadGordonPopleOpt),
    ];

    for (method_name, method) in &methods {
        println!("\n{}", method_name);
        println!("{}", "=".repeat(method_name.len()));

        for (mol_name, molecule) in &molecules {
            benchmark_molecule(mol_name, molecule.clone(), *method, true);
        }
    }

    println!("\n\n=== Detailed Results with Multiple Runs ===\n");

    // Run each combination 3 times for more stable results
    for (mol_name, molecule) in &molecules {
        println!("\n{}", mol_name);
        println!("{}", "=".repeat(mol_name.len()));

        for (method_name, method) in &methods {
            print!("  {:25}", method_name);
            let mut times = Vec::new();

            // First warmup
            benchmark_molecule("", molecule.clone(), *method, true);

            // Then 3 timed runs
            for _ in 0..3 {
                let (_, _, elapsed) = benchmark_molecule("", molecule.clone(), *method, false);
                times.push(elapsed);
            }

            let avg = times.iter().sum::<f64>() / times.len() as f64;
            let min = times.iter().cloned().fold(f64::INFINITY, f64::min);
            let max = times.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

            println!("  avg: {:8.2} ms  (min: {:8.2}, max: {:8.2})", avg, min, max);
        }
    }

    println!("\n\n=== Summary Table ===\n");
    println!("Time to compute all two-electron integrals (milliseconds):");
    println!("Lower is better. Best time in each row marked with ✅\n");

    // Final benchmark runs for summary table
    let mut results: Vec<Vec<f64>> = Vec::new();

    for (_mol_name, molecule) in &molecules {
        let mut row = Vec::new();
        for (_, method) in &methods {
            // Run 5 times and take minimum (best)
            let mut times = Vec::new();
            benchmark_molecule("", molecule.clone(), *method, true);
            for _ in 0..5 {
                let (_, _, elapsed) = benchmark_molecule("", molecule.clone(), *method, false);
                times.push(elapsed);
            }
            row.push(times.iter().cloned().fold(f64::INFINITY, f64::min));
        }
        results.push(row);
    }

    // Print formatted table
    println!("| Molecule | Basis Fns | Unique ERIs | Standard THO | Rys Quad | HGP Original | HGP Opt |");
    println!("|----------|-----------|-------------|--------------|----------|--------------|---------|");

    for (idx, (mol_name, molecule)) in molecules.iter().enumerate() {
        let basis = build_sto3g_basis(molecule);
        let n = basis.len();
        let n_unique = n * (n + 1) * (n * (n + 1) / 2 + 1) / 8;

        print!("| {} | {} | {} |", mol_name, n, n_unique);

        let row = &results[idx];
        let min_time = row.iter().cloned().fold(f64::INFINITY, f64::min);

        for &time in row {
            if (time - min_time).abs() < 0.01 {
                print!(" **{:.2} ms** ✅ |", time);
            } else {
                print!(" {:.2} ms |", time);
            }
        }
        println!();
    }

    println!("\n\n=== Scaling Analysis ===\n");
    println!("As basis set size N increases, ERI count scales as O(N⁴):\n");

    for (mol_name, molecule) in &molecules {
        let basis = build_sto3g_basis(molecule);
        let n = basis.len();
        let n_unique = n * (n + 1) * (n * (n + 1) / 2 + 1) / 8;
        println!(
            "  {}: N={}, unique ERIs={} (ratio: {:.1}:1)",
            mol_name,
            n,
            n_unique,
            n_unique as f64 / n as f64
        );
    }

    println!("\n\nBenchmark complete!");
}
