use rmolints::basis::{build_basis, BasisSet};
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn benchmark_molecule(_name: &str, molecule: &Molecule, method: ERIMethod) -> (usize, usize, f64) {
    let basis = build_basis(molecule, BasisSet::_631GStarStar);
    let n = basis.len();

    // Calculate number of unique ERIs (8-fold symmetry)
    let n_unique = {
        let mut count = 0;
        for i in 0..n {
            for j in 0..=i {
                let ij = i * (i + 1) / 2 + j;
                for k in 0..=i {
                    let l_max = if k == i { j } else { k };
                    for l in 0..=l_max {
                        let kl = k * (k + 1) / 2 + l;
                        if ij >= kl {
                            count += 1;
                        }
                    }
                }
            }
        }
        count
    };

    // Warmup
    let _ = compute_eri_tensor_parallel(&basis, method);

    // Benchmark with 3 runs
    let mut times = Vec::new();
    for _ in 0..3 {
        let start = Instant::now();
        let _eris = compute_eri_tensor_parallel(&basis, method);
        times.push(start.elapsed().as_secs_f64() * 1000.0);
    }

    let min_time = times.iter().cloned().fold(f64::INFINITY, f64::min);
    (n, n_unique, min_time)
}

fn main() {
    println!("\n=== 6-31G(d,p) Basis Set Performance Benchmark ===");
    println!("Testing realistic basis set used by chemists");
    println!("(d polarization on heavy atoms, p polarization on H)\n");

    let molecules = vec![
        ("H2", Molecule::h2(1.4)),
        ("H2O", Molecule::h2o()),
        ("NH3", Molecule::nh3()),
        ("Benzene", Molecule::benzene()),
    ];

    let methods = vec![
        ("HGP", ERIMethod::HeadGordonPople),
        ("Rys", ERIMethod::Rys),
        ("Standard", ERIMethod::Standard),
    ];

    println!("{:<12} {:<10} {:<12} {:<10} {:<10} {:<10}",
             "Molecule", "Basis Fns", "Unique ERIs", "HGP (ms)", "Rys (ms)", "Std (ms)");
    println!("{}", "-".repeat(75));

    for (mol_name, molecule) in &molecules {
        print!("{:<12} ", mol_name);

        let mut results = Vec::new();
        for (_, method) in &methods {
            let (n, n_unique, time) = benchmark_molecule(mol_name, molecule, *method);
            if results.is_empty() {
                print!("{:<10} {:<12} ", n, n_unique);
            }
            results.push(time);
        }

        for time in results {
            print!("{:<10.2} ", time);
        }
        println!();
    }

    println!("\n=== Performance Ratios (vs HGP) ===\n");
    println!("{:<12} {:<15} {:<15}", "Molecule", "Rys / HGP", "Standard / HGP");
    println!("{}", "-".repeat(45));

    for (mol_name, molecule) in &molecules {
        let (_, _, hgp_time) = benchmark_molecule(mol_name, molecule, ERIMethod::HeadGordonPople);
        let (_, _, rys_time) = benchmark_molecule(mol_name, molecule, ERIMethod::Rys);
        let (_, _, std_time) = benchmark_molecule(mol_name, molecule, ERIMethod::Standard);

        println!("{:<12} {:<15.2}x {:<15.2}x",
                 mol_name,
                 rys_time / hgp_time,
                 std_time / hgp_time);
    }

    println!("\n=== Scaling Analysis ===\n");
    println!("Showing how ERI count scales with basis set size (O(N⁴)):\n");

    for (mol_name, molecule) in &molecules {
        let basis = build_basis(molecule, BasisSet::_631GStarStar);
        let n = basis.len();
        let n_unique = {
            let mut count = 0;
            for i in 0..n {
                for j in 0..=i {
                    let ij = i * (i + 1) / 2 + j;
                    for k in 0..=i {
                        let l_max = if k == i { j } else { k };
                        for l in 0..=l_max {
                            let kl = k * (k + 1) / 2 + l;
                            if ij >= kl {
                                count += 1;
                            }
                        }
                    }
                }
            }
            count
        };
        println!("{:<12}: N={:3}, unique ERIs={:8} (ratio: {:6.1}:1)",
                 mol_name, n, n_unique, n_unique as f64 / n as f64);
    }
}
