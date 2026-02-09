//! Example demonstrating the VRR-contracted HGP optimization
//!
//! This shows how to use the optimized electron_repulsion_hgp_contracted()
//! function and the ERIMethod::HeadGordonPopleContracted variant.

use rmolints::basis::{build_basis, BasisSet};
use rmolints::hgp::{electron_repulsion_hgp, electron_repulsion_hgp_contracted};
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use std::time::Instant;

fn main() {
    println!("VRR-Contracted HGP Usage Example");
    println!("=================================\n");

    // Create a water molecule
    let h2o = Molecule::h2o();

    println!("Example 1: Direct function calls");
    println!("---------------------------------");

    // Build basis set
    let basis = build_basis(&h2o, BasisSet::_631GStar);
    println!("H2O with 6-31G* basis: {} functions\n", basis.len());

    // Compare single integral computation
    let bra1 = &basis[0];
    let bra2 = &basis[1];
    let ket1 = &basis[2];
    let ket2 = &basis[3];

    // Original method
    let start = Instant::now();
    let result_original = electron_repulsion_hgp(bra1, bra2, ket1, ket2);
    let time_original = start.elapsed();
    println!("Original HGP:    {:.6} (took {:?})", result_original, time_original);

    // Contracted method
    let start = Instant::now();
    let result_contracted = electron_repulsion_hgp_contracted(bra1, bra2, ket1, ket2);
    let time_contracted = start.elapsed();
    println!("Contracted HGP:  {:.6} (took {:?})", result_contracted, time_contracted);
    println!("Results match:   {}\n", (result_original - result_contracted).abs() < 1e-12);

    println!("Example 2: Parallel computation with ERIMethod");
    println!("-----------------------------------------------");

    // Compute all ERIs using different methods
    let methods = vec![
        ("Standard THO", ERIMethod::Standard),
        ("Rys Quadrature", ERIMethod::Rys),
        ("HGP Original", ERIMethod::HeadGordonPople),
        ("HGP Contracted", ERIMethod::HeadGordonPopleContracted),
    ];

    for (name, method) in methods {
        let start = Instant::now();
        let eris = compute_eri_tensor_parallel(&basis, method);
        let elapsed = start.elapsed();
        println!("{:20} {} integrals in {:.3}s", name, eris.len(), elapsed.as_secs_f64());
    }

    println!("\nExample 3: When to use contracted vs original");
    println!("----------------------------------------------");

    println!("✅ Use CONTRACTED for:");
    println!("  • 6-31G or larger basis sets");
    println!("  • Production calculations");
    println!("  • Full molecule ERI tensors");
    println!("  • Parallel computations\n");

    println!("✅ Use ORIGINAL for:");
    println!("  • Debugging/verification");
    println!("  • Minimal basis sets (STO-3G on small molecules)");
    println!("  • Single integral evaluations\n");

    println!("Example 4: Recommended usage pattern");
    println!("-------------------------------------");
    println!("use rmolints::parallel::{{compute_eri_tensor_parallel, ERIMethod}};");
    println!("use rmolints::basis::{{build_basis, BasisSet}};");
    println!("use rmolints::molecule::Molecule;");
    println!();
    println!("let molecule = Molecule::h2o();");
    println!("let basis = build_basis(&molecule, BasisSet::_631GStar);");
    println!();
    println!("// Use the contracted optimization (recommended)");
    println!("let eris = compute_eri_tensor_parallel(");
    println!("    &basis,");
    println!("    ERIMethod::HeadGordonPopleContracted");
    println!(");");
    println!();
    println!("// Process integrals...");
    println!("for (i, j, k, l, value) in eris {{");
    println!("    // ... your code here");
    println!("}}");
}
