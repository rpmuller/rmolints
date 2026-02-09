use rmolints::basis::{basis_labels, build_sto3g_basis, BasisSet};
use rmolints::molecule::Molecule;
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use rmolints::{one_electron, two_electron};
use std::time::Instant;

fn main() {
    println!("\n=== Benzene (C6H6) Molecular Integrals (STO-3G Basis) ===\n");

    // Create benzene molecule
    let benzene = Molecule::benzene();
    println!("Molecule: C6H6 (Benzene)");
    println!("  Atoms: {} (6 C + 6 H)", benzene.natoms());
    println!("  Electrons: {}", benzene.nelectrons());
    println!("  Occupied orbitals: {}", benzene.nocc());
    println!();

    // Print geometry (first few atoms)
    println!("Geometry (Angstroms) - first 6 atoms:");
    for (i, atom) in benzene.atoms.iter().enumerate().take(6) {
        println!(
            "  {:2} {:2}  {:8.4}  {:8.4}  {:8.4}",
            i + 1,
            atom.element.symbol,
            atom.position.x,
            atom.position.y,
            atom.position.z
        );
    }
    println!("  ... and 6 H atoms");
    println!();

    // Build basis set
    let basis = build_sto3g_basis(&benzene);
    let labels = basis_labels(&basis, &benzene, BasisSet::STO3G);

    println!("Basis set: STO-3G");
    println!("  Basis functions: {}", basis.len());
    println!("  Per C atom: 5 functions (1s, 2s, 2px, 2py, 2pz)");
    println!("  Per H atom: 1 function (1s)");
    println!();

    // Print basis function summary
    println!("Basis function summary:");
    println!("  C atoms (1-6):  functions 1-30");
    println!("  H atoms (7-12): functions 31-36");
    println!();

    let n = basis.len();

    // Compute overlap matrix
    println!("Computing overlap matrix S ({n}x{n})...");
    let start = Instant::now();
    let mut s_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let s_ij = one_electron::overlap(&basis[i], &basis[j]);
            s_matrix[i][j] = s_ij;
            s_matrix[j][i] = s_ij;
        }
    }
    let s_time = start.elapsed();
    println!("  Time: {:.2} ms", s_time.as_secs_f64() * 1000.0);
    println!("  Diagonal elements (normalization check):");
    for i in [0, 5, 6, 11, 30, 35] {
        if i < n {
            println!("    S[{:2},{:2}] = {:8.4}  ({})", i+1, i+1, s_matrix[i][i], labels[i]);
        }
    }
    println!();

    // Compute kinetic energy matrix
    println!("Computing kinetic energy matrix T ({n}x{n})...");
    let start = Instant::now();
    let mut t_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let t_ij = one_electron::kinetic(&basis[i], &basis[j]);
            t_matrix[i][j] = t_ij;
            t_matrix[j][i] = t_ij;
        }
    }
    let t_time = start.elapsed();
    println!("  Time: {:.2} ms", t_time.as_secs_f64() * 1000.0);
    println!();

    // Compute nuclear attraction matrix
    println!("Computing nuclear attraction matrix V ({n}x{n})...");
    let start = Instant::now();
    let mut v_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let mut v_ij = 0.0;
            for atom in &benzene.atoms {
                v_ij += one_electron::nuclear_attraction(&basis[i], &basis[j], atom.position)
                    * atom.element.charge();
            }
            v_matrix[i][j] = v_ij;
            v_matrix[j][i] = v_ij;
        }
    }
    let v_time = start.elapsed();
    println!("  Time: {:.2} ms", v_time.as_secs_f64() * 1000.0);
    println!();

    // Core Hamiltonian
    println!("Computing core Hamiltonian H_core = T + V...");
    let mut h_core = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            h_core[i][j] = t_matrix[i][j] + v_matrix[i][j];
        }
    }
    println!("  Core Hamiltonian diagonal elements (sample):");
    for i in [0, 5, 6, 11, 30, 35] {
        if i < n {
            println!("    H[{:2},{:2}] = {:10.4}  ({})", i+1, i+1, h_core[i][i], labels[i]);
        }
    }
    println!();

    // Compute some two-electron integrals
    println!("Computing sample two-electron integrals...");
    let samples = [
        (0, 0, 0, 0),   // C1 1s with itself
        (1, 1, 1, 1),   // C1 2s with itself
        (2, 2, 2, 2),   // C1 2px with itself
        (5, 5, 10, 10), // C1 2pz with C2 2pz
        (30, 30, 31, 31), // H7 1s with H8 1s
    ];

    for &(i, j, k, l) in &samples {
        if i < n && j < n && k < n && l < n {
            let eri = two_electron::electron_repulsion(&basis[i], &basis[j], &basis[k], &basis[l]);
            println!(
                "  ({:2},{:2}|{:2},{:2}) = {:10.6}",
                i + 1,
                j + 1,
                k + 1,
                l + 1,
                eri
            );
        }
    }
    println!();

    // Count unique two-electron integrals
    let n_unique = n * (n + 1) * (n * (n + 1) / 2 + 1) / 8;
    println!("Total unique two-electron integrals: {}", n_unique);
    println!("  (with 8-fold permutational symmetry)");
    println!();

    // Demonstrate parallel computation
    println!("=== Parallel ERI Computation Demo ===");
    println!("Computing {} unique ERIs using parallel Standard method...", n_unique);
    let start = Instant::now();
    let eris = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);
    let parallel_time = start.elapsed();

    println!("  Time: {:.2} ms", parallel_time.as_secs_f64() * 1000.0);
    println!("  Computed: {} integrals", eris.len());
    println!("  Rate: {:.0} integrals/second", eris.len() as f64 / parallel_time.as_secs_f64());
    println!();

    // Summary
    println!("\n=== Summary ===");
    println!("Molecule:          C6H6 (Benzene)");
    println!("Basis set:         STO-3G");
    println!("Basis functions:   {}", n);
    println!("One-e integrals:   {} unique", n * (n + 1) / 2);
    println!("Two-e integrals:   {} unique", n_unique);
    println!();
    println!("Timings:");
    println!("  Overlap matrix S:      {:.2} ms", s_time.as_secs_f64() * 1000.0);
    println!("  Kinetic matrix T:      {:.2} ms", t_time.as_secs_f64() * 1000.0);
    println!("  Nuclear matrix V:      {:.2} ms", v_time.as_secs_f64() * 1000.0);
    println!("  Two-electron (parallel): {:.2} ms", parallel_time.as_secs_f64() * 1000.0);
    println!();
    println!("These integrals enable:");
    println!("  • Hartree-Fock SCF calculations for benzene");
    println!("  • Analysis of aromatic π-electron system");
    println!("  • Post-HF correlation methods");
    println!("  • Excited state calculations");
}
