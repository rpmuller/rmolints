use rmolints::basis::{basis_labels, build_sto3g_basis};
use rmolints::molecule::Molecule;
use rmolints::{one_electron, two_electron};

fn main() {
    println!("\n=== H2O Molecular Integrals (STO-3G Basis) ===\n");

    // Create H2O molecule
    let h2o = Molecule::h2o();
    println!("Molecule: H2O");
    println!("  Atoms: {}", h2o.natoms());
    println!("  Electrons: {}", h2o.nelectrons());
    println!("  Occupied orbitals: {}", h2o.nocc());
    println!();

    // Print geometry
    println!("Geometry (Angstroms):");
    for (i, atom) in h2o.atoms.iter().enumerate() {
        println!(
            "  {:2} {:2}  {:8.4}  {:8.4}  {:8.4}",
            i + 1,
            atom.element.symbol,
            atom.position.x,
            atom.position.y,
            atom.position.z
        );
    }
    println!();

    // Build basis set
    let basis = build_sto3g_basis(&h2o);
    let labels = basis_labels(&basis, &h2o);

    println!("Basis set: STO-3G");
    println!("  Basis functions: {}", basis.len());
    println!();

    println!("Basis function labels:");
    for (i, label) in labels.iter().enumerate() {
        println!("  {:2}: {}", i + 1, label);
    }
    println!();

    // Compute overlap matrix
    println!("Computing overlap matrix S...");
    let n = basis.len();
    let mut s_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let s_ij = one_electron::overlap(&basis[i], &basis[j]);
            s_matrix[i][j] = s_ij;
            s_matrix[j][i] = s_ij;
        }
    }

    println!("\nOverlap Matrix S:");
    print_matrix(&s_matrix, &labels);

    // Compute kinetic energy matrix
    println!("\n\nComputing kinetic energy matrix T...");
    let mut t_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let t_ij = one_electron::kinetic(&basis[i], &basis[j]);
            t_matrix[i][j] = t_ij;
            t_matrix[j][i] = t_ij;
        }
    }

    println!("\nKinetic Energy Matrix T:");
    print_matrix(&t_matrix, &labels);

    // Compute nuclear attraction matrix
    println!("\n\nComputing nuclear attraction matrix V...");
    let mut v_matrix = vec![vec![0.0; n]; n];

    for i in 0..n {
        for j in 0..=i {
            let mut v_ij = 0.0;
            // Sum over all nuclei
            for atom in &h2o.atoms {
                v_ij += one_electron::nuclear_attraction(&basis[i], &basis[j], atom.position)
                    * atom.element.charge();
            }
            v_matrix[i][j] = v_ij;
            v_matrix[j][i] = v_ij;
        }
    }

    println!("\nNuclear Attraction Matrix V:");
    print_matrix(&v_matrix, &labels);

    // Core Hamiltonian
    println!("\n\nCore Hamiltonian H_core = T + V:");
    let mut h_core = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..n {
            h_core[i][j] = t_matrix[i][j] + v_matrix[i][j];
        }
    }
    print_matrix(&h_core, &labels);

    // Compute some sample two-electron integrals
    println!("\n\nSample Two-Electron Repulsion Integrals:");
    println!("Format: (ij|kl) = ⟨ij|1/r₁₂|kl⟩\n");

    let samples = [
        (0, 0, 0, 0),
        (0, 0, 1, 1),
        (1, 1, 1, 1),
        (2, 2, 2, 2),
        (0, 1, 0, 1),
        (0, 1, 2, 3),
    ];

    for &(i, j, k, l) in &samples {
        if i < n && j < n && k < l && l < n {
            let eri = two_electron::electron_repulsion(&basis[i], &basis[j], &basis[k], &basis[l]);
            println!(
                "  ({:2}{:2}|{:2}{:2}) = ({:6}|{:6}) = {:12.8}",
                i + 1,
                j + 1,
                k + 1,
                l + 1,
                labels[i],
                labels[k],
                eri
            );
        }
    }

    // Count unique two-electron integrals
    let n_unique = n * (n + 1) * (n * (n + 1) / 2 + 1) / 8;
    println!("\n\nTotal unique two-electron integrals: {}", n_unique);
    println!("  (exploiting 8-fold permutational symmetry)");

    // Summary
    println!("\n\n=== Summary ===");
    println!("Molecule:          H2O");
    println!("Basis set:         STO-3G");
    println!("Basis functions:   {}", n);
    println!("One-e integrals:   {} unique ({}x{} matrix)", n * (n + 1) / 2, n, n);
    println!("Two-e integrals:   {} unique", n_unique);
    println!("\nMatrices computed:");
    println!("  ✓ Overlap matrix S");
    println!("  ✓ Kinetic energy matrix T");
    println!("  ✓ Nuclear attraction matrix V");
    println!("  ✓ Core Hamiltonian H_core");
    println!("  ✓ Sample two-electron integrals (ERI)");
    println!("\nThese integrals are the foundation for:");
    println!("  • Hartree-Fock calculations");
    println!("  • Density Functional Theory (DFT)");
    println!("  • Post-HF methods (MP2, CI, CC)");
}

fn print_matrix(matrix: &[Vec<f64>], labels: &[String]) {
    let n = matrix.len();

    // Print column headers
    print!("       ");
    for j in 0..n {
        print!(" {:>12}", format!("{}", j + 1));
    }
    println!();

    // Print matrix with row labels
    for (i, row) in matrix.iter().enumerate().take(n) {
        print!("{:2} {:6}", i + 1, labels[i]);
        for val in row.iter().take(n) {
            if val.abs() < 1e-10 {
                print!(" {:>12}", "0.0");
            } else {
                print!(" {:>12.6}", val);
            }
        }
        println!();
    }
}
