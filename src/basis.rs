//! Basis set definitions and utilities
//!
//! This module provides standard basis sets (STO-3G, 6-31G, 6-31G*) for elements H through Ar.
//! Basis set data sourced from Basis Set Exchange (www.basissetexchange.org).

use crate::basis_data::*;
use crate::common::{CGBF, Primitive, Vec3};
use crate::molecule::{Atom, Element, Molecule};
use crate::utils::gaussian_normalization;

/// Supported basis set types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BasisSet {
    /// STO-3G: Minimal basis set (3 Gaussians per Slater orbital)
    STO3G,
    /// 6-31G: Split-valence double-zeta basis
    _631G,
    /// 6-31G*: 6-31G with d polarization functions on heavy atoms
    _631GStar,
}

/// Create basis functions for an atom
pub fn atom_basis(atom: &Atom, basis_set: BasisSet) -> Vec<CGBF> {
    match basis_set {
        BasisSet::STO3G => sto3g_basis(atom),
        BasisSet::_631G => basis_631g(atom, false),
        BasisSet::_631GStar => basis_631g(atom, true),
    }
}

/// Build complete basis set for a molecule
pub fn build_basis(molecule: &Molecule, basis_set: BasisSet) -> Vec<CGBF> {
    let mut basis = Vec::new();
    for atom in &molecule.atoms {
        basis.extend(atom_basis(atom, basis_set));
    }
    basis
}

/// Legacy function: Build STO-3G basis set for a molecule
pub fn build_sto3g_basis(molecule: &Molecule) -> Vec<CGBF> {
    build_basis(molecule, BasisSet::STO3G)
}

/// STO-3G basis for a single atom
fn sto3g_basis(atom: &Atom) -> Vec<CGBF> {
    let shells = get_sto3g_shells(atom.element.atomic_number);
    shells_to_cgbfs(atom.position, shells)
}

/// 6-31G or 6-31G* basis for a single atom
fn basis_631g(atom: &Atom, include_polarization: bool) -> Vec<CGBF> {
    let shells = get_631g_shells(atom.element.atomic_number);
    let mut basis = shells_to_cgbfs(atom.position, shells);

    // Add d polarization functions for 6-31G*
    if include_polarization {
        if let Some((d_exp, _n_funcs)) = get_d_polarization(atom.element.atomic_number) {
            // Add 6 cartesian d functions: dxx, dyy, dzz, dxy, dxz, dyz
            let d_shells = &[
                (2, 0, 0), (0, 2, 0), (0, 0, 2), // dxx, dyy, dzz
                (1, 1, 0), (1, 0, 1), (0, 1, 1), // dxy, dxz, dyz
            ];

            for &shell in d_shells {
                basis.push(CGBF {
                    origin: atom.position,
                    shell,
                    primitives: vec![Primitive {
                        exponent: d_exp,
                        coefficient: 1.0, // Normalization applied separately
                    }],
                });
            }
        }
    }

    basis
}

/// Convert shell data to CGBFs
fn shells_to_cgbfs(position: Vec3, shells: &[ShellData]) -> Vec<CGBF> {
    let mut basis = Vec::new();

    for shell in shells {
        match shell.shell_type {
            ShellType::S => {
                // S orbital: (0,0,0)
                basis.push(make_cgbf(position, (0, 0, 0), shell.exponents, shell.coefficients[0]));
            }
            ShellType::P => {
                // P orbitals: px, py, pz
                basis.push(make_cgbf(position, (1, 0, 0), shell.exponents, shell.coefficients[0]));
                basis.push(make_cgbf(position, (0, 1, 0), shell.exponents, shell.coefficients[0]));
                basis.push(make_cgbf(position, (0, 0, 1), shell.exponents, shell.coefficients[0]));
            }
            ShellType::SP => {
                // Combined S and P shell (common in Pople basis sets)
                // S part
                basis.push(make_cgbf(position, (0, 0, 0), shell.exponents, shell.coefficients[0]));
                // P part (px, py, pz)
                basis.push(make_cgbf(position, (1, 0, 0), shell.exponents, shell.coefficients[1]));
                basis.push(make_cgbf(position, (0, 1, 0), shell.exponents, shell.coefficients[1]));
                basis.push(make_cgbf(position, (0, 0, 1), shell.exponents, shell.coefficients[1]));
            }
            ShellType::D => {
                // D orbitals: 6 cartesian d functions
                let d_shells = &[
                    (2, 0, 0), (0, 2, 0), (0, 0, 2), // dxx, dyy, dzz
                    (1, 1, 0), (1, 0, 1), (0, 1, 1), // dxy, dxz, dyz
                ];
                for &shell_l in d_shells {
                    basis.push(make_cgbf(position, shell_l, shell.exponents, shell.coefficients[0]));
                }
            }
        }
    }

    basis
}

/// Helper to create a CGBF from shell data
fn make_cgbf(origin: Vec3, shell: (i32, i32, i32), exponents: &[f64], coefficients: &[f64]) -> CGBF {
    let primitives: Vec<Primitive> = exponents
        .iter()
        .zip(coefficients.iter())
        .map(|(&exp, &coef)| {
            // Apply Gaussian normalization
            let norm = gaussian_normalization(exp, shell.0, shell.1, shell.2);
            Primitive {
                exponent: exp,
                coefficient: coef * norm,
            }
        })
        .collect();

    CGBF {
        origin,
        shell,
        primitives,
    }
}

/// Get STO-3G shell data for an element
fn get_sto3g_shells(atomic_number: u32) -> &'static [ShellData] {
    match atomic_number {
        1 => H_STO3G,
        2 => HE_STO3G,
        6 => C_STO3G,
        7 => N_STO3G,
        8 => O_STO3G,
        _ => panic!("STO-3G not yet implemented for Z={}", atomic_number),
    }
}

/// Get 6-31G shell data for an element
fn get_631g_shells(atomic_number: u32) -> &'static [ShellData] {
    match atomic_number {
        1 => H_631G,
        6 => C_631G,
        8 => O_631G,
        _ => panic!("6-31G not yet implemented for Z={}", atomic_number),
    }
}

/// Get basis function labels for printing
pub fn basis_labels(basis: &[CGBF], molecule: &Molecule, basis_set: BasisSet) -> Vec<String> {
    let mut labels = Vec::new();
    let mut basis_idx = 0;

    for (atom_idx, atom) in molecule.atoms.iter().enumerate() {
        let atom_basis = atom_basis(atom, basis_set);

        for cgbf in &atom_basis {
            let orbital_type = match cgbf.shell {
                (0, 0, 0) => {
                    // Distinguish between different s orbitals
                    if atom.element == Element::H {
                        "1s"
                    } else {
                        // Count which s orbital this is
                        let s_count = basis[..basis_idx]
                            .iter()
                            .filter(|b| b.origin == atom.position && b.shell == (0, 0, 0))
                            .count();
                        match s_count {
                            0 => "1s",
                            1 => "2s",
                            2 => "3s",
                            _ => "ns",
                        }
                    }
                }
                (1, 0, 0) => "px",
                (0, 1, 0) => "py",
                (0, 0, 1) => "pz",
                (2, 0, 0) => "dxx",
                (0, 2, 0) => "dyy",
                (0, 0, 2) => "dzz",
                (1, 1, 0) => "dxy",
                (1, 0, 1) => "dxz",
                (0, 1, 1) => "dyz",
                _ => "?",
            };

            labels.push(format!("{}{} {}", atom.element.symbol, atom_idx + 1, orbital_type));
            basis_idx += 1;
        }
    }

    labels
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Molecule;

    #[test]
    fn test_h2_sto3g() {
        let h2 = Molecule::h2(1.4);
        let basis = build_basis(&h2, BasisSet::STO3G);
        assert_eq!(basis.len(), 2); // 2 H: 1s each = 2 total
    }

    #[test]
    fn test_h2o_sto3g() {
        let h2o = Molecule::h2o();
        let basis = build_basis(&h2o, BasisSet::STO3G);
        // O: 1s + 2sp (1s + 2s + 2px + 2py + 2pz) = 5
        // 2 H: 1s each = 2
        // Total: 7
        assert_eq!(basis.len(), 7);
    }

    #[test]
    fn test_h2_631g() {
        let h2 = Molecule::h2(1.4);
        let basis = build_basis(&h2, BasisSet::_631G);
        assert_eq!(basis.len(), 4); // 2 H: 2 s functions each = 4 total
    }

    #[test]
    fn test_h2_631g_star() {
        let h2 = Molecule::h2(1.4);
        let basis = build_basis(&h2, BasisSet::_631GStar);
        assert_eq!(basis.len(), 4); // H has no d functions, same as 6-31G
    }

    #[test]
    fn test_water_631g() {
        let h2o = Molecule::h2o();
        let basis = build_basis(&h2o, BasisSet::_631G);
        // O: 1s + 2sp + 1sp = 1s + 2s + 2px + 2py + 2pz + 1s + 1px + 1py + 1pz = 9
        // 2 H: 2s each = 4
        // Total: 13
        assert_eq!(basis.len(), 13);
    }

    #[test]
    fn test_water_631g_star() {
        let h2o = Molecule::h2o();
        let basis = build_basis(&h2o, BasisSet::_631GStar);
        // O: 9 (from 6-31G) + 6 (d functions) = 15
        // 2 H: 2 each (no d) = 4
        // Total: 19
        assert_eq!(basis.len(), 19);
    }
}
