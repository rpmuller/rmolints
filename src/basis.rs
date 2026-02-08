//! Basis set definitions and utilities
//!
//! This module provides standard basis sets (STO-3G, etc.) and utilities
//! for creating basis functions for molecules.

use crate::common::{CGBF, Primitive, Vec3};
use crate::molecule::{Atom, Element, Molecule};

/// STO-3G basis set data
/// Each entry is (exponent, coefficient) for the 3 Gaussian primitives

/// STO-3G for Hydrogen (1s)
const H_STO3G_1S: [(f64, f64); 3] = [
    (3.425250914, 0.154328967),
    (0.623913730, 0.535328142),
    (0.168855404, 0.444634542),
];

/// STO-3G for Carbon
/// 1s orbital
const C_STO3G_1S: [(f64, f64); 3] = [
    (71.616837, 0.154328967),
    (13.045096, 0.535328142),
    (3.5305122, 0.444634542),
];

/// 2s orbital
const C_STO3G_2S: [(f64, f64); 3] = [
    (2.9412494, -0.099967230),
    (0.6834831, 0.399512826),
    (0.2222899, 0.700115469),
];

/// 2p orbitals
const C_STO3G_2P: [(f64, f64); 3] = [
    (2.9412494, 0.155916275),
    (0.6834831, 0.607683719),
    (0.2222899, 0.391957393),
];

/// STO-3G for Nitrogen
const N_STO3G_1S: [(f64, f64); 3] = [
    (99.106169, 0.154328967),
    (18.052312, 0.535328142),
    (4.8856602, 0.444634542),
];

const N_STO3G_2S: [(f64, f64); 3] = [
    (3.7804559, -0.099967230),
    (0.8784966, 0.399512826),
    (0.2857144, 0.700115469),
];

const N_STO3G_2P: [(f64, f64); 3] = [
    (3.7804559, 0.155916275),
    (0.8784966, 0.607683719),
    (0.2857144, 0.391957393),
];

/// STO-3G for Oxygen
const O_STO3G_1S: [(f64, f64); 3] = [
    (130.70932, 0.154328967),
    (23.808861, 0.535328142),
    (6.4436083, 0.444634542),
];

const O_STO3G_2S: [(f64, f64); 3] = [
    (5.0331513, -0.099967230),
    (1.1695961, 0.399512826),
    (0.3803890, 0.700115469),
];

const O_STO3G_2P: [(f64, f64); 3] = [
    (5.0331513, 0.155916275),
    (1.1695961, 0.607683719),
    (0.3803890, 0.391957393),
];

/// Create basis functions for an atom using STO-3G basis
pub fn sto3g_basis(atom: &Atom) -> Vec<CGBF> {
    let pos = atom.position;
    let mut basis = Vec::new();

    match atom.element {
        Element::H => {
            // 1s orbital
            basis.push(make_cgbf(pos, (0, 0, 0), &H_STO3G_1S));
        }
        Element::C => {
            // 1s orbital
            basis.push(make_cgbf(pos, (0, 0, 0), &C_STO3G_1S));
            // 2s orbital
            basis.push(make_cgbf(pos, (0, 0, 0), &C_STO3G_2S));
            // 2p orbitals (px, py, pz)
            basis.push(make_cgbf(pos, (1, 0, 0), &C_STO3G_2P));
            basis.push(make_cgbf(pos, (0, 1, 0), &C_STO3G_2P));
            basis.push(make_cgbf(pos, (0, 0, 1), &C_STO3G_2P));
        }
        Element::N => {
            basis.push(make_cgbf(pos, (0, 0, 0), &N_STO3G_1S));
            basis.push(make_cgbf(pos, (0, 0, 0), &N_STO3G_2S));
            basis.push(make_cgbf(pos, (1, 0, 0), &N_STO3G_2P));
            basis.push(make_cgbf(pos, (0, 1, 0), &N_STO3G_2P));
            basis.push(make_cgbf(pos, (0, 0, 1), &N_STO3G_2P));
        }
        Element::O => {
            basis.push(make_cgbf(pos, (0, 0, 0), &O_STO3G_1S));
            basis.push(make_cgbf(pos, (0, 0, 0), &O_STO3G_2S));
            basis.push(make_cgbf(pos, (1, 0, 0), &O_STO3G_2P));
            basis.push(make_cgbf(pos, (0, 1, 0), &O_STO3G_2P));
            basis.push(make_cgbf(pos, (0, 0, 1), &O_STO3G_2P));
        }
        _ => panic!("STO-3G basis not implemented for element {:?}", atom.element),
    }

    basis
}

/// Helper to create a CGBF from primitives data
fn make_cgbf(origin: Vec3, shell: (i32, i32, i32), data: &[(f64, f64)]) -> CGBF {
    let primitives = data
        .iter()
        .map(|&(exp, coef)| Primitive {
            exponent: exp,
            coefficient: coef,
        })
        .collect();

    CGBF {
        origin,
        shell,
        primitives,
    }
}

/// Build complete STO-3G basis set for a molecule
pub fn build_sto3g_basis(molecule: &Molecule) -> Vec<CGBF> {
    let mut basis = Vec::new();

    for atom in &molecule.atoms {
        basis.extend(sto3g_basis(atom));
    }

    basis
}

/// Get basis function labels for printing
pub fn basis_labels(_basis: &[CGBF], molecule: &Molecule) -> Vec<String> {
    let mut labels = Vec::new();

    for (idx, atom) in molecule.atoms.iter().enumerate() {
        let atom_basis = sto3g_basis(atom);
        let mut local_count = 0;

        for cgbf in &atom_basis {
            let orbital_type = match cgbf.shell {
                (0, 0, 0) => {
                    if atom.element == Element::H {
                        "1s"
                    } else if local_count == 0 {
                        "1s"
                    } else {
                        "2s"
                    }
                }
                (1, 0, 0) => "2px",
                (0, 1, 0) => "2py",
                (0, 0, 1) => "2pz",
                _ => "unknown",
            };

            labels.push(format!("{}{} {}", atom.element.symbol, idx + 1, orbital_type));
            local_count += 1;
        }
    }

    labels
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::molecule::Molecule;

    #[test]
    fn test_h2_basis() {
        let h2 = Molecule::h2(1.4);
        let basis = build_sto3g_basis(&h2);
        assert_eq!(basis.len(), 2); // 2 1s orbitals
    }

    #[test]
    fn test_h2o_basis() {
        let h2o = Molecule::h2o();
        let basis = build_sto3g_basis(&h2o);
        // O: 1s, 2s, 2px, 2py, 2pz (5)
        // H: 1s (1 each, 2 total)
        // Total: 7 basis functions
        assert_eq!(basis.len(), 7);
    }

    #[test]
    fn test_benzene_basis() {
        let benzene = Molecule::benzene();
        let basis = build_sto3g_basis(&benzene);
        // 6 C atoms: 5 functions each = 30
        // 6 H atoms: 1 function each = 6
        // Total: 36 basis functions
        assert_eq!(basis.len(), 36);
    }
}
