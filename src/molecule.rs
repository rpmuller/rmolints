//! Molecule representation for quantum chemistry calculations
//!
//! This module provides data structures for representing molecules,
//! including atoms, coordinates, and nuclear charges.

use crate::common::Vec3;

/// Atomic number and common element data
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Element {
    pub atomic_number: u32,
    pub symbol: &'static str,
}

impl Element {
    pub const H: Element = Element { atomic_number: 1, symbol: "H" };
    pub const C: Element = Element { atomic_number: 6, symbol: "C" };
    pub const N: Element = Element { atomic_number: 7, symbol: "N" };
    pub const O: Element = Element { atomic_number: 8, symbol: "O" };
    pub const F: Element = Element { atomic_number: 9, symbol: "F" };
    pub const S: Element = Element { atomic_number: 16, symbol: "S" };
    pub const CL: Element = Element { atomic_number: 17, symbol: "Cl" };

    pub fn from_atomic_number(z: u32) -> Option<Element> {
        match z {
            1 => Some(Element::H),
            6 => Some(Element::C),
            7 => Some(Element::N),
            8 => Some(Element::O),
            9 => Some(Element::F),
            16 => Some(Element::S),
            17 => Some(Element::CL),
            _ => None,
        }
    }

    pub fn charge(&self) -> f64 {
        self.atomic_number as f64
    }
}

/// An atom with an element type and position
#[derive(Debug, Clone)]
pub struct Atom {
    pub element: Element,
    pub position: Vec3,
}

impl Atom {
    pub fn new(element: Element, position: Vec3) -> Self {
        Atom { element, position }
    }

    pub fn new_xyz(element: Element, x: f64, y: f64, z: f64) -> Self {
        Atom {
            element,
            position: Vec3::new(x, y, z),
        }
    }
}

/// A molecule consisting of atoms
#[derive(Debug, Clone)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub charge: i32,
    pub multiplicity: i32,
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>, charge: i32, multiplicity: i32) -> Self {
        Molecule {
            atoms,
            charge,
            multiplicity,
        }
    }

    /// Number of atoms in the molecule
    pub fn natoms(&self) -> usize {
        self.atoms.len()
    }

    /// Total number of electrons (sum of atomic numbers minus charge)
    pub fn nelectrons(&self) -> i32 {
        let nuclear_charge: i32 = self.atoms.iter()
            .map(|a| a.element.atomic_number as i32)
            .sum();
        nuclear_charge - self.charge
    }

    /// Number of occupied orbitals (electrons / 2 for closed shell)
    pub fn nocc(&self) -> usize {
        (self.nelectrons() / 2) as usize
    }

    /// Create H2 molecule
    pub fn h2(bond_length: f64) -> Self {
        Molecule::new(
            vec![
                Atom::new_xyz(Element::H, 0.0, 0.0, 0.0),
                Atom::new_xyz(Element::H, 0.0, 0.0, bond_length),
            ],
            0,
            1,
        )
    }

    /// Create H2O molecule (equilibrium geometry)
    /// Bond length: 0.9584 Å, angle: 104.45°
    pub fn h2o() -> Self {
        let r = 0.9584; // O-H bond length in Angstroms
        let theta = 104.45_f64.to_radians() / 2.0; // Half angle

        Molecule::new(
            vec![
                Atom::new_xyz(Element::O, 0.0, 0.0, 0.0),
                Atom::new_xyz(Element::H, r * theta.sin(), r * theta.cos(), 0.0),
                Atom::new_xyz(Element::H, -r * theta.sin(), r * theta.cos(), 0.0),
            ],
            0,
            1,
        )
    }

    /// Create NH3 molecule (equilibrium geometry)
    pub fn nh3() -> Self {
        let r = 1.012; // N-H bond length in Angstroms
        let angle = 106.67_f64.to_radians(); // H-N-H angle

        // Pyramidal geometry
        let h = r * (angle / 2.0).cos();
        let w = r * (angle / 2.0).sin();

        Molecule::new(
            vec![
                Atom::new_xyz(Element::N, 0.0, 0.0, 0.0),
                Atom::new_xyz(Element::H, 0.0, h, w),
                Atom::new_xyz(Element::H, w * (120.0_f64.to_radians()).cos(), h, w * (120.0_f64.to_radians()).sin()),
                Atom::new_xyz(Element::H, w * (240.0_f64.to_radians()).cos(), h, w * (240.0_f64.to_radians()).sin()),
            ],
            0,
            1,
        )
    }

    /// Create benzene molecule (C6H6)
    pub fn benzene() -> Self {
        let r_cc = 1.39; // C-C bond length in Angstroms
        let r_ch = 1.08; // C-H bond length in Angstroms
        let mut atoms = Vec::new();

        // Add carbon atoms in hexagon
        for i in 0..6 {
            let angle = i as f64 * 60.0_f64.to_radians();
            atoms.push(Atom::new_xyz(
                Element::C,
                r_cc * angle.cos(),
                r_cc * angle.sin(),
                0.0,
            ));
        }

        // Add hydrogen atoms
        for i in 0..6 {
            let angle = i as f64 * 60.0_f64.to_radians();
            let r_total = r_cc + r_ch;
            atoms.push(Atom::new_xyz(
                Element::H,
                r_total * angle.cos(),
                r_total * angle.sin(),
                0.0,
            ));
        }

        Molecule::new(atoms, 0, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_h2_molecule() {
        let h2 = Molecule::h2(1.4); // 1.4 Bohr
        assert_eq!(h2.natoms(), 2);
        assert_eq!(h2.nelectrons(), 2);
        assert_eq!(h2.nocc(), 1);
    }

    #[test]
    fn test_h2o_molecule() {
        let h2o = Molecule::h2o();
        assert_eq!(h2o.natoms(), 3);
        assert_eq!(h2o.nelectrons(), 10);
        assert_eq!(h2o.nocc(), 5);
    }

    #[test]
    fn test_benzene_molecule() {
        let benzene = Molecule::benzene();
        assert_eq!(benzene.natoms(), 12); // 6 C + 6 H
        assert_eq!(benzene.nelectrons(), 42); // 6*6 + 6*1
        assert_eq!(benzene.nocc(), 21);
    }
}
