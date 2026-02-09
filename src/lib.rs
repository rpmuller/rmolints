//! # rmolints
//!
//! A high-performance molecular integrals library written in Rust.

#![cfg_attr(feature = "simd", feature(portable_simd))]
//!
//! This library provides implementations of one- and two-electron integrals
//! used in quantum chemistry calculations, ported from the PyQuante2 project.
//!
//! ## Modules
//!
//! - `one_electron`: One-electron integrals (overlap, kinetic energy, nuclear attraction)
//! - `two_electron`: Two-electron repulsion integrals (standard implementation)
//! - `rys`: Two-electron integrals using Rys quadrature
//! - `hgp`: Two-electron integrals using Head-Gordon and Pople method
//! - `parallel`: Parallel computation of integrals using Rayon
//! - `boys`: Boys function for computing auxiliary integrals

pub mod one_electron;
pub mod two_electron;
pub mod rys;
pub mod hgp;
#[cfg(feature = "simd")]
pub mod hgp_simd;
pub mod parallel;
pub mod boys;
pub mod utils;
pub mod molecule;
pub mod basis;
pub mod basis_data;

/// Common types and utilities used across modules
pub mod common {
    /// 3D point or vector
    #[derive(Debug, Clone, Copy, PartialEq)]
    pub struct Vec3 {
        pub x: f64,
        pub y: f64,
        pub z: f64,
    }

    impl Vec3 {
        pub fn new(x: f64, y: f64, z: f64) -> Self {
            Self { x, y, z }
        }

        pub fn distance_squared(&self, other: &Vec3) -> f64 {
            let dx = self.x - other.x;
            let dy = self.y - other.y;
            let dz = self.z - other.z;
            dx * dx + dy * dy + dz * dz
        }

        pub fn distance(&self, other: &Vec3) -> f64 {
            self.distance_squared(other).sqrt()
        }
    }

    /// Gaussian basis function primitive
    #[derive(Debug, Clone, Copy)]
    pub struct Primitive {
        pub exponent: f64,
        pub coefficient: f64,
    }

    /// Contracted Gaussian basis function
    #[derive(Debug, Clone)]
    pub struct CGBF {
        pub origin: Vec3,
        pub shell: (i32, i32, i32), // (l, m, n) - angular momentum quantum numbers
        pub primitives: Vec<Primitive>,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vec3_distance() {
        let a = common::Vec3::new(0.0, 0.0, 0.0);
        let b = common::Vec3::new(1.0, 0.0, 0.0);
        assert!((a.distance(&b) - 1.0).abs() < 1e-10);
    }
}
