//! Basis set coefficient data for elements H through Ar
//!
//! Data sourced from Basis Set Exchange (www.basissetexchange.org)
//! References:
//! - STO-3G: Hehre, Stewart, Pople, J. Chem. Phys. 56, 2657 (1969)
//! - 6-31G: Hehre, Ditchfield, Pople, J. Chem. Phys. 56, 2257 (1972)
//! - 6-31G*: Hariharan and Pople, Theoret. Chimica Acta 28, 213 (1973)

/// Shell data: exponents and coefficients for a single shell
#[derive(Debug, Clone)]
pub struct ShellData {
    pub shell_type: ShellType,
    pub exponents: &'static [f64],
    pub coefficients: &'static [&'static [f64]], // Multiple sets for SP shells
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ShellType {
    S,  // s orbital
    P,  // p orbital
    SP, // Combined s and p (same exponents, different coefficients)
    D,  // d orbital (polarization)
}

// ============================================================================
// STO-3G BASIS SET DATA
// ============================================================================

// Hydrogen (Z=1)
pub const H_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[3.42525091, 0.62391373, 0.1688554],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
];

// Helium (Z=2)
pub const HE_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[6.36242139, 1.158923, 0.31364979],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
];

// Carbon (Z=6) - keeping existing data
pub const C_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[71.616837, 13.045096, 3.5305122],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.9412494, 0.6834831, 0.2222899],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547], // S coefficients
            &[0.15591628, 0.60768372, 0.39195739],  // P coefficients
        ],
    },
];

// Nitrogen (Z=7)
pub const N_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[99.106169, 18.052312, 4.8856602],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[3.7804559, 0.8784966, 0.2857144],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],
            &[0.15591628, 0.60768372, 0.39195739],
        ],
    },
];

// Oxygen (Z=8)
pub const O_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[130.70932, 23.808861, 6.4436083],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[5.0331513, 1.1695961, 0.380389],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],
            &[0.15591628, 0.60768372, 0.39195739],
        ],
    },
];

// ============================================================================
// 6-31G BASIS SET DATA (Selected elements - expand as needed)
// ============================================================================

// Hydrogen 6-31G
pub const H_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[18.731137, 2.8253937, 0.6401217],
        coefficients: &[&[0.03349460, 0.23472695, 0.81375733]],
    },
    ShellData {
        shell_type: ShellType::S,
        exponents: &[0.1612778],
        coefficients: &[&[1.0]],
    },
];

// Carbon 6-31G
pub const C_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[3047.5249, 457.36951, 103.94869, 29.210155, 9.2866630, 3.1639270],
        coefficients: &[&[0.0018347, 0.0140373, 0.0688426, 0.2321844, 0.4679413, 0.3623120]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[7.8682724, 1.8812885, 0.5442493],
        coefficients: &[
            &[-0.1193324, 0.1608542, 1.1434564],  // S
            &[0.0689991, 0.3164240, 0.7443083],   // P
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.1687144],
        coefficients: &[
            &[1.0],  // S
            &[1.0],  // P
        ],
    },
];

// Oxygen 6-31G
pub const O_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[5484.6717, 825.23495, 188.04696, 52.964521, 16.897570, 5.7996353],
        coefficients: &[&[0.0018311, 0.0139501, 0.0684451, 0.2327143, 0.4701930, 0.3585209]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[15.539616, 3.5999336, 1.0137618],
        coefficients: &[
            &[-0.1107775, 0.1480263, 1.1307670],
            &[0.0708743, 0.3397528, 0.7271586],
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.2700058],
        coefficients: &[
            &[1.0],
            &[1.0],
        ],
    },
];

// ============================================================================
// 6-31G* POLARIZATION EXPONENTS
// ============================================================================

/// D polarization exponents for 6-31G*
/// Returns (d_exponent, number_of_d_functions)
pub const fn get_d_polarization(atomic_number: u32) -> Option<(f64, usize)> {
    match atomic_number {
        1 | 2 => None, // H, He: no d functions
        3 => Some((0.2, 6)),      // Li
        4 => Some((0.4, 6)),      // Be
        5 => Some((0.6, 6)),      // B
        6 | 7 | 8 | 9 => Some((0.8, 6)),  // C, N, O, F
        10 => Some((0.8, 6)),     // Ne
        11 | 12 => Some((0.175, 6)),      // Na, Mg
        13 => Some((0.325, 6)),   // Al
        14 => Some((0.45, 6)),    // Si
        15 => Some((0.55, 6)),    // P
        16 => Some((0.65, 6)),    // S
        17 => Some((0.75, 6)),    // Cl
        18 => Some((0.85, 6)),    // Ar
        _ => None,
    }
}
