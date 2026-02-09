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
// POLARIZATION EXPONENTS
// ============================================================================

/// D polarization exponents for 6-31G* and 6-31G**
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

/// P polarization exponents for 6-31G**  (aka 6-31G(d,p))
/// Returns (p_exponent, number_of_p_functions)
pub const fn get_p_polarization(atomic_number: u32) -> Option<(f64, usize)> {
    match atomic_number {
        1 => Some((1.1, 3)),  // H: px, py, pz
        2 => Some((1.1, 3)),  // He: px, py, pz
        _ => None, // Only H and He get p polarization functions
    }
}
// STO-3G Basis Sets

pub const LI_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[16.119575, 2.9362007, 0.7946505],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.6362897, 0.1478601, 0.0480887],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],  // S coefficients
            &[0.15591627, 0.60768372, 0.39195739],  // P coefficients
        ],
    },
];


pub const BE_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[30.167871, 5.4951153, 1.4871927],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.3148331, 0.3055389, 0.0993707],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],  // S coefficients
            &[0.15591627, 0.60768372, 0.39195739],  // P coefficients
        ],
    },
];


pub const B_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[48.791113, 8.8873622, 2.405267],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.2369561, 0.5198205, 0.1690618],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],  // S coefficients
            &[0.15591627, 0.60768372, 0.39195739],  // P coefficients
        ],
    },
];


pub const F_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[166.67913, 30.360812, 8.2168207],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[6.4648032, 1.5022812, 0.4885885],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],  // S coefficients
            &[0.15591627, 0.60768372, 0.39195739],  // P coefficients
        ],
    },
];


pub const NE_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[207.01561, 37.708151, 10.205297],
        coefficients: &[&[0.15432897, 0.53532814, 0.44463454]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[8.2463151, 1.9162662, 0.6232293],
        coefficients: &[
            &[-0.09996723, 0.39951283, 0.70011547],  // S coefficients
            &[0.15591627, 0.60768372, 0.39195739],  // P coefficients
        ],
    },
];


pub const NA_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[250.77243, 45.678511, 12.362388],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[12.040193, 2.7978819, 0.909958],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.4787406, 0.4125649, 0.1614751],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const MG_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[299.2374, 54.50647, 14.75158],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[15.12182, 3.513987, 1.142857],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.395448, 0.389326, 0.15238],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const AL_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[351.4214767, 64.01186067, 17.32410761],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[18.89939621, 4.391813233, 1.42835397],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.395448293, 0.3893265318, 0.1523797659],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const SI_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[407.7975514, 74.28083305, 20.10329229],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[23.19365606, 5.389706871, 1.752899952],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.478740622, 0.4125648801, 0.1614750979],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const P_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[468.3656378, 85.31338559, 23.08913156],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[28.03263958, 6.514182577, 2.118614352],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.743103231, 0.4863213771, 0.1903428909],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const S_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[533.1257359, 97.1095183, 26.28162542],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[33.32975173, 7.745117521, 2.518952599],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.029194274, 0.5661400518, 0.2215833792],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const CL_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[601.3456136, 109.5358542, 29.64467686],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[38.96041889, 9.053563477, 2.944499834],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.129386495, 0.5940934274, 0.232524141],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


pub const AR_STO3G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[674.4465184, 122.8512753, 33.24834945],
        coefficients: &[&[0.1543289673, 0.5353281423, 0.4446345422]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[45.16424392, 10.495199, 3.413364448],
        coefficients: &[
            &[-0.09996722919, 0.3995128261, 0.7001154689],  // S coefficients
            &[0.155916275, 0.6076837186, 0.3919573931],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.621366518, 0.731354605, 0.2862472356],
        coefficients: &[
            &[-0.219620369, 0.2255954336, 0.900398426],  // S coefficients
            &[0.01058760429, 0.5951670053, 0.462001012],  // P coefficients
        ],
    },
];


// 6-31G Basis Sets

pub const HE_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[38.421634, 5.77803, 1.241774],
        coefficients: &[&[0.023766, 0.154679, 0.46963]],
    },
    ShellData {
        shell_type: ShellType::S,
        exponents: &[0.297964],
        coefficients: &[&[1.0]],
    },
];


pub const LI_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[642.41892, 96.798515, 22.091121, 6.2010703, 1.9351177, 0.6367358],
        coefficients: &[&[0.0021426, 0.0162089, 0.0773156, 0.245786, 0.470189, 0.3454708]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.3249184, 0.6324306, 0.0790534],
        coefficients: &[
            &[-0.0350917, -0.1912328, 1.0839878],  // S coefficients
            &[0.0089415, 0.1410095, 0.9453637],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.035962],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const BE_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[1264.5857, 189.93681, 43.159089, 12.098663, 3.8063232, 1.2728903],
        coefficients: &[&[0.0019448, 0.0148351, 0.0720906, 0.2371542, 0.4691987, 0.3565202]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[3.1964631, 0.7478133, 0.2199663],
        coefficients: &[
            &[-0.1126487, -0.2295064, 1.1869167],  // S coefficients
            &[0.0559802, 0.2615506, 0.7939723],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.0823099],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const B_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[2068.8823, 310.64957, 70.683033, 19.86108, 6.2993048, 2.127027],
        coefficients: &[&[0.0018663, 0.0142515, 0.0695516, 0.2325729, 0.4670787, 0.3634314]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[4.727971, 1.1903377, 0.3594117],
        coefficients: &[
            &[-0.1303938, -0.1307889, 1.1309444],  // S coefficients
            &[0.0745976, 0.3078467, 0.7434568],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.1267512],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const N_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[4173.511, 627.4579, 142.9021, 40.23433, 12.82021, 4.390437],
        coefficients: &[&[0.0018348, 0.013995, 0.068587, 0.232241, 0.46907, 0.360455]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[11.626358, 2.71628, 0.772218],
        coefficients: &[
            &[-0.114961, -0.169118, 1.145852],  // S coefficients
            &[0.06758, 0.323907, 0.740895],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.2120313],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const F_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[7001.71309, 1051.36609, 239.28569, 67.3974453, 21.5199573, 7.4031013],
        coefficients: &[&[0.0018196169, 0.0139160796, 0.0684053245, 0.23318576, 0.471267439, 0.356618546]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[20.8479528, 4.80830834, 1.34406986],
        coefficients: &[
            &[-0.108506975, -0.146451658, 1.12868858],  // S coefficients
            &[0.0716287243, 0.345912103, 0.722469957],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.358151393],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const NE_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[8425.85153, 1268.5194, 289.621414, 81.859004, 26.2515079, 9.09472051],
        coefficients: &[&[0.0018843481, 0.0143368994, 0.0701096233, 0.237373266, 0.473007126, 0.348401241]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[26.532131, 6.10175501, 1.69627153],
        coefficients: &[
            &[-0.107118287, -0.146163821, 1.1277735],  // S coefficients
            &[0.0719095885, 0.349513372, 0.719940512],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.4458187],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const NA_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[9993.2, 1499.89, 341.951, 94.6797, 29.7345, 10.0063],
        coefficients: &[&[0.0019377, 0.014807, 0.072706, 0.252629, 0.493242, 0.313169]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[150.963, 35.5878, 11.1683, 3.90201, 1.38177, 0.466382],
        coefficients: &[
            &[-0.0035421, -0.043959, -0.1097521, 0.187398, 0.646699, 0.306058],  // S coefficients
            &[0.0050017, 0.035511, 0.142825, 0.33862, 0.451579, 0.273271],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.497966, 0.084353, 0.066635],
        coefficients: &[
            &[-0.248503, -0.131704, 1.23352],  // S coefficients
            &[-0.023023, 0.950359, 0.059858],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.0259544],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const MG_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[11722.8, 1759.93, 400.846, 112.807, 35.9997, 12.1828],
        coefficients: &[&[0.0019778, 0.015114, 0.073911, 0.249191, 0.487928, 0.319662]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[189.18, 45.2119, 14.3563, 5.13886, 1.90652, 0.705887],
        coefficients: &[
            &[-0.0032372, -0.041008, -0.1126, 0.148633, 0.616497, 0.364829],  // S coefficients
            &[0.0049281, 0.034989, 0.140725, 0.333642, 0.44494, 0.269254],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.92934, 0.269035, 0.117379],
        coefficients: &[
            &[-0.21229, -0.107985, 1.17584],  // S coefficients
            &[-0.022419, 0.19227, 0.846181],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.0421061],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const AL_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[13983.1, 2098.75, 477.705, 134.36, 42.8709, 14.5189],
        coefficients: &[&[0.00194267, 0.0148599, 0.0728494, 0.24683, 0.487258, 0.323496]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[239.668, 57.4419, 18.2859, 6.59914, 2.49049, 0.94454],
        coefficients: &[
            &[-0.00292619, -0.037408, -0.114487, 0.115635, 0.612595, 0.393799],  // S coefficients
            &[0.00460285, 0.033199, 0.136282, 0.330476, 0.449146, 0.265704],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.2779, 0.39759, 0.160095],
        coefficients: &[
            &[-0.227606, 0.00144583, 1.09279],  // S coefficients
            &[-0.017513, 0.244533, 0.804934],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.0556577],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const SI_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[16115.9, 2425.58, 553.867, 156.34, 50.0683, 17.0178],
        coefficients: &[&[0.00195948, 0.0149288, 0.0728478, 0.24613, 0.485914, 0.325002]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[292.718, 69.8731, 22.3363, 8.15039, 3.13458, 1.22543],
        coefficients: &[
            &[-0.00278094, -0.0357146, -0.114985, 0.0935634, 0.603017, 0.418959],  // S coefficients
            &[0.00443826, 0.0326679, 0.134721, 0.328678, 0.44964, 0.261372],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[1.72738, 0.572922, 0.222192],
        coefficients: &[
            &[-0.24463, 0.00431572, 1.09818],  // S coefficients
            &[-0.0177951, 0.253539, 0.800669],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.0778369],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const P_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[19413.3, 2909.42, 661.364, 185.759, 59.1943, 20.031],
        coefficients: &[&[0.0018516, 0.0142062, 0.0699995, 0.240079, 0.484762, 0.3352]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[339.478, 81.0101, 25.878, 9.45221, 3.66566, 1.46746],
        coefficients: &[
            &[-0.00278217, -0.0360499, -0.116631, 0.0968328, 0.614418, 0.403798],  // S coefficients
            &[0.00456462, 0.0336936, 0.139755, 0.339362, 0.450921, 0.238586],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.15623, 0.748997, 0.283145],
        coefficients: &[
            &[-0.252923, 0.0328517, 1.08125],  // S coefficients
            &[-0.0177653, 0.274058, 0.785421],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.0998317],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const S_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[21917.1, 3301.49, 754.146, 212.711, 67.9896, 23.0515],
        coefficients: &[&[0.001869, 0.01423, 0.069696, 0.238487, 0.483307, 0.338074]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[423.735, 100.71, 32.1599, 11.8079, 4.6311, 1.87025],
        coefficients: &[
            &[-0.0023767, -0.031693, -0.113317, 0.05609, 0.592255, 0.455006],  // S coefficients
            &[0.004061, 0.030681, 0.130452, 0.327205, 0.452851, 0.256042],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[2.61584, 0.922167, 0.341287],
        coefficients: &[
            &[-0.250374, 0.066957, 1.05451],  // S coefficients
            &[-0.014511, 0.310263, 0.754483],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.117167],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const CL_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[25180.1, 3780.35, 860.474, 242.145, 77.3349, 26.247],
        coefficients: &[&[0.001833, 0.014034, 0.069097, 0.237452, 0.483034, 0.339856]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[491.765, 116.984, 37.4153, 13.7834, 5.45215, 2.22588],
        coefficients: &[
            &[-0.0022974, -0.030714, -0.112528, 0.045016, 0.589353, 0.465206],  // S coefficients
            &[0.0039894, 0.030318, 0.12988, 0.327951, 0.453527, 0.252154],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[3.18649, 1.14427, 0.420377],
        coefficients: &[
            &[-0.25183, 0.061589, 1.06018],  // S coefficients
            &[-0.014299, 0.323572, 0.743507],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.142657],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


pub const AR_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[28348.3, 4257.62, 969.857, 273.263, 87.3695, 29.6867],
        coefficients: &[&[0.00182526, 0.0139686, 0.0687073, 0.236204, 0.482214, 0.342043]],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[575.891, 136.816, 43.8098, 16.2094, 6.46084, 2.65114],
        coefficients: &[
            &[-0.00215972, -0.0290775, -0.110827, 0.0276999, 0.577613, 0.488688],  // S coefficients
            &[0.00380665, 0.0292305, 0.126467, 0.32351, 0.454896, 0.25663],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[3.86028, 1.41373, 0.516646],
        coefficients: &[
            &[-0.255592, 0.0378066, 1.08056],  // S coefficients
            &[-0.0159197, 0.324646, 0.74399],  // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[0.173888],
        coefficients: &[
            &[1.0],  // S coefficients
            &[1.0],  // P coefficients
        ],
    },
];


