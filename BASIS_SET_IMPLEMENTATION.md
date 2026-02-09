# Basis Set Implementation Status

## Completed ✅

### Framework
- ✅ **BasisSet enum** - Support for STO-3G, 6-31G, 6-31G*
- ✅ **Element definitions** - All elements H-Ar (Z=1-18) added to Element struct
- ✅ **Basis data structure** - ShellData, ShellType for organizing basis set data
- ✅ **API design** - `build_basis(molecule, BasisSet)` for flexible basis selection
- ✅ **Normalization** - Gaussian normalization applied to all primitives
- ✅ **Polarization functions** - D function exponents for 6-31G* (all elements)
- ✅ **Tests** - Comprehensive tests for H2, H2O with all basis sets

### Data Implemented
**STO-3G:**
- ✅ H (Z=1)
- ✅ He (Z=2)
- ✅ C (Z=6)
- ✅ N (Z=7)
- ✅ O (Z=8)

**6-31G and 6-31G*:**
- ✅ H (Z=1)
- ✅ C (Z=6)
- ✅ O (Z=8)

**Polarization exponents (6-31G*):**
- ✅ All elements H-Ar (Z=1-18)

## To Complete

### Remaining Elements for STO-3G
- ⏳ Li (Z=3) - data available in JSON
- ⏳ Be (Z=4) - data available in JSON
- ⏳ B (Z=5) - data available in JSON
- ⏳ F (Z=9) - data available in JSON
- ⏳ Ne (Z=10) - data available in JSON
- ⏳ Na-Ar (Z=11-18) - data available in JSON

### Remaining Elements for 6-31G/6-31G*
- ⏳ He (Z=2) - data available in JSON
- ⏳ Li-B (Z=3-5) - data available in JSON
- ⏳ N, F (Z=7, 9) - data available in JSON
- ⏳ Ne (Z=10) - data available in JSON
- ⏳ Na-Ar (Z=11-18) - data available in JSON

## Data Source

Complete basis set data for all elements H-Ar is available in:
- `basis_sets_h_to_ar.json` (90 KB) - Machine-readable format
- `COMPLETE_BASIS_SETS_H_AR.txt` (52 KB) - Human-readable format
- `BASIS_SETS_SUMMARY.txt` (11 KB) - Quick reference

All data sourced from Basis Set Exchange: https://www.basissetexchange.org/

## Implementation Guide

### Adding a New Element

To add basis set data for an element (e.g., Nitrogen 6-31G):

1. **Find data** in `basis_sets_h_to_ar.json`
2. **Add to basis_data.rs**:
```rust
pub const N_631G: &[ShellData] = &[
    ShellData {
        shell_type: ShellType::S,
        exponents: &[...], // Core S shell (6 primitives)
        coefficients: &[&[...]], // S coefficients
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[...], // Valence SP shell (3 primitives)
        coefficients: &[
            &[...], // S coefficients
            &[...], // P coefficients
        ],
    },
    ShellData {
        shell_type: ShellType::SP,
        exponents: &[...], // Diffuse SP shell (1 primitive)
        coefficients: &[
            &[1.0], // S coefficient
            &[1.0], // P coefficient
        ],
    },
];
```

3. **Add to lookup function** in `basis.rs`:
```rust
fn get_631g_shells(atomic_number: u32) -> &'static [ShellData] {
    match atomic_number {
        1 => H_631G,
        6 => C_631G,
        7 => N_631G, // ADD THIS LINE
        8 => O_631G,
        _ => panic!("6-31G not yet implemented for Z={}", atomic_number),
    }
}
```

4. **Test**:
```rust
#[test]
fn test_nh3_631g() {
    let nh3 = Molecule::nh3();
    let basis = build_basis(&nh3, BasisSet::_631G);
    // N: 9 functions, 3 H: 2 each = 15 total
    assert_eq!(basis.len(), 15);
}
```

## Usage Examples

### Basic Usage
```rust
use rmolints::basis::{build_basis, BasisSet};
use rmolints::molecule::Molecule;

// H2O with STO-3G (7 basis functions)
let h2o = Molecule::h2o();
let basis_sto3g = build_basis(&h2o, BasisSet::STO3G);

// H2O with 6-31G (13 basis functions)
let basis_631g = build_basis(&h2o, BasisSet::_631G);

// H2O with 6-31G* (19 basis functions - adds 6 d functions on O)
let basis_631g_star = build_basis(&h2o, BasisSet::_631GStar);
```

### Computing ERIs with Different Basis Sets
```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};

// Compute ERIs with 6-31G* basis
let basis = build_basis(&molecule, BasisSet::_631GStar);
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPople);
```

## Basis Function Counts

| Molecule | STO-3G | 6-31G | 6-31G* |
|----------|--------|-------|--------|
| H2 | 2 | 4 | 4 |
| H2O | 7 | 13 | 19 |
| NH3 | 9 | 15 | 21 |
| CH4 | 9 | 17 | 23 |
| C6H6 (benzene) | 36 | 66 | 102 |

## Performance Considerations

### 6-31G vs STO-3G
- **Basis functions**: ~2-3x more basis functions
- **ERIs**: Scales as O(N⁴), so ~16-81x more integrals
- **Accuracy**: Significantly better for molecular properties
- **Cost**: Higher computational cost but worth it for production calculations

### 6-31G* vs 6-31G
- **Heavy atoms**: +6 d functions per heavy atom (not H or He)
- **Polarization**: Essential for accurate geometries, energies, properties
- **Standard choice**: 6-31G* is the minimum recommended for most calculations

## References

- **STO-3G**: W.J. Hehre, R.F. Stewart and J.A. Pople, *J. Chem. Phys.* **56**, 2657 (1969)
- **6-31G**: W.J. Hehre, R. Ditchfield and J.A. Pople, *J. Chem. Phys.* **56**, 2257 (1972)
- **6-31G***: P.C. Hariharan and J.A. Pople, *Theoret. Chimica Acta* **28**, 213 (1973)
- **Basis Set Exchange**: https://www.basissetexchange.org/
