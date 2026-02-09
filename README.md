# rmolints - Molecular Integrals Library

A high-performance molecular integrals library in Rust, ported from PyQuante2.

## Features

- **Complete integral library**: One-electron and two-electron integrals
- **Multiple ERI algorithms**: Five methods with different performance characteristics
- **High accuracy**: Matches PyQuante2 to 1e-5 precision
- **Excellent performance**: HGP-Opt method is 2x faster than alternatives
- **Parallel computation**: Multithreaded evaluation using Rayon
- **Well-tested**: 48/48 tests passing

## Quick Start

```rust
use rmolints::common::*;
use rmolints::{one_electron, hgp_opt};

// Create s-orbital
let s = CGBF {
    origin: Vec3::new(0.0, 0.0, 0.0),
    shell: (0, 0, 0),
    primitives: vec![Primitive { exponent: 1.0, coefficient: 1.0 }],
};

// One-electron integrals
let overlap = one_electron::overlap(&s, &s);
let kinetic = one_electron::kinetic(&s, &s);

// Two-electron integral (recommended method)
let eri = hgp_opt::electron_repulsion_hgp_opt(&s, &s, &s, &s);
```

## Performance (Benzene, 36 basis functions)

| Method | Time (ms) | Status |
|--------|-----------|--------|
| **HGP-Opt** | **815** | 🏆 **Fastest - recommended** |
| Rys (optimized) | 1,442 | Second best |
| Standard THO | 1,761 | Baseline |

**Recommendation**: Use HGP-Opt for all production code.

## Installation

```bash
# Standard build (stable Rust)
cargo build --release

# Run tests
cargo test --release

# Run benchmarks
cargo run --release --example molecule_benchmark
```

### Optional: SIMD Support

SIMD provides minimal benefit on ARM and is opt-in:

```bash
# Requires nightly Rust
cargo +nightly build --release --features simd
```

## Capabilities

### One-Electron Integrals
- Overlap: ⟨bra|ket⟩
- Kinetic: ⟨bra|T|ket⟩
- Nuclear: ⟨bra|V|ket⟩

### Two-Electron Integrals

Five methods available:
1. **HGP-Opt** (flat array) - 🏆 Fastest, recommended
2. **Rys** (quadrature) - Good for high angular momentum
3. **Standard** (THO) - Reference implementation
4. **HGP-SIMD** (optional) - Experimental, nightly required
5. **HGP Original** - Deprecated (7x slower)

### Supported Systems
- **Molecules**: H2, H2O, benzene (built-in)
- **Basis sets**: STO-3G
- **Orbitals**: s, p, d (tested); arbitrary angular momentum (theoretical)

## Usage Examples

### Parallel ERI Computation

```rust
use rmolints::parallel::{compute_eri_tensor_parallel, ERIMethod};
use rmolints::molecule::Molecule;
use rmolints::basis::build_sto3g_basis;

let molecule = Molecule::h2o();
let basis = build_sto3g_basis(&molecule);

// Compute all unique ERIs (exploits 8-fold symmetry)
let eris = compute_eri_tensor_parallel(&basis, ERIMethod::HeadGordonPopleOpt);
```

### Benchmarking

```bash
# Comprehensive molecular benchmarks
cargo run --release --example molecule_benchmark

# Performance ratio analysis
cargo run --release --example ratio_analysis

# Profiling with flamegraphs
cargo run --release --example profile_hgp_opt
```

## Performance Insights

### Key Learnings

1. **Memory layout matters**: Flat array storage is 2.4x faster than nested structures
2. **Algorithm selection**: HGP-Opt beats Rys and Standard on real molecules
3. **SIMD limitations**: Minimal benefit for typical quantum chemistry integrals
4. **Profiling results**: VRR stages 6-7 account for 75% of computation time

### Optimization History

- **Nested Vec → Flat Array**: 2.9x speedup
- **HashMap → Flat Array**: 2.4x additional speedup
- **Rys optimizations**: 1.9x speedup on Rys method
- **Total improvement**: 7x faster than original HGP

## Testing

```bash
# All tests
cargo test --release

# Module-specific
cargo test --release hgp_opt
cargo test --release parallel
```

All 48 tests pass with 1e-5 precision against PyQuante2 reference values.

## Documentation

See [`docs/DOCS.md`](docs/DOCS.md) for:
- Complete installation instructions
- Detailed performance analysis
- Algorithm descriptions
- Extended usage examples
- Profiling and optimization guide
- Next steps and roadmap

## References

- **PyQuante2**: https://github.com/rpmuller/pyquante2
- **THO**: Taketa, Huzinaga, O-ohata, *J. Phys. Soc. Japan* **21**, 2313 (1966)
- **Rys**: Augspurger et al., *J. Comp. Chem.* **11**(8), 972-977 (1990)
- **HGP**: Head-Gordon & Pople, *J. Chem. Phys.* **89**(9), 5777 (1988)

## License

Modified BSD license (same as PyQuante2)
