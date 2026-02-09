# Basis Set Scaling Analysis

## Basis Function Counts

### STO-3G vs 6-31G(d,p) Comparison

| Molecule | STO-3G | 6-31G(d,p) | Ratio | Notes |
|----------|--------|------------|-------|-------|
| H2       | 2      | 10         | 5.0x  | Each H: 1s → 2s + 3p |
| H2O      | 7      | 25         | 3.6x  | O: 5 → 15, H: 1 → 5 each |
| NH3      | 9      | 32         | 3.6x  | N: 5 → 15, H: 1 → 5 each |
| Benzene  | 36     | 120        | 3.3x  | C: 5 → 15, H: 1 → 5 each |

### Basis Functions per Atom

| Element | STO-3G | 6-31G | 6-31G(d) | 6-31G(d,p) |
|---------|--------|-------|----------|------------|
| H       | 1      | 2     | 2        | **5** (2s + 3p) |
| C, N, O | 5      | 9     | **15** (9 + 6d) | 15 |

## ERI Scaling

ERIs scale as O(N⁴) with basis set size:

| Molecule | STO-3G ERIs | 6-31G(d,p) ERIs | Ratio |
|----------|-------------|-----------------|-------|
| H2       | ~6          | ~3,000          | ~500x |
| H2O      | ~400        | ~100,000        | ~250x |
| Benzene  | ~220,000    | ~26,000,000     | ~120x |

**Key insight**: Moving from STO-3G to 6-31G(d,p) increases computational cost by 100-500x
for the same molecule!

## Performance Expectations

Based on STO-3G benzene performance (815 ms for 222,111 ERIs):

### Estimated 6-31G(d,p) Benzene Performance
- **Basis functions**: 120 (vs 36 for STO-3G)
- **Unique ERIs**: ~26,000,000 (vs 222,111)
- **Scaling factor**: ~120x more ERIs
- **Estimated time**: 815 ms × 120 = **~98 seconds** (1.6 minutes)

### For H2O (more manageable)
- **Basis functions**: 25 (vs 7 for STO-3G)
- **Unique ERIs**: ~100,000 (vs 406)
- **Scaling factor**: ~250x more ERIs  
- **Estimated time**: 2 ms × 250 = **~500 ms** (half second)

## Recommendations

1. **For testing/development**: Use STO-3G or small molecules with 6-31G(d,p)
2. **For production**: 6-31G(d,p) is the minimum recommended basis set
3. **For larger molecules**: Consider integral screening to skip near-zero ERIs
4. **Benchmarking**: Focus on H2, H2O, NH3 with 6-31G(d,p)

## Why 6-31G(d,p) Matters

- **Standard in quantum chemistry**: Minimal acceptable basis for published work
- **Polarization essential**: d functions on heavy atoms, p on hydrogens
- **Captures bonding accurately**: Allows electron density to polarize properly
- **Industry standard**: Used in drug design, materials science, etc.

The computational cost is necessary for chemical accuracy!
