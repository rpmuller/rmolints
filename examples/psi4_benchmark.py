#!/usr/bin/env python3
"""
Benchmark Psi4 integral computation for comparison with Rust rmolints.
"""
import psi4
import time
import numpy as np

# Suppress Psi4 output
psi4.core.be_quiet()

def benchmark_molecule(name, geometry, basis_name):
    print(f"\n{'='*70}")
    print(f"Molecule: {name} with {basis_name}")
    print('='*70)

    # Create molecule
    mol = psi4.geometry(geometry)

    # Build basis set
    basis = psi4.core.BasisSet.build(mol, "BASIS", basis_name)
    n = basis.nbf()

    print(f"Basis functions: {n}")

    # Calculate number of unique integrals (8-fold symmetry)
    n_unique = n * (n + 1) * (n * (n + 1) // 2 + 1) // 8
    print(f"Unique integrals: {n_unique}")

    # Create integral helper
    mints = psi4.core.MintsHelper(basis)

    print(f"\n{'-'*70}")
    print("Psi4 ERI computation:")

    # Warmup
    _ = mints.ao_eri()

    # Benchmark
    start = time.time()
    eris = mints.ao_eri()
    elapsed = time.time() - start

    # Get the actual numpy array
    eri_array = np.asarray(eris)
    total_sum = eri_array.sum()

    us_per_integral = (elapsed * 1e6) / n_unique

    print(f"  Time:        {elapsed:.3f} s")
    print(f"  Integrals:   {n_unique}")
    print(f"  µs/integral: {us_per_integral:.2f}")
    print(f"  Total sum:   {total_sum:.6f}")

    return elapsed, us_per_integral, n_unique

def main():
    print("Psi4 Integral Benchmark")
    print("="*70)
    print(f"Psi4 version: {psi4.__version__}")
    print()

    results = []

    # H2 molecule
    h2_geom = """
    H  0.000000000000  0.000000000000  0.000000000000
    H  0.000000000000  0.000000000000  1.400000000000
    units bohr
    """
    elapsed, us_per_int, count = benchmark_molecule("H2", h2_geom, "sto-3g")
    results.append(("H2", "STO-3G", count, elapsed, us_per_int))

    # H2O molecule
    h2o_geom = """
    O  0.000000000000  0.000000000000  0.000000000000
    H  0.758602000000  0.504284000000  0.000000000000
    H -0.758602000000  0.504284000000  0.000000000000
    units angstrom
    """
    elapsed, us_per_int, count = benchmark_molecule("H2O", h2o_geom, "sto-3g")
    results.append(("H2O", "STO-3G", count, elapsed, us_per_int))

    # H2O with 6-31G*
    elapsed, us_per_int, count = benchmark_molecule("H2O", h2o_geom, "6-31g*")
    results.append(("H2O", "6-31G*", count, elapsed, us_per_int))

    # H2O with 6-31G**
    elapsed, us_per_int, count = benchmark_molecule("H2O", h2o_geom, "6-31g**")
    results.append(("H2O", "6-31G**", count, elapsed, us_per_int))

    # NH3 molecule
    nh3_geom = """
    N  0.000000000000  0.000000000000  0.000000000000
    H  0.000000000000  0.544277000000  0.893409000000
    H  0.865400000000 -0.282278000000  0.893409000000
    H -0.865400000000 -0.282278000000  0.893409000000
    units angstrom
    """
    elapsed, us_per_int, count = benchmark_molecule("NH3", nh3_geom, "sto-3g")
    results.append(("NH3", "STO-3G", count, elapsed, us_per_int))

    # Summary
    print(f"\n{'='*70}")
    print("SUMMARY - Psi4 Integral Computation")
    print('='*70)
    print(f"{'Molecule':<10} {'Basis':<10} {'Integrals':<12} {'Time (s)':<10} {'µs/int':<10}")
    print('-'*70)
    for mol, basis, count, elapsed, us_per_int in results:
        print(f"{mol:<10} {basis:<10} {count:<12} {elapsed:<10.3f} {us_per_int:<10.2f}")

    print(f"\nNote: Psi4 uses libint2 integral engine (highly optimized C++)")

if __name__ == "__main__":
    main()
