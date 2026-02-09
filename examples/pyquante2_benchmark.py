#!/usr/bin/env python3
"""
Benchmark PyQuante2 C-based integrals for comparison with Rust implementation.
"""
import time
from pyquante2 import molecule, basisset
from pyquante2.ints.integrals import twoe_integrals

def benchmark_molecule(name, mol, basis_name):
    print(f"\n{'='*70}")
    print(f"Molecule: {name} with {basis_name}")
    print('='*70)

    # Build basis set
    basis = basisset(mol, basis_name)
    n = len(basis)

    print(f"Basis functions: {n}")

    # Calculate number of unique integrals (8-fold symmetry)
    n_unique = n * (n + 1) * (n * (n + 1) // 2 + 1) // 8
    print(f"Unique integrals: {n_unique}")

    # Benchmark integral computation
    print(f"\n{'-'*70}")
    print("PyQuante2 C-based integrals (cints):")

    start = time.time()

    # Compute all integrals using PyQuante2's twoe_integrals class
    # This uses the C-based routines by default
    eris = twoe_integrals(basis)

    elapsed = time.time() - start

    # Count actual integrals computed
    eri_count = len(eris._2e_ints)
    total = eris._2e_ints.sum()

    us_per_integral = (elapsed * 1e6) / n_unique

    print(f"  Time:        {elapsed:.3f} s")
    print(f"  Integrals:   {n_unique}")
    print(f"  µs/integral: {us_per_integral:.2f}")
    print(f"  Total sum:   {total:.6f}")

    return elapsed, us_per_integral, n_unique

def main():
    print("PyQuante2 C-Integrals Benchmark")
    print("="*70)
    print("Testing C-based cints routines from PyQuante2\n")

    results = []

    # H2 molecule (using atomic numbers)
    h2 = molecule([(1, 0.0, 0.0, 0.0),
                   (1, 0.0, 0.0, 1.4)],
                  units='Bohr')
    elapsed, us_per_int, count = benchmark_molecule("H2", h2, "sto-3g")
    results.append(("H2", "STO-3G", count, elapsed, us_per_int))

    # H2O molecule (using atomic numbers: O=8, H=1)
    h2o = molecule([(8, 0.0, 0.0, 0.0),
                    (1, 0.758602, 0.504284, 0.0),
                    (1, -0.758602, 0.504284, 0.0)],
                   units='Angstrom')
    elapsed, us_per_int, count = benchmark_molecule("H2O", h2o, "sto-3g")
    results.append(("H2O", "STO-3G", count, elapsed, us_per_int))

    # H2O with 6-31G**
    elapsed, us_per_int, count = benchmark_molecule("H2O", h2o, "6-31g**")
    results.append(("H2O", "6-31G**", count, elapsed, us_per_int))

    # NH3 molecule (using atomic numbers: N=7, H=1)
    nh3 = molecule([(7, 0.0, 0.0, 0.0),
                    (1, 0.0, 0.544277, 0.893409),
                    (1, 0.8654, -0.282278, 0.893409),
                    (1, -0.8654, -0.282278, 0.893409)],
                   units='Angstrom')
    elapsed, us_per_int, count = benchmark_molecule("NH3", nh3, "sto-3g")
    results.append(("NH3", "STO-3G", count, elapsed, us_per_int))

    # Summary table
    print(f"\n{'='*70}")
    print("SUMMARY - PyQuante2 C-based Integrals")
    print('='*70)
    print(f"{'Molecule':<10} {'Basis':<10} {'Integrals':<12} {'Time (s)':<10} {'µs/int':<10}")
    print('-'*70)
    for mol, basis, count, elapsed, us_per_int in results:
        print(f"{mol:<10} {basis:<10} {count:<12} {elapsed:<10.3f} {us_per_int:<10.2f}")

    print("\nNote: PyQuante2 uses C-based routines from cints module by default")

if __name__ == "__main__":
    main()
