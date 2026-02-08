use rmolints::common::*;
use rmolints::parallel::{compute_eri_tensor_parallel, compute_eris_parallel, ERIMethod};
use rmolints::two_electron;
use std::time::Instant;

fn make_s_orbital(alpha: f64, pos: Vec3) -> CGBF {
    CGBF {
        origin: pos,
        shell: (0, 0, 0),
        primitives: vec![Primitive {
            exponent: alpha,
            coefficient: 1.0,
        }],
    }
}

fn make_p_orbital(alpha: f64, pos: Vec3, direction: char) -> CGBF {
    let shell = match direction {
        'x' => (1, 0, 0),
        'y' => (0, 1, 0),
        'z' => (0, 0, 1),
        _ => panic!("Invalid direction"),
    };
    CGBF {
        origin: pos,
        shell,
        primitives: vec![Primitive {
            exponent: alpha,
            coefficient: 1.0,
        }],
    }
}

fn benchmark_serial_vs_parallel(basis: &[CGBF], indices: &[(usize, usize, usize, usize)]) {
    let n_integrals = indices.len();

    // Serial benchmark
    let start = Instant::now();
    let mut serial_results = Vec::new();
    for &(i, j, k, l) in indices {
        let val = two_electron::electron_repulsion(&basis[i], &basis[j], &basis[k], &basis[l]);
        serial_results.push(val);
    }
    let serial_time = start.elapsed();

    // Parallel benchmark
    let start = Instant::now();
    let parallel_results = compute_eris_parallel(basis, indices, ERIMethod::Standard);
    let parallel_time = start.elapsed();

    // Verify results match
    let max_diff = serial_results
        .iter()
        .zip(parallel_results.iter())
        .map(|(s, p)| (s - p).abs())
        .fold(0.0f64, |a, b| a.max(b));

    println!("  Integrals computed: {}", n_integrals);
    println!("  Serial time:        {:7.2} ms", serial_time.as_secs_f64() * 1000.0);
    println!("  Parallel time:      {:7.2} ms", parallel_time.as_secs_f64() * 1000.0);
    println!("  Speedup:            {:.2}x", serial_time.as_secs_f64() / parallel_time.as_secs_f64());
    println!("  Max difference:     {:.2e}", max_diff);
}

fn main() {
    println!("\n=== Parallel Two-Electron Integral Benchmark ===\n");

    // Test 1: Small basis set (2 s-orbitals)
    println!("Test 1: Small basis (2 s-orbitals, 6 unique integrals)");
    println!("{}", "-".repeat(70));
    let s1 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 0.0));
    let s2 = make_s_orbital(1.0, Vec3::new(0.0, 0.0, 1.0));
    let basis = vec![s1, s2];

    // Generate all unique integrals
    let mut indices = Vec::new();
    for i in 0..2 {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..=i {
                let l_max = if k == i { j } else { k };
                for l in 0..=l_max {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        indices.push((i, j, k, l));
                    }
                }
            }
        }
    }

    benchmark_serial_vs_parallel(&basis, &indices);
    println!();

    // Test 2: Medium basis set (4 s-orbitals)
    println!("Test 2: Medium basis (4 s-orbitals, 35 unique integrals)");
    println!("{}", "-".repeat(70));
    let basis: Vec<CGBF> = (0..4)
        .map(|i| make_s_orbital(1.0, Vec3::new(0.0, 0.0, i as f64 * 0.5)))
        .collect();

    let mut indices = Vec::new();
    for i in 0..4 {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..=i {
                let l_max = if k == i { j } else { k };
                for l in 0..=l_max {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        indices.push((i, j, k, l));
                    }
                }
            }
        }
    }

    benchmark_serial_vs_parallel(&basis, &indices);
    println!();

    // Test 3: Larger basis set (6 s-orbitals)
    println!("Test 3: Larger basis (6 s-orbitals, 126 unique integrals)");
    println!("{}", "-".repeat(70));
    let basis: Vec<CGBF> = (0..6)
        .map(|i| make_s_orbital(1.0, Vec3::new(0.0, 0.0, i as f64 * 0.4)))
        .collect();

    let mut indices = Vec::new();
    for i in 0..6 {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..=i {
                let l_max = if k == i { j } else { k };
                for l in 0..=l_max {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        indices.push((i, j, k, l));
                    }
                }
            }
        }
    }

    benchmark_serial_vs_parallel(&basis, &indices);
    println!();

    // Test 4: Even larger basis (8 s-orbitals)
    println!("Test 4: Large basis (8 s-orbitals, 330 unique integrals)");
    println!("{}", "-".repeat(70));
    let basis: Vec<CGBF> = (0..8)
        .map(|i| make_s_orbital(1.0, Vec3::new(0.0, 0.0, i as f64 * 0.3)))
        .collect();

    let mut indices = Vec::new();
    for i in 0..8 {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..=i {
                let l_max = if k == i { j } else { k };
                for l in 0..=l_max {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        indices.push((i, j, k, l));
                    }
                }
            }
        }
    }

    benchmark_serial_vs_parallel(&basis, &indices);
    println!();

    // Test 5: Mixed s and p orbitals (more expensive per integral)
    println!("Test 5: Mixed basis (4 s + 4 p orbitals, 495 unique integrals)");
    println!("{}", "-".repeat(70));
    let mut basis = Vec::new();
    for i in 0..4 {
        basis.push(make_s_orbital(1.0, Vec3::new(0.0, 0.0, i as f64 * 0.5)));
    }
    for i in 0..4 {
        basis.push(make_p_orbital(1.0, Vec3::new(0.0, 0.0, i as f64 * 0.5), 'z'));
    }

    let n = basis.len();
    let mut indices = Vec::new();
    for i in 0..n {
        for j in 0..=i {
            let ij = i * (i + 1) / 2 + j;
            for k in 0..=i {
                let l_max = if k == i { j } else { k };
                for l in 0..=l_max {
                    let kl = k * (k + 1) / 2 + l;
                    if ij >= kl {
                        indices.push((i, j, k, l));
                    }
                }
            }
        }
    }

    benchmark_serial_vs_parallel(&basis, &indices);
    println!();

    // Test 6: Full tensor computation
    println!("Test 6: Full tensor computation (10 s-orbitals)");
    println!("{}", "-".repeat(70));
    let basis: Vec<CGBF> = (0..10)
        .map(|i| make_s_orbital(1.0, Vec3::new(0.0, 0.0, i as f64 * 0.3)))
        .collect();

    let start = Instant::now();
    let tensor = compute_eri_tensor_parallel(&basis, ERIMethod::Standard);
    let parallel_time = start.elapsed();

    println!("  Basis functions:    {}", basis.len());
    println!("  Unique integrals:   {}", tensor.len());
    println!("  Parallel time:      {:7.2} ms", parallel_time.as_secs_f64() * 1000.0);
    println!("  Time per integral:  {:7.2} µs", parallel_time.as_secs_f64() * 1e6 / tensor.len() as f64);
    println!();

    println!("=== Benchmark Complete ===");
    println!("\nNOTE: Speedup depends on:");
    println!("  - Number of CPU cores available");
    println!("  - Number of integrals (more integrals = better parallelization)");
    println!("  - Integral complexity (higher angular momentum = more work per integral)");
}
