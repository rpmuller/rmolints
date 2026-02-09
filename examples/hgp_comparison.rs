use rmolints::common::*;
use rmolints::{hgp, hgp};
use std::time::Instant;

fn s_orbital(alpha: f64, origin: Vec3) -> CGBF {
    CGBF {
        origin,
        shell: (0, 0, 0),
        primitives: vec![Primitive {
            exponent: alpha,
            coefficient: 1.0,
        }],
    }
}

fn p_orbital(alpha: f64, origin: Vec3, direction: char) -> CGBF {
    let shell = match direction {
        'x' => (1, 0, 0),
        'y' => (0, 1, 0),
        'z' => (0, 0, 1),
        _ => panic!("Invalid direction"),
    };
    CGBF {
        origin,
        shell,
        primitives: vec![Primitive {
            exponent: alpha,
            coefficient: 1.0,
        }],
    }
}

fn d_orbital(alpha: f64, origin: Vec3, lmn: (i32, i32, i32)) -> CGBF {
    CGBF {
        origin,
        shell: lmn,
        primitives: vec![Primitive {
            exponent: alpha,
            coefficient: 1.0,
        }],
    }
}

fn benchmark_hgp<F>(name: &str, iterations: usize, f: F) -> f64
where
    F: Fn() -> f64,
{
    // Warmup
    for _ in 0..10 {
        let _ = f();
    }

    let start = Instant::now();
    let mut sum = 0.0;
    for _ in 0..iterations {
        sum += f();
    }
    let elapsed = start.elapsed();

    let time_per_call = elapsed.as_nanos() as f64 / iterations as f64;
    println!(
        "{:25} {:10.2} ns/call  (result: {:.8})",
        name,
        time_per_call,
        sum / iterations as f64
    );

    time_per_call
}

fn main() {
    println!("\n=== HGP Old vs Optimized Comparison ===\n");

    let iterations = 1000;
    let origin = Vec3::new(0.0, 0.0, 0.0);
    let origin2 = Vec3::new(0.0, 0.0, 1.0);

    // Test 1: Four s-orbitals (same position)
    println!("Test 1: Four s-orbitals (same position)");
    println!("{}", "-".repeat(70));
    let s1 = s_orbital(1.0, origin);

    let t_old = benchmark_hgp("HGP Original", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &s1, &s1)
    });

    let t_opt = benchmark_hgp("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &s1, &s1)
    });

    println!("Speedup: {:.2}x\n", t_old / t_opt);

    // Test 2: Four s-orbitals (separated)
    println!("Test 2: Four s-orbitals (separated)");
    println!("{}", "-".repeat(70));
    let s2 = s_orbital(1.0, origin2);

    let t_old = benchmark_hgp("HGP Original", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s2, &s1, &s2)
    });

    let t_opt = benchmark_hgp("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s2, &s1, &s2)
    });

    println!("Speedup: {:.2}x\n", t_old / t_opt);

    // Test 3: Mixed s and p orbitals
    println!("Test 3: Mixed s and p orbitals");
    println!("{}", "-".repeat(70));
    let px = p_orbital(1.0, origin, 'x');

    let t_old = benchmark_hgp("HGP Original", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &px, &px)
    });

    let t_opt = benchmark_hgp("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &px, &px)
    });

    println!("Speedup: {:.2}x\n", t_old / t_opt);

    // Test 4: Four p orbitals
    println!("Test 4: Four p orbitals");
    println!("{}", "-".repeat(70));
    let py = p_orbital(1.0, origin, 'y');

    let t_old = benchmark_hgp("HGP Original", iterations, || {
        hgp::electron_repulsion_hgp(&px, &py, &px, &py)
    });

    let t_opt = benchmark_hgp("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&px, &py, &px, &py)
    });

    println!("Speedup: {:.2}x\n", t_old / t_opt);

    // Test 5: d orbitals (higher angular momentum)
    println!("Test 5: d orbitals (l=2, higher angular momentum)");
    println!("{}", "-".repeat(70));
    let dxx = d_orbital(1.0, origin, (2, 0, 0));
    let dyy = d_orbital(1.0, origin, (0, 2, 0));

    let t_old = benchmark_hgp("HGP Original", iterations, || {
        hgp::electron_repulsion_hgp(&dxx, &s1, &dyy, &s1)
    });

    let t_opt = benchmark_hgp("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&dxx, &s1, &dyy, &s1)
    });

    println!("Speedup: {:.2}x\n", t_old / t_opt);

    // Test 6: Four d orbitals (maximum angular momentum)
    println!("Test 6: Mixed d orbitals (highest complexity)");
    println!("{}", "-".repeat(70));
    let dzz = d_orbital(1.0, origin, (0, 0, 2));

    let t_old = benchmark_hgp("HGP Original", iterations / 2, || {
        hgp::electron_repulsion_hgp(&dxx, &dyy, &dzz, &s1)
    });

    let t_opt = benchmark_hgp("HGP Optimized", iterations / 2, || {
        hgp::electron_repulsion_hgp(&dxx, &dyy, &dzz, &s1)
    });

    println!("Speedup: {:.2}x\n", t_old / t_opt);

    println!("=== Benchmark Complete ===");
}
