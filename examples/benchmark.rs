use rmolints::common::*;
use rmolints::{two_electron, rys, hgp};
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

fn benchmark_method<F>(name: &str, iterations: usize, f: F) -> f64
where
    F: Fn() -> f64,
{
    // Warmup
    for _ in 0..100 {
        let _ = f();
    }

    let start = Instant::now();
    let mut sum = 0.0;
    for _ in 0..iterations {
        sum += f();
    }
    let elapsed = start.elapsed();

    let time_per_call = elapsed.as_nanos() as f64 / iterations as f64;
    println!("{:20} {:8.2} ns/call  (result: {:.6})", name, time_per_call, sum / iterations as f64);

    time_per_call
}

fn main() {
    println!("\n=== Two-Electron Integral Method Benchmark ===\n");

    let iterations = 10000;
    let origin = Vec3::new(0.0, 0.0, 0.0);
    let origin2 = Vec3::new(0.0, 0.0, 1.0);

    // Test 1: Four s-orbitals at same position
    println!("Test 1: Four s-orbitals (same position)");
    println!("{}", "-".repeat(70));
    let s1 = s_orbital(1.0, origin);

    let t_std = benchmark_method("Standard (THO)", iterations, || {
        two_electron::electron_repulsion(&s1, &s1, &s1, &s1)
    });

    let t_rys = benchmark_method("Rys Quadrature", iterations, || {
        rys::electron_repulsion_rys(&s1, &s1, &s1, &s1)
    });

    let t_hgp = benchmark_method("Head-Gordon-Pople", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &s1, &s1)
    });

    let t_hgp = benchmark_method("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &s1, &s1)
    });

    println!("Speedup vs Standard: Rys={:.2}x, HGP={:.2}x, HGP-Opt={:.2}x\n", t_std/t_rys, t_std/t_hgp, t_std/t_hgp);

    // Test 2: Four s-orbitals (separated)
    println!("Test 2: Four s-orbitals (separated)");
    println!("{}", "-".repeat(70));
    let s2 = s_orbital(1.0, origin2);

    let t_std = benchmark_method("Standard (THO)", iterations, || {
        two_electron::electron_repulsion(&s1, &s2, &s1, &s2)
    });

    let t_rys = benchmark_method("Rys Quadrature", iterations, || {
        rys::electron_repulsion_rys(&s1, &s2, &s1, &s2)
    });

    let t_hgp = benchmark_method("Head-Gordon-Pople", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s2, &s1, &s2)
    });

    let t_hgp = benchmark_method("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s2, &s1, &s2)
    });

    println!("Speedup vs Standard: Rys={:.2}x, HGP={:.2}x, HGP-Opt={:.2}x\n", t_std/t_rys, t_std/t_hgp, t_std/t_hgp);

    // Test 3: Mixed s and p orbitals
    println!("Test 3: Mixed s and p orbitals");
    println!("{}", "-".repeat(70));
    let px = p_orbital(1.0, origin, 'x');

    let t_std = benchmark_method("Standard (THO)", iterations, || {
        two_electron::electron_repulsion(&s1, &s1, &px, &px)
    });

    let t_rys = benchmark_method("Rys Quadrature", iterations, || {
        rys::electron_repulsion_rys(&s1, &s1, &px, &px)
    });

    let t_hgp = benchmark_method("Head-Gordon-Pople", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &px, &px)
    });

    let t_hgp = benchmark_method("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&s1, &s1, &px, &px)
    });

    println!("Speedup vs Standard: Rys={:.2}x, HGP={:.2}x, HGP-Opt={:.2}x\n", t_std/t_rys, t_std/t_hgp, t_std/t_hgp);

    // Test 4: Four p orbitals
    println!("Test 4: Four p orbitals");
    println!("{}", "-".repeat(70));
    let py = p_orbital(1.0, origin, 'y');

    let t_std = benchmark_method("Standard (THO)", iterations, || {
        two_electron::electron_repulsion(&px, &py, &px, &py)
    });

    let t_rys = benchmark_method("Rys Quadrature", iterations, || {
        rys::electron_repulsion_rys(&px, &py, &px, &py)
    });

    let t_hgp = benchmark_method("Head-Gordon-Pople", iterations, || {
        hgp::electron_repulsion_hgp(&px, &py, &px, &py)
    });

    let t_hgp = benchmark_method("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&px, &py, &px, &py)
    });

    println!("Speedup vs Standard: Rys={:.2}x, HGP={:.2}x, HGP-Opt={:.2}x\n", t_std/t_rys, t_std/t_hgp, t_std/t_hgp);

    // Test 5: d orbitals (higher angular momentum)
    println!("Test 5: d orbitals (l=2, higher angular momentum)");
    println!("{}", "-".repeat(70));
    let dxx = d_orbital(1.0, origin, (2, 0, 0));
    let dyy = d_orbital(1.0, origin, (0, 2, 0));

    let t_std = benchmark_method("Standard (THO)", iterations, || {
        two_electron::electron_repulsion(&dxx, &s1, &dyy, &s1)
    });

    let t_rys = benchmark_method("Rys Quadrature", iterations, || {
        rys::electron_repulsion_rys(&dxx, &s1, &dyy, &s1)
    });

    let t_hgp = benchmark_method("Head-Gordon-Pople", iterations, || {
        hgp::electron_repulsion_hgp(&dxx, &s1, &dyy, &s1)
    });

    let t_hgp = benchmark_method("HGP Optimized", iterations, || {
        hgp::electron_repulsion_hgp(&dxx, &s1, &dyy, &s1)
    });

    println!("Speedup vs Standard: Rys={:.2}x, HGP={:.2}x, HGP-Opt={:.2}x\n", t_std/t_rys, t_std/t_hgp, t_std/t_hgp);

    // Test 6: Four d orbitals (maximum angular momentum for this test)
    println!("Test 6: Four d orbitals (maximum angular momentum)");
    println!("{}", "-".repeat(70));
    let dzz = d_orbital(1.0, origin, (0, 0, 2));

    let t_std = benchmark_method("Standard (THO)", iterations/2, || {
        two_electron::electron_repulsion(&dxx, &dyy, &dzz, &s1)
    });

    let t_rys = benchmark_method("Rys Quadrature", iterations/2, || {
        rys::electron_repulsion_rys(&dxx, &dyy, &dzz, &s1)
    });

    let t_hgp = benchmark_method("Head-Gordon-Pople", iterations/2, || {
        hgp::electron_repulsion_hgp(&dxx, &dyy, &dzz, &s1)
    });

    let t_hgp = benchmark_method("HGP Optimized", iterations/2, || {
        hgp::electron_repulsion_hgp(&dxx, &dyy, &dzz, &s1)
    });

    println!("Speedup vs Standard: Rys={:.2}x, HGP={:.2}x, HGP-Opt={:.2}x\n", t_std/t_rys, t_std/t_hgp, t_std/t_hgp);

    println!("=== Benchmark Complete ===");
}
