use rmolints::common::{Primitive, Vec3, CGBF};
use rmolints::hgp_opt::electron_repulsion_hgp_opt;
#[cfg(feature = "simd")]
use rmolints::hgp_simd::electron_repulsion_hgp_simd;
use std::time::Instant;

fn create_orbital(origin: Vec3, shell: (i32, i32, i32)) -> CGBF {
    CGBF {
        origin,
        shell,
        primitives: vec![
            Primitive {
                exponent: 0.5,
                coefficient: 0.444635,
            },
            Primitive {
                exponent: 0.3,
                coefficient: 0.535328,
            },
            Primitive {
                exponent: 0.15,
                coefficient: 0.154329,
            },
        ],
    }
}

fn benchmark_orbital_combination(
    label: &str,
    a: &CGBF,
    b: &CGBF,
    c: &CGBF,
    d: &CGBF,
    iterations: usize,
) {
    // Warmup
    for _ in 0..100 {
        let _ = electron_repulsion_hgp_opt(a, b, c, d);
    }

    // Benchmark HGP-Opt (scalar)
    let start = Instant::now();
    for _ in 0..iterations {
        let _ = electron_repulsion_hgp_opt(a, b, c, d);
    }
    let scalar_time = start.elapsed().as_secs_f64() * 1_000_000.0; // microseconds
    let scalar_per_call = scalar_time / iterations as f64;

    #[cfg(feature = "simd")]
    {
        // Warmup SIMD
        for _ in 0..100 {
            let _ = electron_repulsion_hgp_simd(a, b, c, d);
        }

        // Benchmark HGP-SIMD
        let start = Instant::now();
        for _ in 0..iterations {
            let _ = electron_repulsion_hgp_simd(a, b, c, d);
        }
        let simd_time = start.elapsed().as_secs_f64() * 1_000_000.0;
        let simd_per_call = simd_time / iterations as f64;

        let speedup = scalar_per_call / simd_per_call;
        let speedup_str = if speedup > 1.0 {
            format!("\x1b[32m{:.2}x faster\x1b[0m", speedup)
        } else {
            format!("\x1b[31m{:.2}x slower\x1b[0m", 1.0 / speedup)
        };

        println!(
            "{:20} | Scalar: {:8.3} ns | SIMD: {:8.3} ns | {}",
            label,
            scalar_per_call * 1000.0,
            simd_per_call * 1000.0,
            speedup_str
        );
    }

    #[cfg(not(feature = "simd"))]
    {
        println!(
            "{:20} | Scalar: {:8.3} ns | SIMD: N/A (feature disabled)",
            label,
            scalar_per_call * 1000.0
        );
    }
}

fn main() {
    println!("\n=== SIMD Performance Analysis ===\n");
    println!("Comparing HGP-Opt (scalar) vs HGP-SIMD for different orbital combinations\n");

    let origin_a = Vec3::new(0.0, 0.0, 0.0);
    let origin_b = Vec3::new(1.4, 0.0, 0.0);
    let origin_c = Vec3::new(0.7, 1.2, 0.0);
    let origin_d = Vec3::new(0.7, 0.0, 1.0);

    println!("{:20} | {:^18} | {:^18} | Speedup", "Orbital Type", "Scalar Time", "SIMD Time");
    println!("{}", "-".repeat(80));

    // S orbitals (angular momentum = 0)
    let s1 = create_orbital(origin_a, (0, 0, 0));
    let s2 = create_orbital(origin_b, (0, 0, 0));
    benchmark_orbital_combination("(ss|ss)", &s1, &s2, &s1, &s2, 10000);

    // P orbitals (angular momentum = 1)
    let px1 = create_orbital(origin_a, (1, 0, 0));
    let py1 = create_orbital(origin_a, (0, 1, 0));
    let px2 = create_orbital(origin_b, (1, 0, 0));
    benchmark_orbital_combination("(pp|pp)", &px1, &py1, &px2, &px1, 10000);

    // Mixed S and P
    benchmark_orbital_combination("(sp|sp)", &s1, &px1, &s2, &px2, 10000);

    // D orbitals (angular momentum = 2)
    let d_xx = create_orbital(origin_a, (2, 0, 0));
    let d_yy = create_orbital(origin_b, (0, 2, 0));
    benchmark_orbital_combination("(dd|dd)", &d_xx, &d_yy, &d_xx, &d_yy, 5000);

    // Mixed P and D
    benchmark_orbital_combination("(pd|pd)", &px1, &d_xx, &py1, &d_yy, 5000);

    // High angular momentum (sp, pp, dd combinations)
    let s_far = create_orbital(origin_c, (0, 0, 0));
    let p_far = create_orbital(origin_d, (1, 1, 0));
    benchmark_orbital_combination("(sp|sp) distant", &s1, &px1, &s_far, &p_far, 10000);

    println!("\n{}", "-".repeat(80));

    #[cfg(feature = "simd")]
    {
        println!("\n\x1b[1mAnalysis:\x1b[0m");
        println!("- Green values (>1.0x) indicate SIMD is faster");
        println!("- Red values (<1.0x) indicate SIMD is slower");
        println!("- SIMD works best when VRR inner loops have many iterations (high angular momentum)");
        println!("- S orbitals have minimal VRR work, so SIMD overhead dominates");
        println!("- P and D orbitals have more VRR work, better SIMD efficiency");
    }

    #[cfg(not(feature = "simd"))]
    {
        println!("\nSIMD feature is disabled. Run with: cargo run --release --features simd --example simd_benchmark");
    }

    println!("\n=== SIMD Efficiency Test ===\n");

    #[cfg(feature = "simd")]
    {
        // Test with real molecule basis functions
        use rmolints::basis::build_sto3g_basis;
        use rmolints::molecule::Molecule;

        let h2o = Molecule::h2o();
        let basis = build_sto3g_basis(&h2o);

        println!("Testing on H2O molecule with STO-3G basis ({} functions)\n", basis.len());

        // Sample a few representative integrals
        let test_cases = vec![
            (0, 0, 0, 0, "O(s) - O(s)"),
            (1, 1, 1, 1, "O(px) - O(px)"),
            (2, 3, 2, 3, "O(py) - O(pz)"),
            (0, 1, 2, 3, "O(s/px) - O(py/pz)"),
        ];

        for (i, j, k, l, label) in test_cases {
            let a = &basis[i];
            let b = &basis[j];
            let c = &basis[k];
            let d = &basis[l];

            benchmark_orbital_combination(label, a, b, c, d, 10000);
        }
    }

    println!();
}
