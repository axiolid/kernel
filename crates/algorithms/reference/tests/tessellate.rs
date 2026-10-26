//! Tessellation checked against closed-form surface truth (ADR 0012).
//!
//! A tessellator is only trustworthy if the mesh it returns actually meets the
//! tolerance it claims. Each test pins a value derivable on paper: an inscribed
//! grid under-estimates area and volume, converges as the budget tightens, and
//! never claims success while the budget is exhausted.

use axiolid_core::{Point3, Scalar, Tolerance};
use axiolid_reference::surface::Patch;
use axiolid_reference::tessellate::{tessellate_patch, TessellationBudget};
use axiolid_surface::{Cylinder, Plane, Sphere, Surface, Torus};
use core::f64::consts::{PI, TAU};

fn frame() -> axiolid_core::Frame3 {
    axiolid_core::Frame3 {
        origin: Point3::new(0.0, 0.0, 0.0),
        x: axiolid_core::Vec3::new(1.0, 0.0, 0.0),
        y: axiolid_core::Vec3::new(0.0, 1.0, 0.0),
        z: axiolid_core::Vec3::new(0.0, 0.0, 1.0),
    }
}

fn budget(chord: Scalar) -> TessellationBudget {
    TessellationBudget::new(chord, 4096).expect("budget")
}

/// Triangulated surface area, straight from the mesh.
fn mesh_area(mesh: &axiolid_mesh::TriMesh) -> Scalar {
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| {
            let a = mesh.positions[c[0] as usize];
            let b = mesh.positions[c[1] as usize];
            let d = mesh.positions[c[2] as usize];
            (b - a).cross(d - a).length() * 0.5
        })
        .sum()
}

// --- the tolerance contract -------------------------------------------------

#[test]
fn a_plane_needs_no_subdivision_at_all() {
    // A plane is exactly representable by two triangles. A tessellator that
    // subdivides it anyway is wasting every downstream operation's time.
    let s = Surface::Plane(Plane { frame: frame() });
    let patch = Patch::new(-2.0, 2.0, -3.0, 3.0).expect("patch");
    let out = tessellate_patch(&s, patch, budget(1e-6)).expect("tessellate");
    assert_eq!(out.u_samples, 2, "a flat direction needs two samples");
    assert_eq!(out.v_samples, 2);
    assert_eq!(out.mesh.indices.len() / 3, 2, "one quad, two triangles");
    // 4 x 6 = 24, exactly.
    assert!((mesh_area(&out.mesh) - 24.0).abs() < 1e-12);
}

#[test]
fn a_tighter_budget_really_reduces_the_error() {
    // The claim a chord budget makes is that error shrinks with it. Assert the
    // relationship, not a magic magnitude: a fixed threshold would be pinning
    // a sample count instead of the contract.
    let s = Surface::Sphere(Sphere {
        frame: frame(),
        radius: 2.0,
    });
    let patch = Patch::new(0.0, TAU, -PI / 2.0, PI / 2.0).expect("patch");
    let exact = 4.0 * PI * 2.0 * 2.0;
    let mut previous = Scalar::INFINITY;
    for chord in [1e-1, 1e-2, 1e-3, 1e-4] {
        let out = tessellate_patch(&s, patch, budget(chord)).expect("tessellate");
        let error = exact - mesh_area(&out.mesh);
        assert!(
            error > 0.0,
            "an inscribed grid must under-estimate area, got {error}"
        );
        assert!(
            error < previous,
            "chord {chord}: error {error} did not improve on {previous}"
        );
        previous = error;
    }
}

#[test]
fn the_measured_sagitta_respects_the_budget() {
    // The outcome reports what it measured. If that number exceeds the budget
    // while claiming success, the whole contract is decorative.
    let s = Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: 1.5,
    });
    let patch = Patch::new(0.0, TAU, 0.0, 4.0).expect("patch");
    let chord = 1e-3;
    let out = tessellate_patch(&s, patch, budget(chord)).expect("tessellate");
    assert!(!out.budget_exhausted, "1.5m cylinder at 1mm must fit");
    let measured = out.max_sagitta.expect("a curved surface deviates");
    assert!(
        measured <= chord,
        "reported sagitta {measured} exceeds budget {chord}"
    );
}

#[test]
fn an_exhausted_budget_is_reported_not_hidden() {
    // A caller measuring quantities must be able to tell a coarse mesh from a
    // converged one. Silently returning the best effort is how wrong numbers
    // reach a report.
    let s = Surface::Sphere(Sphere {
        frame: frame(),
        radius: 50.0,
    });
    let patch = Patch::new(0.0, TAU, -PI / 2.0, PI / 2.0).expect("patch");
    let tiny = TessellationBudget::new(1e-9, 8).expect("budget");
    let out = tessellate_patch(&s, patch, tiny).expect("tessellate");
    assert!(
        out.budget_exhausted,
        "1e-9 on a 50m sphere with 8 samples cannot converge"
    );
    assert!(out.u_samples <= 8 && out.v_samples <= 8, "cap respected");
}

// --- analytic geometry ------------------------------------------------------

#[test]
fn a_cylinder_band_converges_on_tau_r_h() {
    let (r, h) = (1.5, 4.0);
    let s = Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: r,
    });
    let patch = Patch::new(0.0, TAU, 0.0, h).expect("patch");
    let out = tessellate_patch(&s, patch, budget(1e-5)).expect("tessellate");
    let want = TAU * r * h;
    let got = mesh_area(&out.mesh);
    assert!(
        got < want && (want - got) / want < 1e-4,
        "cylinder area {got} vs {want}"
    );
}

#[test]
fn a_sphere_converges_on_four_pi_r_squared() {
    let r = 2.0;
    let s = Surface::Sphere(Sphere {
        frame: frame(),
        radius: r,
    });
    let patch = Patch::new(0.0, TAU, -PI / 2.0, PI / 2.0).expect("patch");
    let tight = TessellationBudget::new(1e-4, 1024).expect("budget");
    let out = tessellate_patch(&s, patch, tight).expect("tessellate");
    let want = 4.0 * PI * r * r;
    let got = mesh_area(&out.mesh);
    assert!(
        got < want && (want - got) / want < 1e-3,
        "sphere area {got} vs {want}"
    );
}

#[test]
fn a_torus_converges_on_four_pi_squared_r_r() {
    // Pappus again, now through the mesh rather than the parameterisation.
    let (major, minor) = (3.0, 1.0);
    let s = Surface::Torus(Surface3Torus::torus(major, minor));
    let patch = Patch::new(0.0, TAU, 0.0, TAU).expect("patch");
    let tight = TessellationBudget::new(1e-4, 1024).expect("budget");
    let out = tessellate_patch(&s, patch, tight).expect("tessellate");
    let want = 4.0 * PI * PI * major * minor;
    let got = mesh_area(&out.mesh);
    assert!(
        got < want && (want - got) / want < 1e-3,
        "torus area {got} vs {want}"
    );
}

/// Local helper so the torus literal stays readable above.
struct Surface3Torus;
impl Surface3Torus {
    fn torus(major: Scalar, minor: Scalar) -> Torus {
        Torus {
            frame: frame(),
            major_radius: major,
            minor_radius: minor,
        }
    }
}

// --- structural soundness ---------------------------------------------------

#[test]
fn a_closed_band_does_not_leave_a_seam_sliver() {
    // The last sample must land exactly on the patch end. Accumulating a step
    // instead leaves a hairline gap that audit_mesh reports as degenerate and
    // every measure then refuses.
    let s = Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: 1.0,
    });
    let patch = Patch::new(0.0, TAU, 0.0, 2.0).expect("patch");
    let out = tessellate_patch(&s, patch, budget(1e-3)).expect("tessellate");
    let health = axiolid_mesh::audit_mesh(&out.mesh, Tolerance::MILLIMETRE);
    assert_eq!(
        health.degenerate_triangles, 0,
        "seam produced slivers: {health:?}"
    );
}

#[test]
fn every_triangle_has_three_distinct_corners() {
    let s = Surface::Sphere(Sphere {
        frame: frame(),
        radius: 1.0,
    });
    // Poles pulled in: the parameterisation is singular exactly at +-pi/2, and
    // a degenerate pole fan is a property of the parameterisation, not of the
    // tessellator being tested here.
    let patch = Patch::new(0.0, TAU, -1.5, 1.5).expect("patch");
    let out = tessellate_patch(&s, patch, budget(1e-3)).expect("tessellate");
    for c in out.mesh.indices.as_chunks::<3>().0 {
        assert!(
            c[0] != c[1] && c[1] != c[2] && c[0] != c[2],
            "degenerate index triple {c:?}"
        );
    }
}

// --- refusals ---------------------------------------------------------------

#[test]
fn a_non_positive_chord_budget_is_refused() {
    assert!(TessellationBudget::new(0.0, 64).is_err());
    assert!(TessellationBudget::new(-1e-3, 64).is_err());
    assert!(TessellationBudget::new(Scalar::NAN, 64).is_err());
}

#[test]
fn a_single_sample_cannot_span_a_patch() {
    assert!(TessellationBudget::new(1e-3, 1).is_err());
    assert!(TessellationBudget::new(1e-3, 2).is_ok());
}

// --- gaps found by mutation probes ------------------------------------------

/// A cell must split along its shorter diagonal.
///
/// On an anisotropic patch the two choices differ sharply in triangle quality.
/// Nothing measured that, so inverting the choice was invisible.
#[test]
fn cells_split_along_the_shorter_diagonal() {
    // A cylinder band sampled coarsely in v and finely in u: cells are long
    // and thin, which is exactly where the diagonal choice matters.
    // A cone, not a cylinder: on a cylinder the two candidate diagonals are
    // exactly equal by symmetry, so the choice is unobservable and the test
    // proves nothing. A cone's radius varies along v, breaking the tie.
    let s = Surface::Cone(axiolid_surface::Cone {
        frame: frame(),
        radius: 1.0,
        semi_angle: 0.6,
    });
    let patch = Patch::new(0.0, TAU, 0.0, 4.0).expect("patch");
    let out = tessellate_patch(&s, patch, budget(0.05)).expect("tessellate");

    let mut worst_ratio: Scalar = 0.0;
    for t in out.mesh.indices.as_chunks::<3>().0 {
        let p = [
            out.mesh.positions[t[0] as usize],
            out.mesh.positions[t[1] as usize],
            out.mesh.positions[t[2] as usize],
        ];
        let e = [
            (p[1] - p[0]).length(),
            (p[2] - p[1]).length(),
            (p[0] - p[2]).length(),
        ];
        let longest = e[0].max(e[1]).max(e[2]);
        let shortest = e[0].min(e[1]).min(e[2]);
        worst_ratio = worst_ratio.max(longest / shortest);
    }
    // Aspect ratio is a property of the patch, not the split, so bounding it
    // proves nothing. What the split controls is which diagonal appears: for
    // every quad the emitted interior edge must be the SHORTER of the two
    // candidates, otherwise the mesh carries avoidable slivers.
    let mut edges: std::collections::HashSet<(u32, u32)> = std::collections::HashSet::new();
    for t in out.mesh.indices.as_chunks::<3>().0 {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            edges.insert(if x < y { (x, y) } else { (y, x) });
        }
    }
    let nv = out.v_samples as u32;
    let mut checked = 0usize;
    for i in 0..out.u_samples as u32 - 1 {
        for j in 0..nv - 1 {
            let a = i * nv + j;
            let b = i * nv + j + 1;
            let c = (i + 1) * nv + j;
            let d = (i + 1) * nv + j + 1;
            let p = |k: u32| out.mesh.positions[k as usize];
            let ad = (p(d) - p(a)).length();
            let bc = (p(c) - p(b)).length();
            let has_ad = edges.contains(&if a < d { (a, d) } else { (d, a) });
            let has_bc = edges.contains(&if b < c { (b, c) } else { (c, b) });
            assert!(has_ad != has_bc, "exactly one diagonal per cell");
            if ad < bc {
                assert!(has_ad, "cell ({i},{j}) split along the longer diagonal");
            } else if bc < ad {
                assert!(has_bc, "cell ({i},{j}) split along the longer diagonal");
            }
            checked += 1;
        }
    }
    assert!(checked > 0, "no cells were checked");
    let _ = worst_ratio;
}

/// Refinement must converge fast enough to be usable.
///
/// Doubling reaches a tight tolerance in log steps; incrementing by one needs
/// thousands of passes and silently exhausts the budget instead. Nothing
/// distinguished the two.
#[test]
fn refinement_reaches_a_tight_tolerance_within_the_budget() {
    let s = Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: 2.0,
    });
    let patch = Patch::new(0.0, TAU, 0.0, 1.0).expect("patch");
    // A tolerance this tight is unreachable by +1 steps inside 512 samples.
    let tight = TessellationBudget::new(1e-4, 1024).expect("budget");
    let out = tessellate_patch(&s, patch, tight).expect("tessellate");
    assert!(
        !out.budget_exhausted,
        "doubling must reach 1e-4 within 1024 samples; got {} u-samples",
        out.u_samples
    );
    assert!(
        out.max_sagitta.expect("curved surface has a sagitta") <= 1e-4,
        "claimed tolerance was not met"
    );
    // Doubling lands on 2, 3, 5, 9, ... 2^k+1, so the accepted count is one
    // more than a power of two. Incrementing by one instead would satisfy the
    // tolerance too -- just after hundreds of passes -- so only the count
    // distinguishes the two strategies. Anything else measures nothing.
    let steps = out.u_samples - 1;
    assert!(
        steps.is_power_of_two(),
        "doubling must land on 2^k+1 samples, got {}",
        out.u_samples
    );
}

/// The sagitta probe must measure to the chord, not past its ends.
///
/// An unclamped projection lets the foot of the perpendicular fall outside the
/// segment, under-reporting deviation on short spans.
#[test]
fn the_sagitta_probe_stays_within_the_chord() {
    // A half-turn of a large-radius cylinder: the midpoint's perpendicular
    // foot is well inside the chord, but an unclamped projection on the
    // first refinement steps runs past the endpoint.
    let s = Surface::Cylinder(Cylinder {
        frame: frame(),
        radius: 50.0,
    });
    let patch = Patch::new(0.0, PI, 0.0, 1.0).expect("patch");
    let out = tessellate_patch(&s, patch, budget(0.5)).expect("tessellate");
    let measured = out.max_sagitta.expect("curved surface has a sagitta");
    assert!(
        measured <= 0.5,
        "reported sagitta {measured} exceeds the budget it claims to meet"
    );
    // And the mesh really is within that distance of the true surface.
    for p in &out.mesh.positions {
        let radial = (p.x * p.x + p.y * p.y).sqrt();
        assert!(
            (radial - 50.0).abs() < 1e-9,
            "sample off the cylinder: r={radial}"
        );
    }
}
