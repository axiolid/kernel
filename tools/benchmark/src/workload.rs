//! Deterministic workload generation.
//!
//! Inputs must be reproducible byte-for-byte across runs and machines, or an
//! instruction-count comparison measures the input generator rather than the
//! code under test. Nothing here consults the clock, the allocator, or a
//! randomly-seeded hasher.

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_mesh::TriMesh;

/// A counter-based pseudo-random generator with an explicit seed.
///
/// This is SplitMix64. It is chosen over `std`'s hasher-backed randomness for
/// one reason: `std::collections::HashMap`'s `RandomState` is seeded per
/// instance, so anything derived from it differs between runs. A benchmark
/// built on that cannot distinguish a regression from a reseed.
#[derive(Debug, Clone)]
pub struct Rng {
    state: u64,
}

impl Rng {
    /// Start a stream from an explicit seed.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    /// Next raw 64-bit value.
    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Next value in `[0, 1)`.
    ///
    /// Uses the top 53 bits, which is the exact mantissa width of `f64`, so
    /// every representable value in the range is reachable and none is favoured
    /// by truncation.
    pub fn next_unit(&mut self) -> Scalar {
        const SCALE: Scalar = 1.0 / (1u64 << 53) as Scalar;
        (self.next_u64() >> 11) as Scalar * SCALE
    }

    /// Next value in `[low, high)`.
    pub fn next_range(&mut self, low: Scalar, high: Scalar) -> Scalar {
        low + (high - low) * self.next_unit()
    }
}

/// Problem size for a scaling study.
///
/// Named sizes rather than raw counts, so a scaling curve is reported against a
/// stable axis even if the underlying element count of a workload is retuned.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Scale {
    /// Smallest useful size; dominated by fixed overhead.
    Tiny,
    /// Small but past the fixed-cost floor.
    Small,
    /// The size most microbenchmarks should report.
    Medium,
    /// Large enough to leave cache on most machines.
    Large,
}

impl Scale {
    /// Every scale, ascending. The axis of a scaling study.
    #[must_use]
    pub const fn all() -> [Self; 4] {
        [Self::Tiny, Self::Small, Self::Medium, Self::Large]
    }

    /// Stable label for reporting.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Tiny => "tiny",
            Self::Small => "small",
            Self::Medium => "medium",
            Self::Large => "large",
        }
    }

    /// Element count for this scale.
    ///
    /// Powers of four so successive scales are far enough apart to separate a
    /// linear curve from a quadratic one over the four points.
    #[must_use]
    pub const fn count(self) -> usize {
        match self {
            Self::Tiny => 64,
            Self::Small => 256,
            Self::Medium => 1_024,
            Self::Large => 4_096,
        }
    }
}

/// A generated input paired with what it is known to contain.
///
/// The `expected_*` fields are derived from the *construction*, never measured
/// from the produced value, so they are an independent oracle rather than a
/// restatement of whatever the code produced.
#[derive(Debug, Clone)]
pub struct Workload {
    /// Stable name for reporting.
    pub name: &'static str,
    /// The mesh under test.
    pub mesh: TriMesh,
    /// Signed volume the construction guarantees, when one is known.
    pub expected_volume: Option<Scalar>,
    /// Closed-surface Euler characteristic, when the construction guarantees one.
    pub expected_euler: Option<i64>,
}

/// An axis-aligned box of the given extents, centred on the origin.
///
/// Volume is known exactly from the extents, and a box is a closed genus-0
/// shell, so `V - E + F` is 2. Both are stated by construction.
#[must_use]
pub fn box_mesh(size_x: Scalar, size_y: Scalar, size_z: Scalar) -> Workload {
    let (hx, hy, hz) = (size_x / 2.0, size_y / 2.0, size_z / 2.0);
    let positions = vec![
        Point3::new(-hx, -hy, -hz),
        Point3::new(hx, -hy, -hz),
        Point3::new(hx, hy, -hz),
        Point3::new(-hx, hy, -hz),
        Point3::new(-hx, -hy, hz),
        Point3::new(hx, -hy, hz),
        Point3::new(hx, hy, hz),
        Point3::new(-hx, hy, hz),
    ];
    // Outward-facing winding, two triangles per quad.
    let indices = vec![
        0, 2, 1, 0, 3, 2, // -z
        4, 5, 6, 4, 6, 7, // +z
        0, 1, 5, 0, 5, 4, // -y
        3, 7, 6, 3, 6, 2, // +y
        0, 4, 7, 0, 7, 3, // -x
        1, 2, 6, 1, 6, 5, // +x
    ];
    Workload {
        name: "box",
        mesh: TriMesh::new(positions, indices),
        expected_volume: Some(size_x * size_y * size_z),
        expected_euler: Some(2),
    }
}

/// A closed triangulated sphere with `count` roughly-uniform surface samples.
///
/// Built as a UV sphere so the triangulation is deterministic and the shell is
/// closed by construction. The volume of the *polyhedron* is not the volume of
/// the ideal sphere, so no volume oracle is claimed here — only the Euler
/// characteristic, which the construction does guarantee.
#[must_use]
pub fn sphere_mesh(radius: Scalar, count: usize) -> Workload {
    let rings = (count as Scalar).sqrt().max(3.0) as usize;
    let segments = rings.max(3);
    let mut positions = Vec::with_capacity(rings * segments + 2);

    positions.push(Point3::new(0.0, 0.0, radius));
    for ring in 1..rings {
        let theta = core::f64::consts::PI * ring as Scalar / rings as Scalar;
        let (sin_t, cos_t) = theta.sin_cos();
        for segment in 0..segments {
            let phi = core::f64::consts::TAU * segment as Scalar / segments as Scalar;
            let (sin_p, cos_p) = phi.sin_cos();
            positions.push(Point3::new(
                radius * sin_t * cos_p,
                radius * sin_t * sin_p,
                radius * cos_t,
            ));
        }
    }
    positions.push(Point3::new(0.0, 0.0, -radius));

    let south = positions.len() as u32 - 1;
    let mut indices = Vec::new();
    // North cap.
    for segment in 0..segments {
        let next = (segment + 1) % segments;
        indices.extend_from_slice(&[0, 1 + segment as u32, 1 + next as u32]);
    }
    // Quad bands between successive rings.
    for ring in 0..rings.saturating_sub(2) {
        let base = 1 + ring * segments;
        let next_base = base + segments;
        for segment in 0..segments {
            let next = (segment + 1) % segments;
            let (a, b) = ((base + segment) as u32, (base + next) as u32);
            let (c, d) = ((next_base + segment) as u32, (next_base + next) as u32);
            indices.extend_from_slice(&[a, c, d, a, d, b]);
        }
    }
    // South cap.
    let last_base = 1 + rings.saturating_sub(2) * segments;
    for segment in 0..segments {
        let next = (segment + 1) % segments;
        indices.extend_from_slice(&[
            (last_base + next) as u32,
            (last_base + segment) as u32,
            south,
        ]);
    }

    Workload {
        name: "sphere",
        mesh: TriMesh::new(positions, indices),
        expected_volume: None,
        expected_euler: Some(2),
    }
}

/// `count` deterministic points inside the axis-aligned box of the given half-extent.
///
/// For spatial-index workloads, where the geometry is a point set rather than a
/// shell.
#[must_use]
pub fn point_cloud(seed: u64, count: usize, half_extent: Scalar) -> Vec<Point3> {
    let mut rng = Rng::new(seed);
    (0..count)
        .map(|_| {
            Point3::new(
                rng.next_range(-half_extent, half_extent),
                rng.next_range(-half_extent, half_extent),
                rng.next_range(-half_extent, half_extent),
            )
        })
        .collect()
}

/// A wall with `openings` evenly spaced rectangular voids, as separate cutters.
///
/// The shape of a real subtraction workload rather than a synthetic one: a
/// large thin subject and many small disjoint tools. Cutters are disjoint and
/// each passes fully through the wall, so the removed volume is exactly the
/// sum of the cutter volumes -- a *derived* ground truth, not a recording of
/// what the kernel happened to return.
///
/// Deterministic: positions come from the opening index, never from a clock
/// or a hash seed.
#[must_use]
pub fn wall_with_openings(openings: usize) -> WallWorkload {
    // A wall long enough that every opening keeps a solid margin around it.
    let length = 2.0 * openings as Scalar + 2.0;
    let height = 3.0;
    let thickness = 0.4;
    let subject = box_mesh(length, thickness, height).mesh;

    // Each opening is a cuboid that pierces the wall completely, so the
    // subtraction removes its full volume rather than a clipped part of it.
    let opening_width = 1.0;
    let opening_height = 1.5;
    let depth = thickness * 2.0;
    let mut tools = Vec::with_capacity(openings);
    for index in 0..openings {
        let centre = -length / 2.0 + 2.0 * index as Scalar + 2.0;
        let mut cutter = box_mesh(opening_width, depth, opening_height).mesh;
        translate(&mut cutter, Vec3::new(centre, 0.0, 0.0));
        tools.push(cutter);
    }

    let solid = length * thickness * height;
    let removed = opening_width * thickness * opening_height * openings as Scalar;
    WallWorkload {
        subject,
        tools,
        expected_volume: solid - removed,
    }
}

/// A subtraction workload with a ground truth derived from its construction.
pub struct WallWorkload {
    /// The wall being cut.
    pub subject: TriMesh,
    /// Disjoint cutters, each piercing the wall completely.
    pub tools: Vec<TriMesh>,
    /// Solid volume minus the cutter volumes, computed from the dimensions
    /// rather than measured from a result.
    pub expected_volume: Scalar,
}

/// Shift every vertex of `mesh` by `offset`.
fn translate(mesh: &mut TriMesh, offset: Vec3) {
    for position in &mut mesh.positions {
        *position += offset;
    }
}
