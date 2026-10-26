//! Degeneracy-controlled scene generation for benchmarking and testing.
//!
//! The question a predicate suite must answer is not "how fast is it on random
//! data" -- random points in a box are never degenerate, so that measures only
//! the filter's fast path. The useful question is how throughput and the
//! escalation rate move as the fraction of degenerate inputs rises.
//!
//! Scenes are deterministic given a seed, so a measurement is reproducible.

use axiolid_core::{Point2, Point3};

/// Fraction of inputs that are exactly degenerate.
///
/// The tiers span the range that matters: clean data, the incidental
/// degeneracies of authored CAD models, the systematic ones of grid-aligned
/// architecture, and a stress tier well past anything real.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DegeneracyRate {
    /// No degenerate inputs.
    None,
    /// One in ten thousand.
    Rare,
    /// One in a hundred.
    Occasional,
    /// One in ten.
    Frequent,
}

impl DegeneracyRate {
    /// Every tier, ascending.
    pub const ALL: [Self; 4] = [Self::None, Self::Rare, Self::Occasional, Self::Frequent];

    /// The fraction as a ratio.
    #[must_use]
    pub const fn fraction(self) -> f64 {
        match self {
            Self::None => 0.0,
            Self::Rare => 0.0001,
            Self::Occasional => 0.01,
            Self::Frequent => 0.1,
        }
    }

    /// Human-readable label for reports.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::None => "0%",
            Self::Rare => "0.01%",
            Self::Occasional => "1%",
            Self::Frequent => "10%",
        }
    }

    /// Whether the sample at `index` should be degenerate.
    ///
    /// Deterministic striping rather than a random draw, so a given index is
    /// degenerate in every run and a regression is reproducible.
    #[must_use]
    pub fn is_degenerate(self, index: usize) -> bool {
        match self {
            Self::None => false,
            Self::Rare => index.is_multiple_of(10_000),
            Self::Occasional => index.is_multiple_of(100),
            Self::Frequent => index.is_multiple_of(10),
        }
    }
}

/// Deterministic xorshift, so scenes are reproducible across runs and machines.
#[derive(Debug, Clone)]
pub struct SceneRng {
    state: u64,
}

impl SceneRng {
    /// Seed the generator. Zero is remapped: xorshift is stuck at zero.
    #[must_use]
    pub const fn new(seed: u64) -> Self {
        Self {
            state: if seed == 0 {
                0x2545_F491_4F6C_DD1D
            } else {
                seed
            },
        }
    }

    /// Next raw value.
    pub fn next_u64(&mut self) -> u64 {
        self.state ^= self.state << 13;
        self.state ^= self.state >> 7;
        self.state ^= self.state << 17;
        self.state
    }

    /// Integer coordinate in `[-bound, bound)`.
    ///
    /// Integers, not arbitrary floats: they keep the degenerate constructions
    /// below exactly degenerate, which is the entire point of the scene.
    pub fn coordinate(&mut self, bound: i64) -> f64 {
        (((self.next_u64() >> 33) as i64) % (2 * bound) - bound) as f64
    }
}

/// A 2D orientation query: three points.
pub type Orient2Case = [Point2; 3];

/// A 3D orientation query: four points.
pub type Orient3Case = [Point3; 4];

/// Generate `count` orientation queries with the given degeneracy rate.
///
/// Degenerate cases are exactly collinear by construction (`c` is an integer
/// affine combination of `a` and `b`), so they are not merely close to the
/// boundary -- they are on it, which is what forces the exact path.
#[must_use]
pub fn orient2_scene(count: usize, rate: DegeneracyRate, seed: u64) -> Vec<Orient2Case> {
    let mut rng = SceneRng::new(seed);
    (0..count)
        .map(|index| {
            let a = Point2::new(rng.coordinate(1_000), rng.coordinate(1_000));
            let b = Point2::new(rng.coordinate(1_000), rng.coordinate(1_000));
            let c = if rate.is_degenerate(index) {
                let k = ((rng.next_u64() >> 40) as i64 % 7) - 3;
                Point2::new(
                    a.x + (k as f64) * (b.x - a.x),
                    a.y + (k as f64) * (b.y - a.y),
                )
            } else {
                Point2::new(rng.coordinate(1_000), rng.coordinate(1_000))
            };
            [a, b, c]
        })
        .collect()
}

/// Generate `count` 3D orientation queries with the given degeneracy rate.
///
/// Degenerate cases are exactly coplanar: `d = a + i*(b-a) + j*(c-a)` for
/// integers `i, j`, so no rounding can move the point off the plane.
#[must_use]
pub fn orient3_scene(count: usize, rate: DegeneracyRate, seed: u64) -> Vec<Orient3Case> {
    let mut rng = SceneRng::new(seed);
    (0..count)
        .map(|index| {
            let p = |r: &mut SceneRng| {
                Point3::new(
                    r.coordinate(1_000),
                    r.coordinate(1_000),
                    r.coordinate(1_000),
                )
            };
            let (a, b, c) = (p(&mut rng), p(&mut rng), p(&mut rng));
            let d = if rate.is_degenerate(index) {
                let i = (((rng.next_u64() >> 40) as i64 % 5) - 2) as f64;
                let j = (((rng.next_u64() >> 40) as i64 % 5) - 2) as f64;
                Point3::new(
                    a.x + i * (b.x - a.x) + j * (c.x - a.x),
                    a.y + i * (b.y - a.y) + j * (c.y - a.y),
                    a.z + i * (b.z - a.z) + j * (c.z - a.z),
                )
            } else {
                p(&mut rng)
            };
            [a, b, c, d]
        })
        .collect()
}

/// Generate 3D orientation queries that are *nearly* coplanar.
///
/// Distinct from [`orient3_scene`]'s degenerate cases, which lie exactly on the
/// plane and are therefore answered zero. These sit one unit-in-the-last-place
/// off it: the filter must defer, and the exact path must then return a
/// non-zero sign. A suite containing only exact degeneracies cannot tell a
/// correct exact path from one that always answers zero.
#[must_use]
pub fn near_coplanar_scene(count: usize, seed: u64) -> Vec<Orient3Case> {
    let mut rng = SceneRng::new(seed);
    (0..count)
        .map(|_| {
            let p = |r: &mut SceneRng| {
                Point3::new(
                    r.coordinate(1_000),
                    r.coordinate(1_000),
                    r.coordinate(1_000),
                )
            };
            let (a, b, c) = (p(&mut rng), p(&mut rng), p(&mut rng));
            let i = (((rng.next_u64() >> 40) as i64 % 5) - 2) as f64;
            let j = (((rng.next_u64() >> 40) as i64 % 5) - 2) as f64;
            let z = a.z + i * (b.z - a.z) + j * (c.z - a.z);
            [
                a,
                b,
                c,
                Point3::new(
                    a.x + i * (b.x - a.x) + j * (c.x - a.x),
                    a.y + i * (b.y - a.y) + j * (c.y - a.y),
                    // One ULP off the plane: far too small for the filter to
                    // resolve, but a definite side that exact arithmetic must
                    // recover.
                    next_ulp(z),
                ),
            ]
        })
        .collect()
}

/// The next representable f64 above `value`.
///
/// Bit increment rather than `next_after`, which core does not expose. The
/// step is exactly one ULP, the smallest perturbation that still moves the
/// point off the plane.
#[must_use]
fn next_ulp(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    let bits = value.to_bits();
    f64::from_bits(if value >= 0.0 { bits + 1 } else { bits - 1 })
}
