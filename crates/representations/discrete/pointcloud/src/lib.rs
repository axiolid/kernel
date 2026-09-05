//! Point-sampled geometry: an ordered set of 3D points with optional
//! per-point channels.
//!
//! # What this is
//!
//! The value a laser scan or photogrammetry capture becomes once it is
//! inside the kernel. Points carry no topology and no adjacency: a
//! pointcloud is a *sample* of a surface, not a description of one, and
//! pretending otherwise is where scan pipelines usually go wrong.
//!
//! # What this is not
//!
//! Not a file format. LAS, LAZ, E57, PCD and COPC parsing lives outside
//! `crates/`, exactly as STEP and IFC do for B-rep. No wire, vendor, or
//! source-format type may appear here — that boundary is what keeps the
//! kernel portable, and it is recorded in the ingestion-boundary ADR.
//!
//! Not an algorithm. Queries live in `axiolid-spatial`, reconstruction
//! behind the reconstruction contract. This crate owns the value and its
//! validation, nothing else.
//!
//! # Attributes are optional and validated together
//!
//! A scan may carry normals, colour, or intensity — or none of them. Each
//! channel is stored separately so a cloud with no extra data pays nothing,
//! and each is required to have exactly one entry per point. A channel of
//! the wrong length is a construction error rather than a hazard discovered
//! later by an algorithm indexing off the end.

#![forbid(unsafe_code)]

use axiolid_core::{Point3, Scalar, Vec3};

/// Why a pointcloud could not be constructed.
#[derive(Debug, Clone, PartialEq)]
#[non_exhaustive]
pub enum PointCloudError {
    /// A point has a non-finite coordinate.
    ///
    /// Carries the index so a caller can find the offending sample rather
    /// than re-scanning the input to locate it.
    NonFinitePoint {
        /// Index of the offending point.
        index: usize,
        /// The value as supplied.
        point: Point3,
    },
    /// A channel's length does not match the point count.
    ChannelLengthMismatch {
        /// Name of the offending channel.
        channel: &'static str,
        /// Entries the channel supplied.
        found: usize,
        /// Entries the cloud requires.
        expected: usize,
    },
    /// A normal is not finite, or is too short to carry a direction.
    ///
    /// A zero-length normal is refused rather than normalised to an
    /// arbitrary axis: the direction is genuinely unknown, and inventing
    /// one would be indistinguishable downstream from a measured value.
    DegenerateNormal {
        /// Index of the offending normal.
        index: usize,
        /// The value as supplied.
        normal: Vec3,
    },
    /// An intensity is not finite.
    NonFiniteIntensity {
        /// Index of the offending sample.
        index: usize,
        /// The value as supplied.
        intensity: Scalar,
    },
}

impl core::fmt::Display for PointCloudError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NonFinitePoint { index, point } => {
                write!(f, "point {index} is not finite: {point:?}")
            }
            Self::ChannelLengthMismatch {
                channel,
                found,
                expected,
            } => write!(
                f,
                "channel `{channel}` has {found} entries but the cloud has {expected} points"
            ),
            Self::DegenerateNormal { index, normal } => {
                write!(f, "normal {index} is not a usable direction: {normal:?}")
            }
            Self::NonFiniteIntensity { index, intensity } => {
                write!(f, "intensity {index} is not finite: {intensity}")
            }
        }
    }
}

impl core::error::Error for PointCloudError {}

/// Per-point colour, in unsigned 8-bit channels.
///
/// Deliberately not a floating-point triple: capture hardware reports 8- or
/// 16-bit colour, and widening it here would imply a precision the sensor
/// never had.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Colour {
    /// Red channel.
    pub red: u8,
    /// Green channel.
    pub green: u8,
    /// Blue channel.
    pub blue: u8,
}

impl Colour {
    /// A colour from its three channels.
    pub const fn new(red: u8, green: u8, blue: u8) -> Self {
        Self { red, green, blue }
    }
}

/// An ordered set of 3D points with optional per-point channels.
///
/// Order is meaningful: it is the caller's, and every channel is indexed by
/// it. Operations that reorder points must say so, because a consumer
/// holding a parallel array outside the cloud would otherwise be silently
/// desynchronised.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct PointCloud {
    points: Vec<Point3>,
    normals: Option<Vec<Vec3>>,
    colours: Option<Vec<Colour>>,
    intensities: Option<Vec<Scalar>>,
}

impl PointCloud {
    /// A cloud of positions, with no attributes.
    ///
    /// # Errors
    ///
    /// Refuses a non-finite coordinate, naming the offending index.
    pub fn new(points: Vec<Point3>) -> Result<Self, PointCloudError> {
        for (index, point) in points.iter().enumerate() {
            if !point.is_finite() {
                return Err(PointCloudError::NonFinitePoint {
                    index,
                    point: *point,
                });
            }
        }
        Ok(Self {
            points,
            normals: None,
            colours: None,
            intensities: None,
        })
    }

    /// An empty cloud.
    ///
    /// Legitimate rather than exceptional: a query that filters everything
    /// out should return an empty cloud, not an error. Consumers that need
    /// points must check [`is_empty`](Self::is_empty) themselves.
    pub fn empty() -> Self {
        Self::default()
    }

    /// Attach per-point normals.
    ///
    /// # Errors
    ///
    /// Refuses a length mismatch, and a normal that is non-finite or too
    /// short to carry a direction.
    pub fn with_normals(mut self, normals: Vec<Vec3>) -> Result<Self, PointCloudError> {
        if normals.len() != self.points.len() {
            return Err(PointCloudError::ChannelLengthMismatch {
                channel: "normals",
                found: normals.len(),
                expected: self.points.len(),
            });
        }
        for (index, normal) in normals.iter().enumerate() {
            // The threshold is the smallest normal positive value rather
            // than a modelling tolerance: this is a direction, not a
            // length, so it is scale-free.
            if !normal.is_finite() || normal.length() <= Scalar::MIN_POSITIVE {
                return Err(PointCloudError::DegenerateNormal {
                    index,
                    normal: *normal,
                });
            }
        }
        self.normals = Some(normals);
        Ok(self)
    }

    /// Attach per-point colours.
    ///
    /// # Errors
    ///
    /// Refuses a length mismatch. Colour channels cannot be individually
    /// invalid, since every 8-bit value is meaningful.
    pub fn with_colours(mut self, colours: Vec<Colour>) -> Result<Self, PointCloudError> {
        if colours.len() != self.points.len() {
            return Err(PointCloudError::ChannelLengthMismatch {
                channel: "colours",
                found: colours.len(),
                expected: self.points.len(),
            });
        }
        self.colours = Some(colours);
        Ok(self)
    }

    /// Attach per-point intensities.
    ///
    /// # Errors
    ///
    /// Refuses a length mismatch and a non-finite intensity. The range is
    /// not constrained: sensors report intensity on their own scale, and
    /// normalising it here would discard information the caller may need.
    pub fn with_intensities(mut self, intensities: Vec<Scalar>) -> Result<Self, PointCloudError> {
        if intensities.len() != self.points.len() {
            return Err(PointCloudError::ChannelLengthMismatch {
                channel: "intensities",
                found: intensities.len(),
                expected: self.points.len(),
            });
        }
        for (index, intensity) in intensities.iter().enumerate() {
            if !intensity.is_finite() {
                return Err(PointCloudError::NonFiniteIntensity {
                    index,
                    intensity: *intensity,
                });
            }
        }
        self.intensities = Some(intensities);
        Ok(self)
    }

    /// The points, in the caller's order.
    pub fn points(&self) -> &[Point3] {
        &self.points
    }

    /// Per-point normals, if the cloud carries them.
    pub fn normals(&self) -> Option<&[Vec3]> {
        self.normals.as_deref()
    }

    /// Per-point colours, if the cloud carries them.
    pub fn colours(&self) -> Option<&[Colour]> {
        self.colours.as_deref()
    }

    /// Per-point intensities, if the cloud carries them.
    pub fn intensities(&self) -> Option<&[Scalar]> {
        self.intensities.as_deref()
    }

    /// Number of points.
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// Whether the cloud has no points.
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// Whether every point carries a normal.
    ///
    /// Reconstruction algorithms differ sharply on this: Poisson needs
    /// oriented normals, ball pivoting does not. Asking is cheaper than
    /// discovering it mid-operation.
    pub fn has_normals(&self) -> bool {
        self.normals.is_some()
    }

    /// Axis-aligned bounds as `(min, max)`, or `None` when empty.
    ///
    /// Provided here because it needs no algorithm and every consumer wants
    /// it; anything beyond this belongs in a query crate.
    pub fn bounds(&self) -> Option<(Point3, Point3)> {
        let first = *self.points.first()?;
        let mut min = first;
        let mut max = first;
        for p in &self.points[1..] {
            min = Point3::new(min.x.min(p.x), min.y.min(p.y), min.z.min(p.z));
            max = Point3::new(max.x.max(p.x), max.y.max(p.y), max.z.max(p.z));
        }
        Some((min, max))
    }
}
