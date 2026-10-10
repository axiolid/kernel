//! Exact conic curve data.

use axiolid_core::{Frame2, Frame3, Scalar};

/// Circle in a two-dimensional frame.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Circle2 {
    /// Local frame.
    pub frame: Frame2,
    /// Radius. Validation rejects non-positive values.
    pub radius: Scalar,
}

/// Circle in a three-dimensional plane frame.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Circle3 {
    /// Local frame; `z` is the plane normal.
    pub frame: Frame3,
    /// Radius. Validation rejects non-positive values.
    pub radius: Scalar,
}

/// Ellipse in a two-dimensional frame.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Ellipse2 {
    /// Local frame.
    pub frame: Frame2,
    /// Semi-axis along local x.
    pub semi_axis_x: Scalar,
    /// Semi-axis along local y.
    pub semi_axis_y: Scalar,
}

/// Ellipse in a three-dimensional plane frame.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct Ellipse3 {
    /// Local frame; `z` is the plane normal.
    pub frame: Frame3,
    /// Semi-axis along local x.
    pub semi_axis_x: Scalar,
    /// Semi-axis along local y.
    pub semi_axis_y: Scalar,
}
