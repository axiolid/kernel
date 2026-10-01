// SPDX-License-Identifier: MPL-2.0
#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! Exact 3D Delaunay tetrahedralization.
//!
//! [`Delaunay3`] is an incremental Bowyer-Watson triangulation of a point
//! set in space. Its guarantees are stated and tested, not hoped for:
//!
//! - every tetrahedron is positively oriented, decided by the exact
//!   `orient3d`;
//! - no point lies strictly inside any tetrahedron's circumsphere, decided
//!   by the exact `insphere`;
//! - the boundary faces are exactly the convex hull, kept with an infinite
//!   vertex rather than a bounding box, so nothing outside the input is
//!   inserted;
//! - degenerate input -- cospherical points such as a cubic grid, points in
//!   the plane of a hull face, exact duplicates -- is handled by a symbolic
//!   perturbation (simulation of simplicity) that depends only on the
//!   coordinates, so the triangulation is the same for every insertion
//!   order.
//!
//! # What is refused
//!
//! - Points spanning fewer than three dimensions have no tetrahedra:
//!   [`Delaunay3::from_points`] reports [`Delaunay3Error::NotFullDimensional`]
//!   (a coplanar set belongs to `axiolid-triangulate`).
//! - Non-finite coordinates, and non-zero coordinates whose magnitude lies
//!   outside `[2^-100, 2^100]`, where the degree-5 and degree-6 exact
//!   predicates would leave binary64's normal range.
//! - Weighted (regular) triangulation, and constraints: neither input
//!   segments nor input faces are recovered (#126).
//!
//! # Example
//!
//! ```
//! use axiolid_core::Point3;
//! use axiolid_tetrahedralize::Delaunay3;
//!
//! // The eight corners of a cube are cospherical: several Delaunay
//! // tetrahedralizations exist, and the perturbation picks one of them.
//! let corners: Vec<Point3> = (0..8)
//!     .map(|i| Point3::new(f64::from(i & 1), f64::from((i >> 1) & 1), f64::from(i >> 2)))
//!     .collect();
//! let delaunay = Delaunay3::from_points(&corners)?;
//! let mesh = delaunay.tetrahedra();
//! assert!(mesh.tetrahedra.len() == 5 || mesh.tetrahedra.len() == 6);
//! assert_eq!(delaunay.hull_triangles().len(), 12);
//! # Ok::<(), axiolid_tetrahedralize::Delaunay3Error>(())
//! ```

mod delaunay;
mod hilbert;
mod sos;

pub use delaunay::{
    Delaunay3, Delaunay3Error, Insertion, Tetrahedra, MAX_COORDINATE, MIN_COORDINATE,
};
