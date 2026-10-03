//! Closed conversions from built-in representations into graph nodes.
//!
//! Backend/provider traits remain intentionally open. This trait is sealed
//! because `GeometryNode` is a closed built-in vocabulary: downstream types
//! must first translate into one of these neutral representations.

use axiolid_core::{Aabb, Frame2, Frame3, Transform3};
use axiolid_curve::{Curve2, Curve3};
use axiolid_mesh::{PolygonMesh, TriMesh};
use axiolid_primitive::{HalfSpace, Primitive};
use axiolid_profile::Profile;
use axiolid_surface::Surface;
use axiolid_topology::BRep;

use crate::{
    CurveRelation, CurveStation, GeometryNode, Instance, NodeId, OpenProfile, PointOnCurve,
    PointOnSurface, SolidOperation, SurfaceRelation,
};

mod private {
    pub trait Sealed {}
}

/// A built-in neutral value that can be inserted into a [`crate::GeometryGraph`].
///
/// This is sealed deliberately. Third-party kernels extend execution through
/// the open traits in operation-contract packages; source adapters extend input
/// by translating their values into one of the canonical node families.
///
/// ```compile_fail
/// use axiolid_model::{BuiltInNode, GeometryNode};
///
/// struct ForeignNode;
///
/// impl Into<GeometryNode> for ForeignNode {
///     fn into(self) -> GeometryNode {
///         panic!("compile-only example")
///     }
/// }
///
/// // Rejected because the private sealing trait cannot be implemented here.
/// impl BuiltInNode for ForeignNode {}
/// ```
pub trait BuiltInNode: private::Sealed + Into<GeometryNode> {}

macro_rules! built_in_node {
    ($value:ty, $variant:ident) => {
        impl From<$value> for GeometryNode {
            fn from(value: $value) -> Self {
                Self::$variant(value)
            }
        }

        impl private::Sealed for $value {}
        impl BuiltInNode for $value {}
    };
}

built_in_node!(Frame2, Frame2);
built_in_node!(Frame3, Frame3);
built_in_node!(Transform3, Transform);
built_in_node!(Curve2, Curve2);
built_in_node!(Curve3, Curve3);
built_in_node!(CurveRelation, CurveRelation);
built_in_node!(PointOnCurve, PointOnCurve);
built_in_node!(Surface, Surface);
built_in_node!(SurfaceRelation, SurfaceRelation);
built_in_node!(PointOnSurface, PointOnSurface);
built_in_node!(Profile, Profile);
built_in_node!(OpenProfile, OpenProfile);
built_in_node!(Primitive, Primitive);
built_in_node!(HalfSpace, HalfSpace);
built_in_node!(SolidOperation, SolidOperation);
built_in_node!(BRep<NodeId>, BRep);
built_in_node!(PolygonMesh, PolygonMesh);
built_in_node!(TriMesh, TriMesh);
built_in_node!(Aabb, BoundingBox);
built_in_node!(Instance, Instance);
built_in_node!(CurveStation, CurveStation);

#[cfg(test)]
mod tests {
    use axiolid_primitive::Primitive;

    use super::*;

    fn accepts_built_in<T: BuiltInNode>(value: T) -> GeometryNode {
        value.into()
    }

    #[test]
    fn built_in_values_convert_without_format_knowledge() {
        let node = accepts_built_in(Primitive::Sphere { radius: 2.0 });
        assert!(matches!(
            node,
            GeometryNode::Primitive(Primitive::Sphere { radius: 2.0 })
        ));
    }
}
