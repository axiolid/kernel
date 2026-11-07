//! Capability handshake for the feature set compiled into the facade.

use axiolid_contracts::{
    BoundaryContract, IntegrationDescriptor, IntegrationProfile, Representation,
};

macro_rules! compiled_feature {
    ($target:expr, $name:literal) => {
        if cfg!(feature = $name) {
            $target.push($name);
        }
    };
}

macro_rules! compiled_representation {
    ($target:expr, $feature:literal, $representation:expr) => {
        if cfg!(feature = $feature) {
            $target.push($representation);
        }
    };
}

/// Describe the representations and features compiled into this facade build.
///
/// Operation features expose portable request/response traits; they are not
/// executable capability claims. A later provider-registration step adds a
/// `CapabilityDescriptor` only after a concrete operation provider is linked.
pub fn descriptor() -> IntegrationDescriptor {
    let mut descriptor = IntegrationDescriptor::empty(IntegrationProfile::RustFacade);
    descriptor.boundary = BoundaryContract::rust();
    descriptor.enabled_features = enabled_features();
    descriptor.representations = enabled_representations();
    descriptor
}

fn enabled_features() -> Vec<&'static str> {
    let mut features = Vec::new();
    compiled_feature!(features, "integration");
    compiled_feature!(features, "mesh");
    compiled_feature!(features, "linear");
    compiled_feature!(features, "predicates");
    compiled_feature!(features, "linear-intersection");
    compiled_feature!(features, "profiles");
    compiled_feature!(features, "curves");
    compiled_feature!(features, "surfaces");
    compiled_feature!(features, "topology");
    compiled_feature!(features, "brep");
    compiled_feature!(features, "primitives");
    compiled_feature!(features, "model");
    compiled_feature!(features, "evaluate");
    compiled_feature!(features, "nurbs");
    compiled_feature!(features, "tessellation");
    compiled_feature!(features, "spatial");
    compiled_feature!(features, "ray-mesh");
    compiled_feature!(features, "measure");
    compiled_feature!(features, "overlay");
    compiled_feature!(features, "field");
    compiled_feature!(features, "field-ops");
    compiled_feature!(features, "field-navigation");
    compiled_feature!(features, "heal");
    compiled_feature!(features, "contracts");
    compiled_feature!(features, "mesh-contracts");
    compiled_feature!(features, "mesh-boolean");
    compiled_feature!(features, "mesh-section");
    compiled_feature!(features, "mesh-hidden-lines");
    compiled_feature!(features, "graph-compile");
    compiled_feature!(features, "dispatch-mesh-boolean");
    compiled_feature!(features, "dispatch-mesh-section");
    compiled_feature!(features, "dispatch-mesh-hidden-lines");
    compiled_feature!(features, "hidden-lines-provider");
    compiled_feature!(features, "generate");
    compiled_feature!(features, "cpu");
    compiled_feature!(features, "parallel");
    compiled_feature!(features, "simd");
    compiled_feature!(features, "gpu");
    compiled_feature!(features, "discrete");
    compiled_feature!(features, "parametric");
    compiled_feature!(features, "advanced");
    compiled_feature!(features, "application");
    compiled_feature!(features, "portable-provider");
    compiled_feature!(features, "full");
    features.sort_unstable();
    features
}

fn enabled_representations() -> Vec<Representation> {
    let mut representations = vec![Representation::Scalar];
    compiled_representation!(representations, "linear", Representation::Linear);
    compiled_representation!(representations, "profiles", Representation::Profile2d);
    compiled_representation!(representations, "curves", Representation::AnalyticCurve);
    compiled_representation!(representations, "surfaces", Representation::AnalyticSurface);
    compiled_representation!(representations, "topology", Representation::Topology);
    compiled_representation!(representations, "brep", Representation::ExactBrep);
    compiled_representation!(representations, "mesh", Representation::TriangleMesh);
    compiled_representation!(representations, "model", Representation::ModelGraph);
    compiled_representation!(representations, "field", Representation::SampledField);
    representations.sort_unstable();
    representations.dedup();
    representations
}

#[cfg(test)]
mod tests {
    use super::*;
    use axiolid_contracts::{INTEGRATION_API_VERSION, MINIMUM_RUST_VERSION};

    #[test]
    fn compiled_descriptor_is_honest_about_default_surface() {
        let descriptor = descriptor();
        assert_eq!(descriptor.api_version, INTEGRATION_API_VERSION);
        assert_eq!(descriptor.minimum_rust_version, Some(MINIMUM_RUST_VERSION));
        assert!(descriptor.enabled_features.contains(&"integration"));
        assert!(descriptor.representations.contains(&Representation::Scalar));
        assert!(descriptor.capabilities.is_empty());
    }

    #[test]
    fn feature_and_representation_rows_are_unique() {
        let descriptor = descriptor();
        let mut features = descriptor.enabled_features.clone();
        features.sort_unstable();
        features.dedup();
        assert_eq!(features, descriptor.enabled_features);

        let mut representations = descriptor.representations.clone();
        representations.sort_unstable();
        representations.dedup();
        assert_eq!(representations, descriptor.representations);
    }
}
