//! Scenes that binary64 rounding used to refuse as `Degenerate` (#410).
//!
//! Two seeded scenes of rotated boxes, in an oblique and an isometric view
//! with a near plane, frozen from the provider's own oracle tests. In each,
//! some region of the exact arrangement has vertices that round to a single
//! point or an edge that passes through the rounding cell of another
//! vertex. They have no hand-derived rings: a provider must answer them
//! (ADR 0100), and the answer must show the snapping in its evidence.

use axiolid_core::{Frame3, Point3, Vec3};
use axiolid_mesh::TriMesh;

use super::FaceFixture;
use crate::fixtures::box_mesh;
use crate::{MeshRole, OrthographicView, VisibleFaceLimits};

/// A scene with no expected rings: only that it is answered and snapped.
pub(crate) struct SnapCase {
    pub name: &'static str,
    pub meshes: Vec<FaceFixture>,
    pub view: OrthographicView,
    pub limits: VisibleFaceLimits,
}

/// A frozen scene: the frame, the near plane and the corners of each box,
/// all as binary64 bit patterns.
struct SnapScene {
    origin: [u64; 3],
    axes: [[u64; 3]; 3],
    near: u64,
    boxes: &'static [[[u64; 3]; 8]],
}

/// Seed 1 of the oracle scenes, frozen as binary64 bit patterns
/// (the generator uses `sin_cos`, which is not portable).
const SEED_1_OBLIQUE: SnapScene = SnapScene {
    origin: [0xbfd999999999999a, 0x3fe6666666666666, 0x3fb999999999999a],
    axes: [
        [0x3feb7095010f9355, 0x3fe076bfcd6fbecd, 0x8000000000000000],
        [0xbfda9c103e6b7186, 0x3fe62cb8340433ef, 0x3fe2d9362c36c5be],
        [0x3fd3651a0eb63341, 0xbfe02995b6ed2ab6, 0x3fe9dc22be484457],
    ],
    near: 0x3fef98fa73b3804d,
    boxes: &[
        [
            [0x3fd8000000000000, 0x3ff2000000000000, 0x3fe0000000000000],
            [0x3fe8000000000000, 0x3ff2000000000000, 0x3fe0000000000000],
            [0x3fe8000000000000, 0x4001000000000000, 0x3fe0000000000000],
            [0x3fd8000000000000, 0x4001000000000000, 0x3fe0000000000000],
            [0x3fd8000000000000, 0x3ff2000000000000, 0x3ff6000000000000],
            [0x3fe8000000000000, 0x3ff2000000000000, 0x3ff6000000000000],
            [0x3fe8000000000000, 0x4001000000000000, 0x3ff6000000000000],
            [0x3fd8000000000000, 0x4001000000000000, 0x3ff6000000000000],
        ],
        [
            [0x3ff6000000000000, 0x3fd0000000000000, 0x3ff0000000000000],
            [0x3ff8000000000000, 0x3fd0000000000000, 0x3ff0000000000000],
            [0x3ff8000000000000, 0x3fe4000000000000, 0x3ff0000000000000],
            [0x3ff6000000000000, 0x3fe4000000000000, 0x3ff0000000000000],
            [0x3ff6000000000000, 0x3fd0000000000000, 0x3ffc000000000000],
            [0x3ff8000000000000, 0x3fd0000000000000, 0x3ffc000000000000],
            [0x3ff8000000000000, 0x3fe4000000000000, 0x3ffc000000000000],
            [0x3ff6000000000000, 0x3fe4000000000000, 0x3ffc000000000000],
        ],
        [
            [0x3ff248d225efef18, 0x3fc3dc6cbf64c634, 0x3f9f548b0f33ef40],
            [0x3ff98ea6560e04ee, 0x3fef9dce504907e2, 0xbfd241457b9c59d6],
            [0x3feb7ccfa9781721, 0x3ff89e54aeec7aa2, 0x3fbaabd5b63c3658],
            [0x3fd9e24e9277d6e6, 0x3fe695f63d691eee, 0x3fdae1839a1ea65e],
            [0x3ffa41982b43f470, 0x3fcb0d5a889c2af0, 0x3fecaa8549387935],
            [0x4000c3b62db10523, 0x3ff0b504e14b7089, 0x3fe28f3e32f0acd1],
            [0x3ff5b72dda1010e8, 0x3ff984726813673a, 0x3fef055ba7866086],
            [0x3fece2b353e3f623, 0x3fe86231afb6f81e, 0x3ff490515ee71676],
        ],
        [
            [0x3fe4000000000000, 0x0000000000000000, 0x3fe0000000000000],
            [0x3ff4000000000000, 0x0000000000000000, 0x3fe0000000000000],
            [0x3ff4000000000000, 0x3fd0000000000000, 0x3fe0000000000000],
            [0x3fe4000000000000, 0x3fd0000000000000, 0x3fe0000000000000],
            [0x3fe4000000000000, 0x0000000000000000, 0x3fe4000000000000],
            [0x3ff4000000000000, 0x0000000000000000, 0x3fe4000000000000],
            [0x3ff4000000000000, 0x3fd0000000000000, 0x3fe4000000000000],
            [0x3fe4000000000000, 0x3fd0000000000000, 0x3fe4000000000000],
        ],
        [
            [0x3fe00b9b6f609aa8, 0x3fec77fee3f89cae, 0x3fe0686ffdfaae10],
            [0x3fe30b3a9ba3da0a, 0x3ff9b7841c2c195b, 0x3fd43a30d6418a04],
            [0x3fd4e9d1f4eb2432, 0x3ffb563e0dee39ec, 0x3fe1d821d29b0ea7],
            [0x3fcdd52738c94ae0, 0x3fefb572c77cddd1, 0x3fe823796574f7b5],
            [0x3fe58b17058a6de7, 0x3fed5383e4238c27, 0x3fe627de2d64f159],
            [0x3fe88ab631cdad48, 0x3ffa25469c419118, 0x3fdfb90d35161096],
            [0x3fdfe8c9213ecab0, 0x3ffbc4008e03b1a9, 0x3fe79790020551f0],
            [0x3fd9e98ac8b84bed, 0x3ff0487be3d3e6a5, 0x3fede2e794df3afe],
        ],
    ],
};

/// Seed 9 of the oracle scenes, frozen as binary64 bit patterns
/// (the generator uses `sin_cos`, which is not portable).
const SEED_9_ISOMETRIC: SnapScene = SnapScene {
    origin: [0x3fb999999999999a, 0x3fc999999999999a, 0x3fd3333333333333],
    axes: [
        [0x3fe6a09e667f3bcd, 0xbfe6a09e667f3bcd, 0x0000000000000000],
        [0x3fda20bd700c2c3f, 0x3fda20bd700c2c3f, 0xbfea20bd700c2c3f],
        [0x3fe279a74590331d, 0x3fe279a74590331d, 0x3fe279a74590331d],
    ],
    near: 0x3ff3dc609dfb03bf,
    boxes: &[
        [
            [0x3fe0000000000000, 0x3fd8000000000000, 0x0000000000000000],
            [0x3ff6000000000000, 0x3fd8000000000000, 0x0000000000000000],
            [0x3ff6000000000000, 0x3ff6000000000000, 0x0000000000000000],
            [0x3fe0000000000000, 0x3ff6000000000000, 0x0000000000000000],
            [0x3fe0000000000000, 0x3fd8000000000000, 0x3fe4000000000000],
            [0x3ff6000000000000, 0x3fd8000000000000, 0x3fe4000000000000],
            [0x3ff6000000000000, 0x3ff6000000000000, 0x3fe4000000000000],
            [0x3fe0000000000000, 0x3ff6000000000000, 0x3fe4000000000000],
        ],
        [
            [0x3fe0000000000000, 0x3ff2000000000000, 0x3fe0000000000000],
            [0x3fec000000000000, 0x3ff2000000000000, 0x3fe0000000000000],
            [0x3fec000000000000, 0x4000000000000000, 0x3fe0000000000000],
            [0x3fe0000000000000, 0x4000000000000000, 0x3fe0000000000000],
            [0x3fe0000000000000, 0x3ff2000000000000, 0x3ff2000000000000],
            [0x3fec000000000000, 0x3ff2000000000000, 0x3ff2000000000000],
            [0x3fec000000000000, 0x4000000000000000, 0x3ff2000000000000],
            [0x3fe0000000000000, 0x4000000000000000, 0x3ff2000000000000],
        ],
        [
            [0x3fc228e97dfbc7fa, 0x3ff297713d1c6835, 0x3febcded320c3faf],
            [0x3fd5139da655e149, 0x3ff4683f7acd1f3a, 0x3fe7f3f6e440bb78],
            [0x3fd094bb45397e06, 0x3ff7ea75010e3e58, 0x3feb117e41662865],
            [0x3fb25649778602e8, 0x3ff619a6c35d8753, 0x3feeeb748f31ac9c],
            [0x3fe7b5a25d6340fd, 0x3ff0158afef1c1a8, 0x3ffa7740df4cebce],
            [0x3fedb536d10f3fa3, 0x3ff1e6593ca278ad, 0x3ff88a45b86729b2],
            [0x3feb75c5a0810e02, 0x3ff5688ec2e397cb, 0x3ffa190966f9e028],
            [0x3fe576312cd50f5c, 0x3ff397c08532e0c6, 0x3ffc06048ddfa244],
        ],
        [
            [0x3fd8000000000000, 0x3ff6000000000000, 0x0000000000000000],
            [0x3ff2000000000000, 0x3ff6000000000000, 0x0000000000000000],
            [0x3ff2000000000000, 0x3ffa000000000000, 0x0000000000000000],
            [0x3fd8000000000000, 0x3ffa000000000000, 0x0000000000000000],
            [0x3fd8000000000000, 0x3ff6000000000000, 0x3ff0000000000000],
            [0x3ff2000000000000, 0x3ff6000000000000, 0x3ff0000000000000],
            [0x3ff2000000000000, 0x3ffa000000000000, 0x3ff0000000000000],
            [0x3fd8000000000000, 0x3ffa000000000000, 0x3ff0000000000000],
        ],
        [
            [0x3ff85c9106996d72, 0x3fed08f77ffb729a, 0x3fec587aa904f514],
            [0x3ff3d43f445a715d, 0x3ff3ae5672a17cb5, 0x3fe4263a6eadeffe],
            [0x3fe43c01859b3f76, 0x3fea189f7a7fb3b3, 0x3fe8eada45498782],
            [0x3fed4ca50a1937a0, 0x3fdf89d42a7059c6, 0x3ff08e8d3fd0464b],
            [0x3ff5e1ff3d326045, 0x3ff6f3b042c02626, 0x3ffb8a92dd5b3c3f],
            [0x3ff159ad7af36430, 0x3ffc1d8af563e98e, 0x3ff77172c02fb9b5],
            [0x3fde8dbbe59a4a38, 0x3ff57b84400246b3, 0x3ff9d3c2ab7d8576],
            [0x3fe85781774b1d46, 0x3ff051a98d5e834b, 0x3ffdece2c8a90801],
        ],
        [
            [0x3fe4000000000000, 0x3fd0000000000000, 0x3fe0000000000000],
            [0x3fe8000000000000, 0x3fd0000000000000, 0x3fe0000000000000],
            [0x3fe8000000000000, 0x3fe8000000000000, 0x3fe0000000000000],
            [0x3fe4000000000000, 0x3fe8000000000000, 0x3fe0000000000000],
            [0x3fe4000000000000, 0x3fd0000000000000, 0x3ff0000000000000],
            [0x3fe8000000000000, 0x3fd0000000000000, 0x3ff0000000000000],
            [0x3fe8000000000000, 0x3fe8000000000000, 0x3ff0000000000000],
            [0x3fe4000000000000, 0x3fe8000000000000, 0x3ff0000000000000],
        ],
        [
            [0x3ff8000000000000, 0x3ff0000000000000, 0x3fe4000000000000],
            [0x3ffc000000000000, 0x3ff0000000000000, 0x3fe4000000000000],
            [0x3ffc000000000000, 0x3ff6000000000000, 0x3fe4000000000000],
            [0x3ff8000000000000, 0x3ff6000000000000, 0x3fe4000000000000],
            [0x3ff8000000000000, 0x3ff0000000000000, 0x3ff8000000000000],
            [0x3ffc000000000000, 0x3ff0000000000000, 0x3ff8000000000000],
            [0x3ffc000000000000, 0x3ff6000000000000, 0x3ff8000000000000],
            [0x3ff8000000000000, 0x3ff6000000000000, 0x3ff8000000000000],
        ],
        [
            [0x3fe8000000000000, 0x3ff2000000000000, 0x3fe8000000000000],
            [0x3fec000000000000, 0x3ff2000000000000, 0x3fe8000000000000],
            [0x3fec000000000000, 0x3ff6000000000000, 0x3fe8000000000000],
            [0x3fe8000000000000, 0x3ff6000000000000, 0x3fe8000000000000],
            [0x3fe8000000000000, 0x3ff2000000000000, 0x3ff4000000000000],
            [0x3fec000000000000, 0x3ff2000000000000, 0x3ff4000000000000],
            [0x3fec000000000000, 0x3ff6000000000000, 0x3ff4000000000000],
            [0x3fe8000000000000, 0x3ff6000000000000, 0x3ff4000000000000],
        ],
        [
            [0x3fffacd9283a3c66, 0xbfd298f138ef3d54, 0x3fe25acb84d76a1e],
            [0x400148ccf907300b, 0x3fe0a151ce772f08, 0x3fd305daa2b12e8c],
            [0x3ffd65f05c0cf94b, 0x3fe488300f2826ab, 0x3fdc5e70ba217554],
            [0x3ffa812f9238d59b, 0xbfc596696f1a9c1e, 0x3fe70716908f8d82],
            [0x40004d07d1f9835b, 0xbfd110601e504d56, 0x3fe5d0c7a2ef4556],
            [0x4001bf6836e39533, 0x3fe1659a5bc6a708, 0x3fd9f1d2dee0e4fb],
            [0x3ffe5326d7c5c39a, 0x3fe54c789c779eaa, 0x3fe1a5347b2895e2],
            [0x3ffb6e660df19fea, 0xbfc2854739dcbc20, 0x3fea7d12aea768ba],
        ],
        [
            [0x3fe8000000000000, 0x3ff0000000000000, 0x3fe4000000000000],
            [0x3ff4000000000000, 0x3ff0000000000000, 0x3fe4000000000000],
            [0x3ff4000000000000, 0x3ffc000000000000, 0x3fe4000000000000],
            [0x3fe8000000000000, 0x3ffc000000000000, 0x3fe4000000000000],
            [0x3fe8000000000000, 0x3ff0000000000000, 0x3ffa000000000000],
            [0x3ff4000000000000, 0x3ff0000000000000, 0x3ffa000000000000],
            [0x3ff4000000000000, 0x3ffc000000000000, 0x3ffa000000000000],
            [0x3fe8000000000000, 0x3ffc000000000000, 0x3ffa000000000000],
        ],
        [
            [0x3fd8000000000000, 0x3ff4000000000000, 0x3fe4000000000000],
            [0x3fe8000000000000, 0x3ff4000000000000, 0x3fe4000000000000],
            [0x3fe8000000000000, 0x3ff8000000000000, 0x3fe4000000000000],
            [0x3fd8000000000000, 0x3ff8000000000000, 0x3fe4000000000000],
            [0x3fd8000000000000, 0x3ff4000000000000, 0x3fe8000000000000],
            [0x3fe8000000000000, 0x3ff4000000000000, 0x3fe8000000000000],
            [0x3fe8000000000000, 0x3ff8000000000000, 0x3fe8000000000000],
            [0x3fd8000000000000, 0x3ff8000000000000, 0x3fe8000000000000],
        ],
        [
            [0x3fe4000000000000, 0x3fd0000000000000, 0x3ff0000000000000],
            [0x3ff2000000000000, 0x3fd0000000000000, 0x3ff0000000000000],
            [0x3ff2000000000000, 0x3ff4000000000000, 0x3ff0000000000000],
            [0x3fe4000000000000, 0x3ff4000000000000, 0x3ff0000000000000],
            [0x3fe4000000000000, 0x3fd0000000000000, 0x3ff8000000000000],
            [0x3ff2000000000000, 0x3fd0000000000000, 0x3ff8000000000000],
            [0x3ff2000000000000, 0x3ff4000000000000, 0x3ff8000000000000],
            [0x3fe4000000000000, 0x3ff4000000000000, 0x3ff8000000000000],
        ],
    ],
};

fn vec3(bits: [u64; 3]) -> Vec3 {
    Vec3::new(
        f64::from_bits(bits[0]),
        f64::from_bits(bits[1]),
        f64::from_bits(bits[2]),
    )
}

impl SnapScene {
    fn case(&self, name: &'static str) -> SnapCase {
        // Every box has the topology of `box_mesh`.
        let indices = box_mesh([0.0; 3], [1.0; 3]).indices;
        let meshes = self
            .boxes
            .iter()
            .enumerate()
            .map(|(i, corners)| {
                let positions = corners
                    .iter()
                    .map(|&c| {
                        let v = vec3(c);
                        Point3::new(v.x, v.y, v.z)
                    })
                    .collect();
                FaceFixture::new(
                    i as u64 + 1,
                    MeshRole::Solid,
                    TriMesh::new(positions, indices.clone()),
                )
            })
            .collect();
        let origin = vec3(self.origin);
        let frame = Frame3 {
            origin: Point3::new(origin.x, origin.y, origin.z),
            x: vec3(self.axes[0]),
            y: vec3(self.axes[1]),
            z: vec3(self.axes[2]),
        };
        SnapCase {
            name,
            meshes,
            view: OrthographicView::new(frame).with_near(f64::from_bits(self.near)),
            limits: VisibleFaceLimits::new(16, 100_000, 100_000, 10_000_000, 100_000, 1_000_000),
        }
    }
}

/// The scenes the conformance suite requires a provider to answer.
pub(crate) fn snap_cases() -> Vec<SnapCase> {
    vec![
        SEED_1_OBLIQUE.case("snap_seed_1_oblique"),
        SEED_9_ISOMETRIC.case("snap_seed_9_isometric"),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_frozen_scenes_are_closed_boxes() {
        for case in snap_cases() {
            assert!(!case.meshes.is_empty());
            for fixture in &case.meshes {
                assert_eq!(fixture.mesh.positions.len(), 8, "{}", case.name);
                assert_eq!(fixture.mesh.triangle_count(), 12, "{}", case.name);
            }
        }
    }
}
