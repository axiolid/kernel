//! Construction and validation of point-sampled values.

use axiolid_core::{Point3, Scalar, Vec3};
use axiolid_pointcloud::{Colour, PointCloud, PointCloudError};

fn sample_points() -> Vec<Point3> {
    vec![
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ]
}

/// A cloud of finite points is accepted and keeps the caller's order.
///
/// Order matters because every attribute channel is indexed by it, so it
/// is asserted rather than assumed.
#[test]
fn finite_points_are_accepted_in_the_given_order() {
    let cloud = PointCloud::new(sample_points()).expect("finite points construct");

    assert_eq!(cloud.len(), 3);
    assert!(!cloud.is_empty());
    assert_eq!(cloud.points()[1], Point3::new(1.0, 0.0, 0.0));
    assert!(!cloud.has_normals());
}

/// An empty cloud is a value, not an error.
///
/// A filter that removes everything must be able to say so without
/// inventing a failure.
#[test]
fn an_empty_cloud_is_legitimate() {
    let cloud = PointCloud::empty();
    assert!(cloud.is_empty());
    assert_eq!(cloud.len(), 0);
    assert_eq!(cloud.bounds(), None);

    let also_empty = PointCloud::new(Vec::new()).expect("an empty vector constructs");
    assert!(also_empty.is_empty());
}

/// A non-finite coordinate is refused, naming the point.
#[test]
fn a_non_finite_point_is_refused_by_index() {
    let mut points = sample_points();
    points.push(Point3::new(0.0, Scalar::NAN, 0.0));

    let error = PointCloud::new(points).expect_err("NaN must be refused");
    assert!(
        matches!(error, PointCloudError::NonFinitePoint { index: 3, .. }),
        "the offending index must be reported, got {error:?}"
    );

    let mut infinite = sample_points();
    infinite[0] = Point3::new(Scalar::INFINITY, 0.0, 0.0);
    assert!(matches!(
        PointCloud::new(infinite),
        Err(PointCloudError::NonFinitePoint { index: 0, .. })
    ));
}

/// A channel of the wrong length is refused at construction.
///
/// Otherwise the mismatch surfaces later as an out-of-bounds index inside
/// whichever algorithm happened to touch it first.
#[test]
fn a_channel_of_the_wrong_length_is_refused() {
    let cloud = PointCloud::new(sample_points()).expect("constructs");

    let error = cloud
        .clone()
        .with_normals(vec![Vec3::new(0.0, 0.0, 1.0)])
        .expect_err("one normal for three points must be refused");
    assert!(
        matches!(
            error,
            PointCloudError::ChannelLengthMismatch {
                channel: "normals",
                found: 1,
                expected: 3
            }
        ),
        "the mismatch must name the channel and both counts, got {error:?}"
    );

    assert!(matches!(
        cloud.clone().with_colours(vec![Colour::new(1, 2, 3)]),
        Err(PointCloudError::ChannelLengthMismatch {
            channel: "colours",
            ..
        })
    ));
    assert!(matches!(
        cloud.with_intensities(vec![0.5]),
        Err(PointCloudError::ChannelLengthMismatch {
            channel: "intensities",
            ..
        })
    ));
}

/// A zero-length normal is refused rather than normalised.
///
/// Normalising it would pick an arbitrary direction that is
/// indistinguishable downstream from a measured one.
#[test]
fn a_degenerate_normal_is_refused_rather_than_invented() {
    let cloud = PointCloud::new(sample_points()).expect("constructs");

    let error = cloud
        .clone()
        .with_normals(vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::ZERO,
            Vec3::new(1.0, 0.0, 0.0),
        ])
        .expect_err("a zero normal must be refused");
    assert!(
        matches!(error, PointCloudError::DegenerateNormal { index: 1, .. }),
        "the offending normal must be named, got {error:?}"
    );

    assert!(matches!(
        cloud.with_normals(vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(Scalar::NAN, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        ]),
        Err(PointCloudError::DegenerateNormal { index: 1, .. })
    ));
}

/// Valid channels attach and read back.
#[test]
fn attributes_attach_and_read_back() {
    let cloud = PointCloud::new(sample_points())
        .expect("constructs")
        .with_normals(vec![
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
        ])
        .expect("normals attach")
        .with_colours(vec![
            Colour::new(255, 0, 0),
            Colour::new(0, 255, 0),
            Colour::new(0, 0, 255),
        ])
        .expect("colours attach")
        .with_intensities(vec![0.1, 0.5, 0.9])
        .expect("intensities attach");

    assert!(cloud.has_normals());
    assert_eq!(
        cloud.normals().expect("normals")[2],
        Vec3::new(1.0, 0.0, 0.0)
    );
    assert_eq!(cloud.colours().expect("colours")[0], Colour::new(255, 0, 0));
    assert_eq!(cloud.intensities().expect("intensities")[1], 0.5);
}

/// A cloud with no attributes reports none, rather than empty channels.
///
/// The distinction matters: "no normals were captured" and "normals were
/// captured and are all zero" are different claims about the sensor.
#[test]
fn absent_channels_are_absent_not_empty() {
    let cloud = PointCloud::new(sample_points()).expect("constructs");
    assert_eq!(cloud.normals(), None);
    assert_eq!(cloud.colours(), None);
    assert_eq!(cloud.intensities(), None);
}

/// An unnormalised normal is kept as given.
///
/// Length may carry confidence in some capture pipelines, so it is the
/// caller's to interpret; only a direction that does not exist is refused.
#[test]
fn normals_are_not_silently_renormalised() {
    let cloud = PointCloud::new(sample_points())
        .expect("constructs")
        .with_normals(vec![
            Vec3::new(0.0, 0.0, 5.0),
            Vec3::new(0.0, 2.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
        ])
        .expect("normals attach");

    assert_eq!(
        cloud.normals().expect("normals")[0],
        Vec3::new(0.0, 0.0, 5.0),
        "a normal must be stored as supplied"
    );
}

/// Bounds cover every point.
#[test]
fn bounds_enclose_every_point() {
    let cloud = PointCloud::new(vec![
        Point3::new(-1.0, 2.0, 0.5),
        Point3::new(3.0, -4.0, 0.0),
        Point3::new(0.0, 0.0, 7.0),
    ])
    .expect("constructs");

    let (min, max) = cloud.bounds().expect("a non-empty cloud has bounds");
    assert_eq!(min, Point3::new(-1.0, -4.0, 0.0));
    assert_eq!(max, Point3::new(3.0, 2.0, 7.0));
}

/// Duplicate points are preserved, not silently collapsed.
///
/// Repeated returns are real in scan data and their multiplicity carries
/// information; deduplicating here would destroy it before any consumer
/// could decide whether it mattered.
#[test]
fn duplicate_points_are_preserved() {
    let repeated = Point3::new(1.0, 1.0, 1.0);
    let cloud = PointCloud::new(vec![repeated, repeated, repeated]).expect("constructs");
    assert_eq!(cloud.len(), 3, "duplicates must not be collapsed");
}
