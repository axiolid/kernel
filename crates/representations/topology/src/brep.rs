//! Append-only typed arenas for B-rep topology.

use crate::{
    Edge, EdgeId, Face, FaceId, Loop, LoopId, Shell, ShellId, Solid, SolidId, Vertex, VertexId,
};

/// Owned B-rep. Generic geometry handles avoid a dependency cycle with the
/// model graph that stores exact curves and surfaces.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(serde::Serialize, serde::Deserialize),
    serde(deny_unknown_fields)
)]
pub struct BRep<Curve3, Curve2 = Curve3, Surface = Curve3> {
    vertices: Vec<Vertex>,
    edges: Vec<Edge<Curve3>>,
    loops: Vec<Loop<Curve2>>,
    faces: Vec<Face<Surface>>,
    shells: Vec<Shell>,
    solids: Vec<Solid>,
}

impl<Curve3, Curve2, Surface> Default for BRep<Curve3, Curve2, Surface> {
    fn default() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            loops: Vec::new(),
            faces: Vec::new(),
            shells: Vec::new(),
            solids: Vec::new(),
        }
    }
}

impl<Curve3, Curve2, Surface> BRep<Curve3, Curve2, Surface> {
    /// Fallibly reserve append-only arena capacity for bounded construction.
    pub fn try_reserve(
        &mut self,
        vertices: usize,
        edges: usize,
        loops: usize,
        faces: usize,
        shells: usize,
        solids: usize,
    ) -> Result<(), std::collections::TryReserveError> {
        self.vertices.try_reserve_exact(vertices)?;
        self.edges.try_reserve_exact(edges)?;
        self.loops.try_reserve_exact(loops)?;
        self.faces.try_reserve_exact(faces)?;
        self.shells.try_reserve_exact(shells)?;
        self.solids.try_reserve_exact(solids)?;
        Ok(())
    }

    /// Add a vertex and return its typed handle.
    pub fn add_vertex(&mut self, value: Vertex) -> VertexId {
        let id = VertexId::from_index(self.vertices.len());
        self.vertices.push(value);
        id
    }

    /// Add an edge and return its typed handle.
    pub fn add_edge(&mut self, value: Edge<Curve3>) -> EdgeId {
        let id = EdgeId::from_index(self.edges.len());
        self.edges.push(value);
        id
    }

    /// Add a loop and return its typed handle.
    pub fn add_loop(&mut self, value: Loop<Curve2>) -> LoopId {
        let id = LoopId::from_index(self.loops.len());
        self.loops.push(value);
        id
    }

    /// Add a face and return its typed handle.
    pub fn add_face(&mut self, value: Face<Surface>) -> FaceId {
        let id = FaceId::from_index(self.faces.len());
        self.faces.push(value);
        id
    }

    /// Add a shell and return its typed handle.
    pub fn add_shell(&mut self, value: Shell) -> ShellId {
        let id = ShellId::from_index(self.shells.len());
        self.shells.push(value);
        id
    }

    /// Add a solid and return its typed handle.
    pub fn add_solid(&mut self, value: Solid) -> SolidId {
        let id = SolidId::from_index(self.solids.len());
        self.solids.push(value);
        id
    }

    /// Vertices in stable insertion order.
    pub fn vertices(&self) -> &[Vertex] {
        &self.vertices
    }

    /// Edge handle at a dense arena index, if present.
    pub fn edge_id_at(&self, index: usize) -> Option<EdgeId> {
        (index < self.edges.len()).then(|| EdgeId::from_index(index))
    }

    /// Loop handle at a dense arena index, if present.
    pub fn loop_id_at(&self, index: usize) -> Option<LoopId> {
        (index < self.loops.len()).then(|| LoopId::from_index(index))
    }

    /// All edges in insertion order.
    pub fn edges(&self) -> &[Edge<Curve3>] {
        &self.edges
    }

    /// Loops in stable insertion order.
    pub fn loops(&self) -> &[Loop<Curve2>] {
        &self.loops
    }

    /// Face handle at a dense arena index, if present.
    pub fn face_id_at(&self, index: usize) -> Option<FaceId> {
        (index < self.faces.len()).then(|| FaceId::from_index(index))
    }

    /// All faces in insertion order.
    pub fn faces(&self) -> &[Face<Surface>] {
        &self.faces
    }

    /// Shells in stable insertion order.
    pub fn shells(&self) -> &[Shell] {
        &self.shells
    }

    /// Solids in stable insertion order.
    pub fn solids(&self) -> &[Solid] {
        &self.solids
    }
}
