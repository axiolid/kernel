//! Complete graph compilation capability.

use crate::CompileOutcome;
use axiolid_mesh::TriMesh;
use axiolid_model::{GeometryGraph, NodeId};

use axiolid_contracts::{Backend, ExecutionOptions, GeomResult, OutputBound, ScratchRequirement};

/// Backend/orchestrator capable of lowering any advertised graph node to mesh.
///
/// Implementations must return [`axiolid_contracts::GeomError::Unsupported`] for a node they
/// do not support. They must not silently omit it or approximate exact geometry
/// unless the execution policy explicitly permits that precision.
pub trait MeshCompiler: Backend {
    /// Scratch this compiler needs beyond the graph and the produced meshes.
    ///
    /// Defaults to [`ScratchRequirement::Unbounded`]: an unaudited compiler is
    /// treated as unbudgetable rather than silently assumed cheap.
    fn scratch_requirement(&self) -> ScratchRequirement {
        ScratchRequirement::Unbounded
    }

    /// Outputs produced per requested root.
    ///
    /// Compilation is one mesh per root, so the destination size is known
    /// before the batch runs and workers can write disjoint slots.
    fn output_bound(&self) -> OutputBound {
        OutputBound::OneToOne
    }

    /// Compile one root.
    fn compile_mesh(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<TriMesh>;

    /// Compile one root and report what happened to each attribute channel.
    ///
    /// The default wraps [`Self::compile_mesh`] and reports
    /// [`CompileOutcome::attribute_fates`] as `None` (not tracked). A
    /// compiler that carries channels overrides this so a caller can tell a
    /// dropped channel from one that was never there.
    fn compile_mesh_reported(
        &self,
        graph: &GeometryGraph,
        root: NodeId,
        options: &ExecutionOptions,
    ) -> GeomResult<CompileOutcome> {
        self.compile_mesh(graph, root, options)
            .map(CompileOutcome::untracked)
    }

    /// Compile roots into a caller-provided buffer.
    ///
    /// This is the seam a batching implementation should override: the caller
    /// owns the destination, so a provider can reserve once from
    /// [`Self::output_bound`] and have workers write disjoint slots instead of
    /// growing a vector under a lock. `destination` is appended to, never
    /// cleared, so results can be accumulated across calls.
    ///
    /// The default is a serial loop over [`Self::compile_mesh`], which stays the
    /// only required primitive.
    fn compile_mesh_batch_into(
        &self,
        graph: &GeometryGraph,
        roots: &[NodeId],
        options: &ExecutionOptions,
        destination: &mut Vec<TriMesh>,
    ) -> GeomResult<()> {
        destination.reserve(roots.len());
        for &root in roots {
            destination.push(self.compile_mesh(graph, root, options)?);
        }
        Ok(())
    }

    /// Compile roots as a batch.
    ///
    /// Convenience over [`Self::compile_mesh_batch_into`]; overriding that instead
    /// gives both call shapes the batched behaviour.
    fn compile_mesh_batch(
        &self,
        graph: &GeometryGraph,
        roots: &[NodeId],
        options: &ExecutionOptions,
    ) -> GeomResult<Vec<TriMesh>> {
        let mut destination = Vec::with_capacity(roots.len());
        self.compile_mesh_batch_into(graph, roots, options, &mut destination)?;
        Ok(destination)
    }
}
