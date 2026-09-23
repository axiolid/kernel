//! Capability ledger: what Axiolid can do, measured against OCCT and CGAL.
//!
//! `architecture/capability-ledger.toml` is the durable record of a one-off
//! like-for-like audit. It exists so the next contributor does not rerun
//! that audit: `cargo xtask gaps` answers "what is missing, where does the
//! reference implementation live, which issue tracks it" in one command.
//!
//! `check` is the half that keeps the answer true. A ledger that nothing
//! verifies rots into a wish list, so the gate enforces what can be checked
//! mechanically: every evidence path exists, every cited symbol is defined
//! in its file, every issue key resolves, and a row that claims capability
//! cites evidence. Whether the prose is right stays a review question.

mod model;
mod query;
mod verify;

pub use query::{list, next, show};
pub use verify::check;
