//! NEC-2 antenna modeling engine: the method of moments for thin wires.
//!
//! This crate is a Rust port of the public-domain Numerical Electromagnetics
//! Code, version 2 (NEC-2). It solves the electric-field integral equation for
//! the currents on a structure of thin wire segments, and the magnetic-field
//! integral equation for the currents on closed surfaces made of patches. From
//! the currents it computes input impedances, radiation patterns and near
//! fields, in free space or over a perfect or lossy ground.
//!
//! The crate only computes: it reads no files and starts no threads. Decks of
//! NEC cards are handled by the `nec2-deck` crate.
//!
//! ```
//! use nec2::{Engine, Excitation, Geometry};
//! use num_complex::Complex64;
//!
//! // Half-wave dipole of 21 segments, fed at the centre.
//! let mut geometry = Geometry::new();
//! geometry.wire(1, 21, [0.0, 0.0, -0.25], [0.0, 0.0, 0.25], 0.001, 1.0, 1.0)?;
//! let structure = geometry.connect(0)?;
//!
//! let mut engine = Engine::new(&structure);
//! engine.set_frequency(299.8);
//! engine.fill_and_factor(1.0, false)?;
//! engine.set_sources(&[(11, Complex64::new(1.0, 0.0))], &[]);
//! engine.excite(&Excitation::Sources)?;
//! let report = engine.solve(&[], false, false)?;
//!
//! let z = report.sources[0].impedance;
//! assert!((z.re - 80.0).abs() < 5.0 && z.im.abs() < 60.0);
//! # Ok::<(), nec2::Error>(())
//! ```

// The numerical routines keep the constants of NEC-2 as published (rounded
// values of pi and of other quantities) and its index-based loops, so that
// they can be compared line by line with the original.
#![allow(clippy::approx_constant, clippy::excessive_precision, clippy::needless_range_loop)]

pub mod consts;
mod engine;
mod error;
mod geometry;

pub use engine::{
    Asymmetry, AverageGain, Engine, Excitation, FarFieldRow, GainSample, GroundInfo, GroundParams, GroundType,
    GroundWaveRow, Load, LoadLine, LoadingReport, Network, NetworkKind, NormalizedGain, PatternRequest, PortRecord,
    RadiationPattern, SolveReport,
};
pub use error::{Error, Result};
pub use geometry::{Geometry, HelixInfo, Overlap, OverlapCheck, Patch, Structure, isegno};
