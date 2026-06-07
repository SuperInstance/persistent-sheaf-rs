//! # Persistent Sheaf
//!
//! Topological data analysis via persistent sheaf cohomology.
//! Combines persistent homology with sheaf theory for multi-modal data fusion.
//!
//! # Key Concepts
//!
//! - **Simplicial complex**: A collection of simplices closed under taking faces.
//!   Supports Vietoris-Rips construction from point clouds and distance matrices.
//! - **Filtration**: A nested sequence of simplicial complexes indexed by scale.
//! - **Cellular sheaf**: Assigns data (stalks) to cells of a complex with linear
//!   restriction maps between faces.
//! - **Sheaf Laplacian**: Generalizes the graph Laplacian to encode both geometric
//!   and non-geometric information.
//! - **Persistence diagram**: Tracks how topological features (H₀, H₁) appear and
//!   disappear across filtration scales.
//! - **Bottleneck distance**: Distance between persistence diagrams.
//!
//! # Quick Start
//!
//! ```rust
//! use persistent_sheaf::{
//!     SimplicialComplex, Filtration, PersistenceDiagram, CellularSheaf, SheafLaplacian,
//! };
//!
//! // Point cloud in 2D
//! let points = vec![
//!     vec![0.0, 0.0],
//!     vec![1.0, 0.0],
//!     vec![0.5, 0.866],
//! ];
//!
//! // H₀: connected components at a given threshold
//! let h0 = SimplicialComplex::h0_betti(&points, 0.5);
//! println!("Connected components at ε=0.5: {h0}");
//!
//! // Vietoris-Rips complex
//! let complex = SimplicialComplex::vietoris_rips(&points, 1.0);
//!
//! // Constant sheaf on the complex
//! let sheaf = CellularSheaf::constant(complex, 2);
//!
//! // Sheaf Laplacian
//! let lap = SheafLaplacian::from_sheaf(&sheaf);
//! let eigenvalues = lap.eigenvalues();
//! println!("Sheaf Laplacian eigenvalues: {:?}", eigenvalues);
//!
//! // Persistent homology across a filtration
//! let filt = Filtration::from_point_cloud(&points, 20);
//! let diagram = filt.compute_persistence();
//! let tp = diagram.total_persistence(1.0);
//! println!("Total persistence: {tp}");
//! ```

mod filtration;
mod laplacian;
mod persistence;
mod pushforward;
mod sheaf;
mod sheaf_laplacian;
mod simplicial;

pub use filtration::Filtration;
pub use laplacian::SheafLaplacian;
pub use persistence::{BottleneckDistance, PersistenceDiagram, PersistencePair};
pub use pushforward::{PushforwardSheaf, SimplicialMap};
pub use sheaf::CellularSheaf;
pub use sheaf_laplacian::{HodgeLaplacian, harmonic_0_sections, harmonic_1_sections};
pub use simplicial::SimplicialComplex;
