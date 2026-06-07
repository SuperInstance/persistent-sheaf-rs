# Integration Guide: persistent-sheaf

## What This Crate Provides

- **`SimplicialComplex`** — Simplicial complex with Vietoris-Rips construction, Betti numbers, boundary operators
- **`Filtration`** — Nested sequence of simplicial complexes indexed by scale; built from point clouds
- **`PersistenceDiagram`** — Tracks how topological features (H₀, H₁) appear/disappear across filtration scales
- **`PersistencePair`** — A single (birth, death) pair in a persistence diagram
- **`CellularSheaf`** — Assigns data (stalks) to cells with linear restriction maps; supports constant, custom, and product sheaves
- **`SheafLaplacian`** — Generalizes graph Laplacian to encode geometric + non-geometric information; computes eigenvalues
- **`HodgeLaplacian`** — Hodge decomposition: harmonic 0-sections and 1-sections
- **`PushforwardSheaf`** / **`SimplicialMap`** — Push a sheaf forward along a simplicial map
- **`BottleneckDistance`** — Distance between persistence diagrams

This crate provides topological data analysis via persistent sheaf cohomology — combining persistent homology with sheaf theory for multi-modal data fusion. It discovers the shape of data, not just its statistics.

## How to Add This Crate

```bash
cargo add persistent-sheaf
```

```rust
use persistent_sheaf::{SimplicialComplex, Filtration, CellularSheaf, SheafLaplacian};

let points = vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.5, 0.866]];
let h0 = SimplicialComplex::h0_betti(&points, 0.5);
println!("Connected components: {h0}");

let complex = SimplicialComplex::vietoris_rips(&points, 1.0);
let sheaf = CellularSheaf::constant(complex, 2);
let lap = SheafLaplacian::from_sheaf(&sheaf);
println!("Eigenvalues: {:?}", lap.eigenvalues());
```

## Integration Points

### room-topology

- **Why**: room-topology defines the spatial structure of agent rooms/environments; persistent-sheaf discovers the topological structure of point cloud data. Together they provide a complete topological framework: room-topology for known structure, persistent-sheaf for discovered structure.
- **How**: Use room-topology's spatial graph as a `SimplicialComplex`, then construct a `CellularSheaf` that assigns agent state data to rooms with restriction maps along corridors.

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};

// Room layout as a simplicial complex
// Vertices = rooms, edges = corridors, triangles = shared spaces
let rooms = vec![vec![0.0, 0.0], vec![5.0, 0.0], vec![2.5, 4.0]];
let complex = SimplicialComplex::vietoris_rips(&rooms, 6.0);

// Assign agent data to rooms via sheaf
let sheaf = CellularSheaf::constant(complex, 3); // 3D agent state per room
let lap = SheafLaplacian::from_sheaf(&sheaf);
println!("Topology-aware eigenvalues: {:?}", lap.eigenvalues());
```

### wasserstein-agents

- **Why**: wasserstein-agents moves distributions optimally; persistent-sheaf reveals the topological structure that constrains WHERE agents can move. Transport plans should respect topology — agents can't move through walls.
- **How**: Build a `Filtration` from agent positions, compute persistence to identify connected components, then constrain wasserstein transport to stay within topologically connected regions.

```rust
use persistent_sheaf::{SimplicialComplex, Filtration, PersistenceDiagram};
use wasserstein_agents::AgentDistribution;

// Discover topology from agent positions
let agents = AgentDistribution::uniform(vec![
    vec![0.0, 0.0], vec![0.1, 0.1], vec![5.0, 5.0],
]);

let positions = agents.positions.clone();
let filt = Filtration::from_point_cloud(&positions, 20);
let diagram = filt.compute_persistence();
println!("Total persistence: {}", diagram.total_persistence(1.0));

// Connected components at scale ε define transport neighborhoods
let h0 = SimplicialComplex::h0_betti(&positions, 1.0);
println!("Transport neighborhoods: {h0}");
```

## For AI Agents

- **Context needed**: Point cloud data (agent positions, sensor readings, etc.), scale parameter ε, desired stalk dimension
- **Key imports**: `persistent_sheaf::{SimplicialComplex, Filtration, PersistenceDiagram, CellularSheaf, SheafLaplacian}`
- **Integration pattern**: Build `SimplicialComplex` from data → construct `CellularSheaf` with stalks → compute `SheafLaplacian` → extract eigenvalues for coordination analysis
- **Error handling**: Empty point clouds return empty complexes. Betti numbers may be approximate for large datasets. Filtration resolution (number of steps) controls accuracy vs. performance.

## For Humans

- **Prerequisites**: Basic topology (simplicial complexes, homology), understanding of persistence diagrams, sheaf theory basics
- **Learning path**: Start with `simplicial.rs` (complexes), then `filtration.rs` + `persistence.rs` (persistent homology), then `sheaf.rs` (sheaves on complexes), then `sheaf_laplacian.rs` (spectral analysis)
- **Common pitfalls**:
  - Vietoris-Rips scale ε too large → everything is connected (one component); too small → nothing is connected (isolated points)
  - Filtration resolution (number of steps) too low → misses topological features; too high → slow computation
  - Constant sheaves (same stalk everywhere) are simplest but least informative — use custom restriction maps for domain-specific structure
  - `total_persistence` is sensitive to the power parameter p — use p=1 for robustness, p=2 for sensitivity
