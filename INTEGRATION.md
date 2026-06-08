# Integration Guide: persistent-sheaf

## What This Crate Provides

- **`SimplicialComplex`** — Vietoris-Rips construction from point clouds / distance matrices, flag complexes, Betti numbers (β₀, β₁), Euler characteristic, connected components
- **`Filtration`** — Nested sequence of complexes indexed by scale; persistent homology computation
- **`PersistenceDiagram`** / **`PersistencePair`** — Birth-death tracking, bottleneck distance, Wasserstein distance, Betti curves, pruning
- **`CellularSheaf`** — Data stalks with linear restriction maps; constant and weighted sheaves; cohomology dimensions
- **`SheafLaplacian`** — Eigenvalue decomposition, Fiedler value, power iteration, edge/triangle local Laplacians
- **`HodgeLaplacian`** — Coboundary operators, harmonic 0/1-sections, sheaf Betti numbers, spectral gap
- **`SimplicialMap`** / **`PushforwardSheaf`** — Direct image sheaves along simplicial maps
- **`RoomTopologySheaf`** — Floor-plan sheaf with per-room stalks and door restriction maps
- **`PersistenceReport`** — Connected-component analysis, isolated rooms, loop detection

This crate combines persistent homology with sheaf theory for multi-modal topological data analysis: tracking how shape and data coherence evolve across scales.

## How to Add This Crate

```bash
cargo add persistent-sheaf
```

```rust
use persistent_sheaf::{
    SimplicialComplex, Filtration, PersistenceDiagram,
    CellularSheaf, SheafLaplacian,
};
```

## Cross-Repo Connections

### With `wasserstein-agents`: Optimal Transport Between Persistence Diagrams

Compare the topological signatures of two agent formations using bottleneck distance, then compute the optimal transport plan:

```rust
use persistent_sheaf::{Filtration, PersistenceDiagram};
use wasserstein_agents::transport::SinkhornSolver;

fn formation_distance(formation_a: &[Vec<f64>], formation_b: &[Vec<f64>]) -> f64 {
    let filt_a = Filtration::from_point_cloud(formation_a, 20);
    let filt_b = Filtration::from_point_cloud(formation_b, 20);
    let diag_a = filt_a.compute_persistence();
    let diag_b = filt_b.compute_persistence();
    diag_a.bottleneck_distance(&diag_b)
}
```

### With `spectral-fleet`: Spectral Clustering on Sheaf Laplacians

Use the sheaf Laplacian eigenvalues to rank agents by their topological influence:

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};
use spectral_fleet::spectral_clustering::spectral_clustering;
use rand::thread_rng;

fn sheaf_based_clustering(points: &[Vec<f64>], k: usize) -> Vec<usize> {
    let complex = SimplicialComplex::vietoris_rips(points, 1.5);
    let sheaf = CellularSheaf::constant(complex, 1);
    let lap = SheafLaplacian::from_sheaf(&sheaf);
    let eigenvalues = lap.eigenvalues();
    println!("Fiedler value (connectivity): {:.4}", lap.fiedler_value());

    // Build affinity from eigenvalue gaps
    let mut affinity: Vec<Vec<f64>> = vec![vec![0.0; points.len()]; points.len()];
    for i in 0..points.len() {
        for j in 0..points.len() {
            affinity[i][j] = (-eigenvalues.get(i.min(j)).copied().unwrap_or(0.0)).exp();
        }
    }
    let mut rng = thread_rng();
    spectral_clustering(&affinity, k, 50, 1e-6, &mut rng).unwrap().labels
}
```

### With `hodge-consensus`: Hodge Decomposition of Agent Disagreements

Model a fleet as a simplicial complex and use sheaf cohomology to detect irreducible disagreement cycles:

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, HodgeLaplacian, sheaf_betti_numbers};
use hodge_consensus::{disagreement_matrix, decompose_disagreement};

fn fleet_disagreement_topology(agents: &[usize], edges: &[(usize, usize)]) -> (usize, usize) {
    let mut complex = SimplicialComplex::new();
    for &(a, b) in edges {
        complex.add_edge(a, b);
    }
    let sheaf = CellularSheaf::constant(complex, 1);
    let delta = persistent_sheaf::sheaf_laplacian::coboundary_matrix(
        agents.len(), edges, 1, &sheaf.restriction_maps,
    );
    let hodge = HodgeLaplacian::from_coboundary(delta);
    sheaf_betti_numbers(&hodge, 1e-8)
}
```

## Design Patterns

### Pattern: Persistent Topology Health Check

Monitor fleet cohesion by tracking how β₀ (connected components) changes over time:

```rust
use persistent_sheaf::SimplicialComplex;

fn topology_health(positions: &[Vec<f64>], threshold: f64) -> bool {
    let beta0 = SimplicialComplex::h0_betti(positions, threshold);
    beta0 == 1 // fleet should be a single connected component
}
```

### Pattern: Room Topology Sheaf for Building Navigation

Model a building floor plan as a sheaf and query room-to-room restrictions:

```rust
use persistent_sheaf::spatial_analysis::RoomTopologySheaf;

fn building_sheaf(adjacency: &[Vec<usize>]) -> RoomTopologySheaf {
    let mut sheaf = RoomTopologySheaf::from_adjacency(adjacency, 3);
    sheaf.set_stalk(0, vec![22.5, 45.0, 60.0]);
    sheaf.set_stalk(1, vec![21.0, 50.0, 55.0]);
    sheaf
}
```

### With `ga-core`: Conformal Point Cloud Embedding Before Topology

Embed agent positions into conformal space for rotation-invariant persistent homology:

```rust
use persistent_sheaf::SimplicialComplex;
use ga_core::conformal::Conformal;

fn conformal_vr_complex(positions: &[[f64; 3]], epsilon: f64) -> SimplicialComplex {
    let embedded: Vec<Vec<f64>> = positions.iter()
        .map(|p| {
            let c = Conformal::embed_point(*p);
            vec![c[0], c[1], c[2], c[3], c[4]]
        })
        .collect();
    SimplicialComplex::vietoris_rips(&embedded, epsilon)
}
```

### With `categorical-agents`: Functorial Sheaf Pushforward

Treat simplicial maps as categorical morphisms and pushforward sheaves as functor images:

```rust
use persistent_sheaf::{CellularSheaf, SimplicialComplex, SimplicialMap, PushforwardSheaf};
use categorical_agents::monad::ListMonad;

fn pushforward_as_functor(sheaf: &CellularSheaf, map: &SimplicialMap, target: SimplicialComplex) -> PushforwardSheaf {
    PushforwardSheaf::compute(sheaf, map, target)
}
```

### With `agent-homeostasis`: Sheaf Connectivity as Homeostatic Sensor

Use the Fiedler value (algebraic connectivity) of the sheaf Laplacian as a sensor input for PID regulation:

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};
use agent_homeostasis::sensor::SensorReading;

fn connectivity_sensor(positions: &[Vec<f64>], threshold: f64) -> SensorReading {
    let complex = SimplicialComplex::vietoris_rips(positions, threshold);
    let sheaf = CellularSheaf::constant(complex, 1);
    let lap = SheafLaplacian::from_sheaf(&sheaf);
    let fiedler = lap.fiedler_value();
    SensorReading::new("connectivity", fiedler)
}
```

## Design Patterns

### Pattern: Multi-Scale Topology Pipeline

Run persistent homology at multiple resolutions to detect both fine and coarse structure:

```rust
use persistent_sheaf::{SimplicialComplex, Filtration, PersistenceDiagram};

fn multiscale_topology(points: &[Vec<f64>]) -> Vec<PersistenceDiagram> {
    let mut diagrams = vec![];
    for steps in [10, 20, 50] {
        let filt = Filtration::from_point_cloud(points, steps);
        diagrams.push(filt.compute_persistence());
    }
    diagrams
}
```

### Pattern: Sheaf Filtration for Time-Varying Data

Build a sheaf on each complex in a filtration and track how global sections evolve:

```rust
use persistent_sheaf::{Filtration, CellularSheaf};

fn track_global_sections(filt: &Filtration, stalk_dim: usize) -> Vec<usize> {
    let sheaves = filt.make_sheaf_filtration(stalk_dim);
    sheaves.iter()
        .map(|s| s.global_section_dimension())
        .collect()
}
```

### Pattern: Pushforward for Fleet Aggregation

Collapse a detailed agent graph onto a coarse command structure using simplicial maps:

```rust
use persistent_sheaf::{CellularSheaf, SimplicialComplex, SimplicialMap, PushforwardSheaf};

fn aggregate_fleet(detail_sheaf: &CellularSheaf, squad_map: &[usize], target: SimplicialComplex) -> PushforwardSheaf {
    let map = SimplicialMap::new(squad_map.to_vec());
    PushforwardSheaf::compute(detail_sheaf, &map, target)
}
```
