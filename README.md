# persistent-sheaf

**Persistent sheaf cohomology and cellular sheaf Laplacians in pure Rust — multi-modal data fusion through topology.**

Combines persistent homology with sheaf theory. A cellular sheaf assigns data spaces (stalks) to cells of a simplicial complex with linear restriction maps between them. The sheaf Laplacian generalizes the graph Laplacian by encoding both geometric and non-geometric information. This library builds filtrations, computes sheaf cohomology, and constructs the sheaf Laplacian for spectral analysis.

## What This Gives You

- **Simplicial complexes** — vertices, edges, triangles, tetrahedra with Vietoris-Rips construction
- **Cellular sheaves** — constant sheaves, weighted sheaves, custom restriction maps
- **Sheaf Laplacian** — L_F generalizes graph Laplacian with stalk/restriction information
- **Persistent homology** — birth-death pairs across filtration scales
- **Persistence diagrams** — total persistence, essential features, bottleneck distance
- **Zero dependencies** — pure Rust

## Quick Start

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian, Filtration};

// Build distance matrix
let distances = vec![/* N×N */];

// Vietoris-Rips filtration
let filt = Filtration::from_distance_matrix(&distances, 20);

// Constant sheaf on the complex
let complex = SimplicialComplex::vietoris_rips(&distances, 0.5);
let sheaf = CellularSheaf::constant(complex, 3); // 3D stalks

// Sheaf Laplacian
let lap = SheafLaplacian::from_sheaf(&sheaf);
```

## API Reference

| Module | Key Types |
|--------|-----------|
| `simplicial` | `SimplicialComplex` |
| `sheaf` | `CellularSheaf` — constant, weighted, custom |
| `laplacian` | `SheafLaplacian` |
| `filtration` | `Filtration` |
| `persistence` | `PersistenceDiagram`, `PersistencePair` |

## Testing

```bash
cargo test
```

## Installation

```toml
[dependencies]
persistent-sheaf = { git = "https://github.com/SuperInstance/persistent-sheaf" }
```

## How It Fits

Part of the SuperInstance ecosystem:

- **[persistent-social](https://github.com/SuperInstance/persistent-social)** — Social network TDA in Go
- **[gpu-sheaf-laplacian](https://github.com/SuperInstance/gpu-sheaf-laplacian)** — CUDA sheaf Laplacian
- **persistent-sheaf** — Rust sheaf cohomology library (this repo)

## License

MIT
