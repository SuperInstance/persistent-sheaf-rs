# persistent-sheaf-rs

Topological data analysis via persistent sheaf cohomology.

Combines persistent homology with sheaf theory for multi-modal data fusion. Point clouds become simplicial complexes. Simplicial complexes become filtrations. Filtrations produce persistence diagrams tracking how topological features appear and disappear across scales. Sheaves add structure — each cell carries data, each face map is a linear transformation, and the sheaf Laplacian reveals the geometry of consistency.

Part of the **sunset-ecosystem**: agent state distributions from `wasserstein-agents-rs` become point clouds, which become filtrations, which become persistence diagrams. `conservation-law` enforces that topological invariants are preserved under fleet reconfiguration. `si-fleet-api` uses Betti numbers to detect when the fleet's connectivity topology has changed.

## The Math

### Simplicial Complexes

A **simplicial complex** $K$ is a collection of simplices closed under taking faces. A $k$-simplex is the convex hull of $k+1$ vertices:
- 0-simplex: vertex
- 1-simplex: edge
- 2-simplex: triangle
- 3-simplex: tetrahedron

The **Euler characteristic** is:

$$\chi(K) = \sum_{k=0}^{n} (-1)^k |K_k| = V - E + F - T$$

### Vietoris-Rips Complex

Given a point cloud $X = \{x_1, \ldots, x_n\}$ and threshold $\varepsilon > 0$:

$$\text{VR}_\varepsilon(X) = \{\sigma \subseteq X : d(x_i, x_j) \leq \varepsilon \text{ for all } x_i, x_j \in \sigma\}$$

As $\varepsilon$ increases, the complex grows — this is a **filtration**.

### Betti Numbers

The $k$-th Betti number $\beta_k$ counts the number of $k$-dimensional "holes":
- $\beta_0$: connected components
- $\beta_1$: independent cycles (loops)
- $\beta_2$: enclosed cavities

By the Euler-Poincaré formula: $\chi = \sum_k (-1)^k \beta_k$.

### Persistent Homology

A filtration $K_0 \subseteq K_1 \subseteq \cdots \subseteq K_n$ induces homology groups $H_k(K_i)$. A **birth-death pair** $(b, d)$ records when a topological feature appears at scale $b$ and disappears at scale $d$. The collection of all pairs is the **persistence diagram**.

### Bottleneck Distance

The distance between two persistence diagrams $D_1$ and $D_2$:

$$d_B(D_1, D_2) = \inf_{\text{matchings } \phi} \sup_{p \in D_1} \|p - \phi(p)\|_\infty$$

### Cellular Sheaves

A **cellular sheaf** $F$ on a complex $K$ assigns:
- A vector space $F(\sigma)$ (the **stalk**) to each cell $\sigma$
- A linear map $F_{\sigma \to \tau}: F(\sigma) \to F(\tau)$ (**restriction map**) for each face $\tau < \sigma$

### Sheaf Laplacian

The **sheaf Laplacian** $L_F$ generalizes the graph Laplacian:

$$(L_F f)(v) = \sum_{v \sim w} F_{v \leftarrow w}^\top (F_{v \leftarrow w} f(v) - F_{w \leftarrow v} f(w))$$

When all stalks are $\mathbb{R}$ and all restriction maps are $1$, this reduces to the standard graph Laplacian $L = D - A$.

The **Fiedler value** (second-smallest eigenvalue of $L_F$) measures algebraic connectivity. The **harmonic sections** (kernel of $L_F$) are the global sections of the sheaf — the consistent assignments across all cells.

## Installation

```toml
[dependencies]
persistent-sheaf-rs = { git = "https://github.com/SuperInstance/persistent-sheaf-rs" }
```

## Usage

### Building Simplicial Complexes

```rust
use persistent_sheaf::SimplicialComplex;

let mut complex = SimplicialComplex::new();
complex.add_edge(0, 1);
complex.add_edge(1, 2);
complex.add_edge(2, 0); // triangle (edges only)

println!("Vertices: {}", complex.num_simplices(0)); // 3
println!("Edges: {}", complex.num_simplices(1));    // 3
println!("Triangles: {}", complex.num_simplices(2)); // 0

let betti = complex.betti_numbers();
println!("β₀ = {} (one component)", betti[0]);
println!("β₁ = {} (one loop)", betti[1]);

let euler = complex.euler_characteristic();
println!("χ = {} (3 - 3 + 0)", euler);

// Fill the triangle: β₁ drops to 0
let mut filled = SimplicialComplex::new();
filled.add_triangle(0, 1, 2);
let filled_betti = filled.betti_numbers();
println!("After filling: β₀={}, β₁={}", filled_betti[0], filled_betti[1]); // 1, 0
```

### Vietoris-Rips Complexes from Point Clouds

```rust
use persistent_sheaf::SimplicialComplex;

// Three points in 2D: an equilateral triangle
let points = vec![
    vec![0.0, 0.0],
    vec![1.0, 0.0],
    vec![0.5, 0.866],
];

// At ε = 0.5: no edges connected
let small = SimplicialComplex::vietoris_rips(&points, 0.5);
println!("ε=0.5: edges={}", small.num_simplices(1)); // 0

// At ε = 1.0: all edges connected (sides ≈ 1.0)
let medium = SimplicialComplex::vietoris_rips(&points, 1.0);
println!("ε=1.0: edges={}, triangles={}",
    medium.num_simplices(1), medium.num_simplices(2)); // 3, 1

// H₀: how many components at a given scale?
let h0_05 = SimplicialComplex::h0_betti(&points, 0.5); // 3 (isolated)
let h0_10 = SimplicialComplex::h0_betti(&points, 1.0); // 1 (connected)
println!("Components at ε=0.5: {}, at ε=1.0: {}", h0_05, h0_10);

// H₁: how many independent cycles?
let h1 = SimplicialComplex::h1_betti(&points, 1.0);
println!("Cycles at ε=1.0: {}", h1); // 0 (triangle is filled)
```

### Filtrations and Persistent Homology

```rust
use persistent_sheaf::{SimplicialComplex, Filtration, PersistenceDiagram};

// Point cloud: two clusters that merge at a certain scale
let points = vec![
    vec![0.0, 0.0],
    vec![0.1, 0.0],
    vec![2.0, 0.0],
    vec![2.1, 0.0],
];

// Build a filtration with 10 steps
let filtration = Filtration::from_point_cloud(&points, 10);
println!("Filtration: {} steps", filtration.len());

// Compute persistent homology
let diagram = filtration.compute_persistence();
println!("Persistence pairs: {}", diagram.len());

// Inspect H₀ features (component merges)
let h0 = diagram.filter_dimension(0);
println!("H₀ features: {}", h0.len());
for pair in h0 {
    if pair.is_essential() {
        println!("  Essential: born at {:.2}, survives to ∞", pair.birth);
    } else {
        println!("  Transient: born at {:.2}, dies at {:.2} (persistence={:.2})",
            pair.birth, pair.death, pair.persistence());
    }
}

// Total persistence: sum of persistence values
let tp = diagram.total_persistence(1.0);
println!("Total persistence (p=1): {:.4}", tp);

// Most persistent feature
if let Some(mp) = diagram.most_persistent() {
    println!("Most persistent: born {:.2}, died {:.2}, dim {}",
        mp.birth, mp.death, mp.dimension);
}

// Betti curve: how many features alive at each threshold?
let thresholds: Vec<f64> = (0..20).map(|i| i as f64 * 0.15).collect();
let curve = diagram.betti_curve(&thresholds);
println!("Betti curve: {:?}", curve);

// Prune short-lived features (noise)
let mut pruned = diagram.clone();
pruned.prune(0.3);
println!("After pruning: {} features remain", pruned.len());
```

### Bottleneck Distance Between Diagrams

```rust
use persistent_sheaf::PersistenceDiagram;

let mut d1 = PersistenceDiagram::new();
d1.add(0.0, 1.0, 0);
d1.add(0.5, 2.5, 0);
d1.add(0.0, 0.3, 1);

let mut d2 = PersistenceDiagram::new();
d2.add(0.1, 1.1, 0);
d2.add(0.6, 2.4, 0);
d2.add(0.0, 0.35, 1);

let bn = d1.bottleneck_distance(&d2);
println!("Bottleneck distance: {:.4}", bn);

// Identical diagrams → distance 0
let d3 = d1.clone();
let self_dist = d1.bottleneck_distance(&d3);
println!("Self-distance: {:.6}", self_dist); // ≈ 0
```

### Filtration from Distance Matrices

```rust
use persistent_sheaf::Filtration;

// Precomputed distance matrix (3 nodes in a line)
let distances = vec![
    vec![0.0, 1.0, 3.0],
    vec![1.0, 0.0, 1.0],
    vec![3.0, 1.0, 0.0],
];

let filtration = Filtration::from_distance_matrix(&distances, 10);

println!("Thresholds: {:?}", filtration.thresholds());
println!("Max threshold: {:.2}", filtration.max_threshold());

let diagram = filtration.compute_persistence();
println!("Features detected: {}", diagram.len());
println!("Essential features: {}", diagram.num_essential());

// At each step, check the complex
for i in 0..filtration.len() {
    if let Some((threshold, complex)) = filtration.get(i) {
        if i % 3 == 0 {
            println!("  ε={:.2}: V={}, E={}, T={}",
                threshold,
                complex.num_simplices(0),
                complex.num_simplices(1),
                complex.num_simplices(2));
        }
    }
}
```

### Cellular Sheaves

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf};

// Build a complex
let mut complex = SimplicialComplex::new();
complex.add_edge(0, 1);
complex.add_edge(1, 2);
complex.add_edge(2, 0);

// Constant sheaf: all stalks are R³, all restriction maps are identity
let sheaf = CellularSheaf::constant(complex.clone(), 3);
println!("Stalk dimension: {}", sheaf.stalk_dimension);
println!("Edge restriction maps: {}", sheaf.restriction_maps.len());

// Sheaf cohomology
let h0 = sheaf.cohomology_dimension(0);
println!("H⁰ dimension: {}", h0); // 3 (matches stalk dim for connected complex)

// Weighted sheaf: edges have scalar weights
let weighted = CellularSheaf::from_weights(complex.clone(), &[1.0, 0.5, 2.0]);
println!("Weighted sheaf stalk dim: {}", weighted.stalk_dim());

// Access restriction maps
let map_a = sheaf.edge_restriction_a(0);
let map_b = sheaf.edge_restriction_b(0);
println!("Edge 0 restriction maps:");
println!("  to vertex a: {:?}", map_a);
println!("  to vertex b: {:?}", map_b);
```

### Sheaf Laplacian

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};
use nalgebra::DVector;

// Build a triangle
let mut complex = SimplicialComplex::new();
complex.add_edge(0, 1);
complex.add_edge(1, 2);

// Constant sheaf with R² stalks
let sheaf = CellularSheaf::constant(complex, 2);

// Build the sheaf Laplacian
let lap = SheafLaplacian::from_sheaf(&sheaf);
println!("Laplacian dimension: {} (3 vertices × 2 stalk)", lap.dimension);

// Eigenvalues reveal the spectral structure
let eigenvalues = lap.eigenvalues();
println!("Eigenvalues: {:?}", eigenvalues);

// Fiedler value: algebraic connectivity
let fiedler = lap.fiedler_value();
println!("Fiedler value: {:.4}", fiedler);

// Multiply by a cochain
let cochain = DVector::from_vec(vec![1.0, 0.0, 0.0, 1.0, 0.0, 1.0]);
let result = lap.mul_vec(&cochain);
println!("L·f = {:?}", result);

// Trace and Frobenius norm
println!("Trace: {:.4}", lap.trace());
println!("Frobenius norm: {:.4}", lap.frobenius_norm());

// Compare with standard graph Laplacian
let graph_lap = SheafLaplacian::graph_laplacian(3, &[(0, 1), (1, 2)]);
let graph_eigenvalues = graph_lap.eigenvalues();
println!("Graph Laplacian eigenvalues: {:?}", graph_eigenvalues);
```

### Edge and Triangle Laplacians

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};

// Build a triangle
let mut complex = SimplicialComplex::new();
complex.add_triangle(0, 1, 2);

let sheaf = CellularSheaf::constant(complex, 2);

// Per-edge Laplacian
let edge_lap = SheafLaplacian::edge_laplacian(&sheaf, 0);
println!("Edge 0 Laplacian: {}×{}", edge_lap.dimension, edge_lap.dimension);

// Per-triangle Laplacian
let tri_lap = SheafLaplacian::triangle_laplacian(&sheaf, 0);
println!("Triangle 0 Laplacian: {}×{}", tri_lap.dimension, tri_lap.dimension);

// Eigenvalues of local Laplacians reveal local consistency
let edge_eigenvalues = edge_lap.eigenvalues();
println!("Edge Laplacian eigenvalues: {:?}", edge_eigenvalues);
```

### Sheaf Filtration

```rust
use persistent_sheaf::Filtration;

// Build a filtration from a point cloud
let points = vec![
    vec![0.0, 0.0],
    vec![1.0, 0.0],
    vec![0.0, 1.0],
    vec![3.0, 0.0],
];

let filtration = Filtration::from_point_cloud(&points, 8);

// Create a sheaf at every filtration step
let sheaf_filtration = filtration.make_sheaf_filtration(2);
println!("Sheaf filtration: {} steps", sheaf_filtration.len());

for (i, sheaf) in sheaf_filtration.iter().enumerate() {
    println!("Step {}: {} vertices, {} edges, stalk_dim={}",
        i,
        sheaf.complex.num_simplices(0),
        sheaf.complex.num_simplices(1),
        sheaf.stalk_dimension);
}
```

### Simplicial Maps and Pushforward

```rust
use persistent_sheaf::pushforward::SimplicialMap;

// Map from a 4-vertex complex to a 3-vertex complex
let map = SimplicialMap::new(vec![0, 1, 2, 1]); // vertex 3 maps to 1

// Check if this is simplicial
let source_edges = &[(0, 1), (1, 2), (2, 3)];
let target_edges = &[(0, 1), (1, 2)];
let is_simplicial = map.is_simplicial(source_edges, target_edges);
println!("Is simplicial: {}", is_simplicial);

// Image of an edge
let image = map.edge_image(2, 3);
println!("Edge (2,3) maps to: {:?}", image); // Some((1, 2))

// Collapsing edge: both vertices map to same target
let collapsing = SimplicialMap::new(vec![0, 0]);
let collapsed = collapsing.edge_image(0, 1);
println!("Collapsing edge: {:?}", collapsed); // None

// Identity map
let identity = SimplicialMap::identity(5);
assert_eq!(identity.vertex_map, vec![0, 1, 2, 3, 4]);
```

### Spatial Analysis — Room Topology

```rust
use persistent_sheaf::spatial_analysis::{
    analyze_connectivity, connected_components, room_distance, find_isolated_rooms,
    RoomTopologySheaf,
};

// Room adjacency: 4 rooms, connected as a square
let adjacency = vec![
    vec![1, 3],    // room 0 connects to 1 and 3
    vec![0, 2],    // room 1 connects to 0 and 2
    vec![1, 3],    // room 2 connects to 1 and 3
    vec![0, 2],    // room 3 connects to 0 and 2
];

// Full connectivity analysis
let report = analyze_connectivity(adjacency.clone());
println!("Rooms: {}", report.n_rooms);
println!("Doors: {}", report.n_doors);
println!("Components: {}", report.n_components);
println!("Loops: {}", report.n_loops); // 1 (the square)

// Isolated rooms
let isolated = find_isolated_rooms(adjacency.clone());
println!("Isolated rooms: {:?}", isolated); // none

// Connected components
let components = connected_components(&adjacency);
println!("Components: {:?}", components); // [[0, 1, 2, 3]]

// Shortest path between rooms
let dist = room_distance(&adjacency, 0, 2);
println!("Distance room 0 → 2: {:?}", dist); // Some(2)

// Create a sheaf on the room topology
let mut sheaf = RoomTopologySheaf::from_adjacency(&adjacency, 2);
sheaf.set_stalk(0, vec![1.0, 0.0]); // room 0 has data [1, 0]
sheaf.set_stalk(1, vec![0.0, 1.0]); // room 1 has data [0, 1]

// Restrict stalk from room 0 to room 1 (identity map)
let restricted = sheaf.restrict(0, 1);
println!("Restricted from room 0 to 1: {:?}", restricted); // Some([1.0, 0.0])
```

### Hodge Laplacian and Harmonic Sections

```rust
use persistent_sheaf::{SimplicialComplex, CellularSheaf, SheafLaplacian};
use persistent_sheaf::sheaf_laplacian::{HodgeLaplacian, harmonic_0_sections, harmonic_1_sections};

// Build a sheaf on a line graph
let mut complex = SimplicialComplex::new();
complex.add_edge(0, 1);
complex.add_edge(1, 2);
complex.add_edge(2, 3);

let sheaf = CellularSheaf::constant(complex, 1);

// Hodge Laplacian
let hodge = HodgeLaplacian::from_sheaf(&sheaf);
println!("Hodge Laplacian eigenvalues: {:?}", hodge.eigenvalues());

// Harmonic 0-sections: global sections of the sheaf
let h0 = harmonic_0_sections(&sheaf);
println!("Harmonic 0-sections dimension: {}", h0);

// Harmonic 1-sections: cocycles
let h1 = harmonic_1_sections(&sheaf);
println!("Harmonic 1-sections dimension: {}", h1);
```

## API Reference

### SimplicialComplex

| Method | Description |
|--------|-------------|
| `new()` | Empty complex |
| `add_vertex(v)` | Add 0-simplex |
| `add_edge(a, b)` | Add 1-simplex |
| `add_triangle(a, b, c)` | Add 2-simplex |
| `add_tetrahedron(a, b, c, d)` | Add 3-simplex |
| `vietoris_rips(points, ε)` | VR complex from point cloud |
| `betti_numbers()` | $\beta_0, \beta_1$ |
| `euler_characteristic()` | $\chi = V - E + F - T$ |
| `connected_components()` | Number of components |
| `h0_betti(points, ε)` | $\beta_0$ at threshold $\varepsilon$ |
| `h1_betti(points, ε)` | $\beta_1$ at threshold $\varepsilon$ |

### Filtration

| Method | Description |
|--------|-------------|
| `from_point_cloud(points, steps)` | VR filtration |
| `from_distance_matrix(dist, steps)` | Filtration from distances |
| `compute_persistence()` | Persistence diagram |
| `make_sheaf_filtration(dim)` | Sheaves at each step |
| `thresholds()` | All threshold values |

### PersistenceDiagram

| Method | Description |
|--------|-------------|
| `add(birth, death, dim)` | Add a birth-death pair |
| `filter_dimension(dim)` | Filter by homology dimension |
| `bottleneck_distance(other)` | $d_B$ between diagrams |
| `total_persistence(p)` | $\sum_i (d_i - b_i)^p$ |
| `most_persistent()` | Longest-lived feature |
| `betti_curve(thresholds)` | Alive count at each $t$ |
| `prune(min_persistence)` | Remove short-lived features |

### CellularSheaf

| Method | Description |
|--------|-------------|
| `constant(complex, dim)` | Identity restriction maps |
| `from_weights(complex, w)` | Scalar weight maps |
| `global_section_dimension()` | Dim of global sections |
| `cohomology_dimension(deg)` | $\dim H^k(F)$ |

### SheafLaplacian

| Method | Description |
|--------|-------------|
| `from_sheaf(sheaf)` | Build from sheaf |
| `graph_laplacian(n, edges)` | Standard $L = D - A$ |
| `edge_laplacian(sheaf, idx)` | Per-edge $2d \times 2d$ |
| `triangle_laplacian(sheaf, idx)` | Per-triangle $3d \times 3d$ |
| `eigenvalues()` | Full spectrum |
| `fiedler_value()` | Second-smallest eigenvalue |
| `largest_eigenvalue(iter)` | Power iteration |

## Why This Matters for Agent Systems

1. **Connectivity detection**: Betti numbers tell you when the fleet is fragmented ($\beta_0 > 1$) or has communication loops ($\beta_1 > 0$).
2. **Scale-invariant analysis**: Persistence diagrams reveal structure at all scales simultaneously — no need to pick a single threshold.
3. **Data-aware topology**: Sheaves carry actual data (sensor readings, agent states) on topological structures. The sheaf Laplacian reveals inconsistencies.
4. **Harmonic sections**: The kernel of the sheaf Laplacian gives the globally consistent state assignments — essential for consensus algorithms.
5. **Drift detection**: The bottleneck distance between persistence diagrams detects when the fleet's topological structure has qualitatively changed.

## Integration

### With `wasserstein-agents-rs`

```rust
// Agent distributions become point clouds
// Vietoris-Rips complexes reveal fleet topology
// Persistence diagrams track topological changes over time
```

### With `conservation-law`

```rust
// Topological invariants (Betti numbers, Euler characteristic) are conserved
// under valid fleet reconfigurations
```

### With `si-fleet-api`

```rust
// Betti number monitoring triggers alerts when connectivity changes
// Sheaf Laplacian eigenvalues track network health
```

## License

MIT
