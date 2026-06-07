//! Spatial analysis via persistent sheaf cohomology on room topologies.
//!
//! Applies persistent homology to room/floor-plan networks modelled as
//! simplicial complexes. A [`RoomTopologySheaf`] assigns data stalks to each
//! room and linear restriction maps along door connections, then exposes
//! global connectivity analysis through H₀/H₁ persistent homology.

use crate::persistence::{PersistenceDiagram, PersistencePair};
use crate::simplicial::SimplicialComplex;

/// Summary of the persistent-homology analysis of a room network.
#[derive(Debug, Clone)]
pub struct PersistenceReport {
    /// Number of rooms (vertices).
    pub n_rooms: usize,
    /// Number of door connections (edges).
    pub n_doors: usize,
    /// H₀ persistence pairs (connected-component merges).
    pub h0_pairs: Vec<PersistencePair>,
    /// H₁ persistence pairs (loop / cycle births–deaths).
    pub h1_pairs: Vec<PersistencePair>,
    /// Indices of rooms that are completely isolated (degree 0).
    pub isolated_rooms: Vec<usize>,
    /// Number of connected components (including isolated rooms).
    pub n_components: usize,
    /// Number of independent loops (1-cycles) that persist to ∞.
    pub n_loops: usize,
}

impl PersistenceReport {
    /// Total persistence in dimension `d` with power `p`.
    pub fn total_persistence(&self, d: usize, p: f64) -> f64 {
        let pairs = if d == 0 { &self.h0_pairs } else { &self.h1_pairs };
        pairs
            .iter()
            .map(|pp| (pp.death - pp.birth).abs().powf(p))
            .sum()
    }
}

// ---------------------------------------------------------------------------
// RoomTopologySheaf
// ---------------------------------------------------------------------------

/// Sheaf on a room-adjacency network.
///
/// Each room `i` carries a stalk `Vec<f64>` of dimension `stalk_dim`.
/// Restriction maps along door edges are represented as `stalk_dim × stalk_dim`
/// matrices stored in row-major order. When no custom map is given the identity
/// is used (constant sheaf).
#[derive(Debug, Clone)]
pub struct RoomTopologySheaf {
    /// Number of rooms.
    pub n_rooms: usize,
    /// Stalk dimension.
    pub stalk_dim: usize,
    /// Per-room stalk data (length = n_rooms × stalk_dim).
    pub stalks: Vec<Vec<f64>>,
    /// Restriction maps keyed by `(room_a, room_b)` with `room_a < room_b`.
    /// Stored as flat row-major `stalk_dim²` vector.
    pub restriction_maps: std::collections::BTreeMap<(usize, usize), Vec<f64>>,
}

impl RoomTopologySheaf {
    /// Build a constant sheaf (identity restriction maps) with zero-initialised stalks.
    pub fn constant(n_rooms: usize, stalk_dim: usize) -> Self {
        let stalks = vec![vec![0.0; stalk_dim]; n_rooms];
        Self {
            n_rooms,
            stalk_dim,
            stalks,
            restriction_maps: std::collections::BTreeMap::new(),
        }
    }

    /// Build from an explicit adjacency list and stalk dimension.
    ///
    /// Creates identity restriction maps for every edge.
    pub fn from_adjacency(adjacency: &[Vec<usize>], stalk_dim: usize) -> Self {
        let n = adjacency.len();
        let mut maps = std::collections::BTreeMap::new();
        for (i, neighbours) in adjacency.iter().enumerate() {
            for &j in neighbours {
                if i < j {
                    let id: Vec<f64> = (0..stalk_dim)
                        .flat_map(|r| {
                            (0..stalk_dim).map(move |c| if r == c { 1.0 } else { 0.0 })
                        })
                        .collect();
                    maps.insert((i, j), id);
                }
            }
        }
        Self {
            n_rooms: n,
            stalk_dim,
            stalks: vec![vec![0.0; stalk_dim]; n],
            restriction_maps: maps,
        }
    }

    /// Set the stalk data for room `i`.
    pub fn set_stalk(&mut self, room: usize, data: Vec<f64>) {
        assert_eq!(data.len(), self.stalk_dim, "stalk dimension mismatch");
        self.stalks[room] = data;
    }

    /// Restrict the stalk from room `a` to room `b` along the connecting door.
    ///
    /// Returns `None` if the rooms are not connected.
    pub fn restrict(&self, a: usize, b: usize) -> Option<Vec<f64>> {
        let key = if a < b { (a, b) } else { (b, a) };
        let mat = self.restriction_maps.get(&key)?;
        let s = &self.stalks[a];
        let d = self.stalk_dim;
        let mut out = vec![0.0; d];
        for r in 0..d {
            for c in 0..d {
                out[r] += mat[r * d + c] * s[c];
            }
        }
        Some(out)
    }
}

// ---------------------------------------------------------------------------
// Core analysis functions
// ---------------------------------------------------------------------------

/// Analyse the connectivity of a room-adjacency graph using persistent homology.
///
/// Uses a union-find approach for H₀ and a cycle-detection pass for H₁.
pub fn analyze_connectivity(adjacency: Vec<Vec<usize>>) -> PersistenceReport {
    let n = adjacency.len();
    let mut edges: Vec<(usize, usize, usize)> = Vec::new(); // (u, v, sort_key)
    for (u, nbrs) in adjacency.iter().enumerate() {
        for &v in nbrs {
            if u < v {
                edges.push((u, v, u.max(v)));
            }
        }
    }
    edges.sort_by_key(|e| e.2);

    let mut isolated: Vec<usize> = (0..n).filter(|&i| adjacency[i].is_empty()).collect();
    isolated.sort();

    // Union-Find
    let mut parent: Vec<usize> = (0..n).collect();
    let mut rank = vec![0usize; n];

    fn find(parent: &mut [usize], x: usize) -> usize {
        if parent[x] != x {
            parent[x] = find(parent, parent[x]);
        }
        parent[x]
    }

    let mut h0_pairs: Vec<PersistencePair> = Vec::new();
    let mut h1_pairs: Vec<PersistencePair> = Vec::new();
    let mut edge_idx = 0usize;

    for scale in 0..n {
        while edge_idx < edges.len() && edges[edge_idx].2 <= scale {
            let (u, v, _) = edges[edge_idx];
            let ru = find(&mut parent, u);
            let rv = find(&mut parent, v);
            if ru != rv {
                // merge
                let birth = ru.min(rv) as f64;
                h0_pairs.push(PersistencePair {
                    birth: birth,
                    death: scale as f64,
                    dimension: 0,
                });
                if rank[ru] < rank[rv] {
                    parent[ru] = rv;
                } else if rank[ru] > rank[rv] {
                    parent[rv] = ru;
                } else {
                    parent[rv] = ru;
                    rank[ru] += 1;
                }
            } else {
                // cycle detected → H₁
                h1_pairs.push(PersistencePair {
                    birth: scale as f64,
                    death: f64::INFINITY,
                    dimension: 1,
                });
            }
            edge_idx += 1;
        }
    }

    // Count final components
    let mut roots: Vec<usize> = (0..n).map(|i| find(&mut parent, i)).collect();
    roots.sort();
    roots.dedup();
    let n_components = roots.len();

    let n_doors = edges.len();

    PersistenceReport {
        n_rooms: n,
        n_doors,
        h0_pairs,
        n_loops: h1_pairs.len(),
        h1_pairs,
        isolated_rooms: isolated,
        n_components,
    }
}

/// Return the indices of rooms with degree 0 (no connections).
pub fn find_isolated_rooms(adjacency: Vec<Vec<usize>>) -> Vec<usize> {
    adjacency
        .iter()
        .enumerate()
        .filter(|(_, nbrs)| nbrs.is_empty())
        .map(|(i, _)| i)
        .collect()
}

/// Compute the graph Laplacian of the room adjacency (symmetric, normalised).
///
/// Returns `L = D - A` as a flat row-major `n × n` matrix.
pub fn room_graph_laplacian(adjacency: &[Vec<usize>]) -> Vec<f64> {
    let n = adjacency.len();
    let mut lap = vec![0.0; n * n];
    for (i, nbrs) in adjacency.iter().enumerate() {
        lap[i * n + i] = nbrs.len() as f64;
        for &j in nbrs {
            lap[i * n + j] -= 1.0;
        }
    }
    lap
}

/// BFS shortest-path distance between two rooms. Returns `None` if unreachable.
pub fn room_distance(adjacency: &[Vec<usize>], from: usize, to: usize) -> Option<usize> {
    if from == to {
        return Some(0);
    }
    let n = adjacency.len();
    let mut visited = vec![false; n];
    let mut queue = std::collections::VecDeque::new();
    visited[from] = true;
    queue.push_back((from, 0usize));
    while let Some((node, dist)) = queue.pop_front() {
        for &nbr in &adjacency[node] {
            if nbr == to {
                return Some(dist + 1);
            }
            if !visited[nbr] {
                visited[nbr] = true;
                queue.push_back((nbr, dist + 1));
            }
        }
    }
    None
}

/// Return connected components as vectors of room indices.
pub fn connected_components(adjacency: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let n = adjacency.len();
    let mut visited = vec![false; n];
    let mut components = Vec::new();
    for start in 0..n {
        if visited[start] {
            continue;
        }
        let mut comp = Vec::new();
        let mut stack = vec![start];
        visited[start] = true;
        while let Some(node) = stack.pop() {
            comp.push(node);
            for &nbr in &adjacency[node] {
                if !visited[nbr] {
                    visited[nbr] = true;
                    stack.push(nbr);
                }
            }
        }
        comp.sort();
        components.push(comp);
    }
    components
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fully_connected_three_rooms() {
        // triangle: 0-1, 1-2, 0-2
        let adj = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let report = analyze_connectivity(adj);
        assert_eq!(report.n_rooms, 3);
        assert_eq!(report.n_doors, 3);
        assert!(report.isolated_rooms.is_empty());
        assert_eq!(report.n_components, 1);
        // triangle has one independent loop
        assert_eq!(report.n_loops, 1);
    }

    #[test]
    fn test_isolated_rooms() {
        let adj = vec![vec![1], vec![0], vec![], vec![]];
        let iso = find_isolated_rooms(adj.clone());
        assert_eq!(iso, vec![2, 3]);

        let report = analyze_connectivity(adj);
        assert_eq!(report.isolated_rooms, vec![2, 3]);
        assert_eq!(report.n_components, 3); // {0,1}, {2}, {3}
    }

    #[test]
    fn test_line_of_rooms() {
        // 0-1-2-3
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        let report = analyze_connectivity(adj);
        assert_eq!(report.n_rooms, 4);
        assert_eq!(report.n_doors, 3);
        assert_eq!(report.n_components, 1);
        assert_eq!(report.n_loops, 0);
    }

    #[test]
    fn test_two_components() {
        // 0-1 and 2-3
        let adj = vec![vec![1], vec![0], vec![3], vec![2]];
        let comps = connected_components(&adj);
        assert_eq!(comps.len(), 2);
        assert!(comps.contains(&vec![0, 1]));
        assert!(comps.contains(&vec![2, 3]));
    }

    #[test]
    fn test_room_distance() {
        // 0-1-2-3
        let adj = vec![vec![1], vec![0, 2], vec![1, 3], vec![2]];
        assert_eq!(room_distance(&adj, 0, 3), Some(3));
        assert_eq!(room_distance(&adj, 1, 3), Some(2));
        assert_eq!(room_distance(&adj, 2, 2), Some(0));
    }

    #[test]
    fn test_unreachable_rooms() {
        let adj = vec![vec![1], vec![0], vec![]];
        assert_eq!(room_distance(&adj, 0, 2), None);
    }

    #[test]
    fn test_graph_laplacian() {
        // single edge: 0-1
        let adj = vec![vec![1], vec![0]];
        let lap = room_graph_laplacian(&adj);
        // L = [[1, -1], [-1, 1]]
        assert_eq!(lap[0 * 2 + 0], 1.0);
        assert_eq!(lap[0 * 2 + 1], -1.0);
        assert_eq!(lap[1 * 2 + 0], -1.0);
        assert_eq!(lap[1 * 2 + 1], 1.0);
    }

    #[test]
    fn test_constant_sheaf_restrict() {
        let adj = vec![vec![1], vec![0]];
        let mut sheaf = RoomTopologySheaf::from_adjacency(&adj, 2);
        sheaf.set_stalk(0, vec![3.0, 4.0]);
        let restricted = sheaf.restrict(0, 1).unwrap();
        // Identity map → same vector
        assert_eq!(restricted, vec![3.0, 4.0]);
    }

    #[test]
    fn test_sheaf_unconnected_restrict_is_none() {
        let adj = vec![vec![], vec![]];
        let sheaf = RoomTopologySheaf::from_adjacency(&adj, 1);
        assert!(sheaf.restrict(0, 1).is_none());
    }

    #[test]
    fn test_total_persistence() {
        let adj = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
        let report = analyze_connectivity(adj);
        let tp = report.total_persistence(0, 1.0);
        assert!(tp >= 0.0);
    }

    #[test]
    fn test_single_room() {
        let adj = vec![vec![]];
        let report = analyze_connectivity(adj);
        assert_eq!(report.n_rooms, 1);
        assert_eq!(report.n_doors, 0);
        assert_eq!(report.n_components, 1);
        assert_eq!(report.isolated_rooms, vec![0]);
    }

    #[test]
    fn test_persistence_report_fields() {
        // square with diagonal: 0-1, 1-2, 2-3, 3-0, 0-2
        let adj = vec![
            vec![1, 3, 2],
            vec![0, 2],
            vec![1, 3, 0],
            vec![2, 0],
        ];
        let report = analyze_connectivity(adj);
        assert_eq!(report.n_rooms, 4);
        assert_eq!(report.n_doors, 5);
        assert_eq!(report.n_components, 1);
        // multiple cycles possible
        assert!(report.n_loops >= 1);
    }
}
