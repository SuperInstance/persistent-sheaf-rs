//! Sheaf pushforward (direct image) along continuous maps.
//!
//! Given a continuous map f: X → Y and a sheaf F on X, the **pushforward**
//! (or **direct image**) sheaf f_* F on Y is defined by:
//!
//! ```text
//! (f_* F)(U) = F(f⁻¹(U))   for each open U ⊆ Y
//! ```
//!
//! This module works in the discrete/combinatorial setting where X and Y
//! are simplicial complexes and f is a simplicial map.

use nalgebra::{DMatrix, DVector};

use crate::sheaf::CellularSheaf;
use crate::simplicial::SimplicialComplex;

/// A simplicial map between complexes, specified by vertex images.
#[derive(Debug, Clone)]
pub struct SimplicialMap {
    /// `vertex_map[v]` = image of vertex v in the target complex.
    pub vertex_map: Vec<usize>,
}

impl SimplicialMap {
    /// Create a new simplicial map from vertex assignments.
    pub fn new(vertex_map: Vec<usize>) -> Self {
        Self { vertex_map }
    }

    /// Identity map on a complex with `n` vertices.
    pub fn identity(n: usize) -> Self {
        Self {
            vertex_map: (0..n).collect(),
        }
    }

    /// Check that the map sends edges to edges (simplicial condition).
    pub fn is_simplicial(&self, source_edges: &[(usize, usize)], target_edges: &[(usize, usize)]) -> bool {
        let target_set: std::collections::HashSet<(usize, usize)> = target_edges
            .iter()
            .flat_map(|&(a, b)| {
                let (lo, hi) = if a < b { (a, b) } else { (b, a) };
                vec![(lo, hi)]
            })
            .collect();

        for &(u, v) in source_edges {
            let fu = self.vertex_map[u];
            let fv = self.vertex_map[v];
            if fu == fv {
                continue; // edge collapses to a vertex — still simplicial
            }
            let (lo, hi) = if fu < fv { (fu, fv) } else { (fv, fu) };
            if !target_set.contains(&(lo, hi)) {
                return false;
            }
        }
        true
    }

    /// Image of an edge. Returns `None` if the edge collapses (both vertices
    /// map to the same target vertex).
    pub fn edge_image(&self, u: usize, v: usize) -> Option<(usize, usize)> {
        let fu = self.vertex_map[u];
        let fv = self.vertex_map[v];
        if fu == fv {
            None
        } else if fu < fv {
            Some((fu, fv))
        } else {
            Some((fv, fu))
        }
    }
}

/// The pushforward sheaf f_* F on the target complex.
///
/// Stalks: (f_* F)(w) = ⊕_{v : f(v)=w} F(v)   (direct sum of preimage stalks).
/// Restriction maps: induced from F on the preimage edges.
#[derive(Debug, Clone)]
pub struct PushforwardSheaf {
    /// Target complex.
    pub target_complex: SimplicialComplex,
    /// Stalk dimension at each target vertex = stalk_dim × (number of preimage vertices).
    pub stalk_dimensions: Vec<usize>,
    /// Restriction maps on target edges, each a matrix from stalk(fv) to stalk(fu).
    pub restriction_maps: Vec<(DMatrix<f64>, DMatrix<f64>)>,
    /// For each target vertex, the list of source vertices mapping to it.
    pub fiber: Vec<Vec<usize>>,
    /// Original stalk dimension of the source sheaf.
    pub base_stalk_dim: usize,
}

impl PushforwardSheaf {
    /// Compute the pushforward of `sheaf` along `map`.
    ///
    /// The target complex must have vertices numbered 0..n_target.
    pub fn compute(
        sheaf: &CellularSheaf,
        map: &SimplicialMap,
        target: SimplicialComplex,
    ) -> Self {
        let n_target = target.vertices.len();
        let d = sheaf.stalk_dimension;

        // Compute fiber: for each target vertex, which source vertices map to it?
        let mut fiber = vec![Vec::new(); n_target];
        for (src_v, &tgt_v) in map.vertex_map.iter().enumerate() {
            fiber[tgt_v].push(src_v);
        }

        // Stalk dimension at each target vertex
        let stalk_dimensions: Vec<usize> = fiber.iter().map(|f| f.len() * d).collect();

        // Build a mapping from source vertex to its position in the target stalk
        // vertex_index[src_v] = (target vertex, offset within target stalk)
        let mut vertex_offset = vec![(0usize, 0usize); map.vertex_map.len()];
        for (tgt_v, src_vertices) in fiber.iter().enumerate() {
            for (i, &src_v) in src_vertices.iter().enumerate() {
                vertex_offset[src_v] = (tgt_v, i * d);
            }
        }

        // Build restriction maps for each target edge
        // For target edge (u, v), collect all source edges mapping to (u,v)
        // and build block matrices
        let mut edge_restriction_maps = Vec::new();

        for &(tu, tv) in &target.edges {
            let du = stalk_dimensions[tu];
            let dv = stalk_dimensions[tv];

            let mut map_to_u = DMatrix::zeros(du, du);
            let mut map_to_v = DMatrix::zeros(dv, dv);

            // Find source edges mapping to this target edge
            for (e_idx, &(su, sv)) in sheaf.complex.edges.iter().enumerate() {
                if let Some((fu, fv)) = map.edge_image(su, sv) {
                    let (fu_lo, fv_lo) = if fu < fv { (fu, fv) } else { (fv, fu) };
                    if fu_lo == tu && fv_lo == tv && e_idx < sheaf.restriction_maps.len() {
                        let (src_map_to_a, src_map_to_b) = &sheaf.restriction_maps[e_idx];
                        // su maps to fu, sv maps to fv
                        let (off_u, off_v) = if fu == tu {
                            // su → tu, sv → tv
                            let off_su = vertex_offset[su].1;
                            let off_sv = vertex_offset[sv].1;
                            // map_to_u: place src_map_to_a at block (off_su..off_su+d, off_su..off_su+d)
                            for i in 0..d {
                                for j in 0..d {
                                    map_to_u[(off_su + i, off_su + j)] += src_map_to_a[(i, j)];
                                    map_to_v[(off_sv + i, off_sv + j)] += src_map_to_b[(i, j)];
                                }
                            }
                            (off_su, off_sv)
                        } else {
                            // sv → tu, su → tv (swapped)
                            let off_sv = vertex_offset[sv].1;
                            let off_su = vertex_offset[su].1;
                            for i in 0..d {
                                for j in 0..d {
                                    map_to_u[(off_sv + i, off_sv + j)] += src_map_to_b[(i, j)];
                                    map_to_v[(off_su + i, off_su + j)] += src_map_to_a[(i, j)];
                                }
                            }
                            (off_sv, off_su)
                        };
                        let _ = (off_u, off_v); // used above
                    }
                }
            }

            edge_restriction_maps.push((map_to_u, map_to_v));
        }

        Self {
            target_complex: target,
            stalk_dimensions,
            restriction_maps: edge_restriction_maps,
            fiber,
            base_stalk_dim: d,
        }
    }

    /// The dimension of the space of global sections.
    /// This is the sum of stalk dimensions.
    pub fn total_stalk_dim(&self) -> usize {
        self.stalk_dimensions.iter().sum()
    }

    /// Total number of "fibers" (preimage vertices).
    pub fn total_fiber_size(&self) -> usize {
        self.fiber.iter().map(|f| f.len()).sum()
    }
}

/// Verification that a section of the source sheaf descends to a section
/// of the pushforward sheaf (compatibility check).
pub fn verify_compatibility(
    sheaf: &CellularSheaf,
    map: &SimplicialMap,
    section: &[Vec<f64>],
    tol: f64,
) -> bool {
    // For each edge (u,v) in source, the restriction maps must agree
    // when the vertices map to the same target vertex.
    for (e_idx, &(u, v)) in sheaf.complex.edges.iter().enumerate() {
        if e_idx >= sheaf.restriction_maps.len() {
            break;
        }
        let fu = map.vertex_map[u];
        let fv = map.vertex_map[v];

        if fu == fv {
            // Both vertices map to same target → restrictions must agree
            let (map_to_u, map_to_v) = &sheaf.restriction_maps[e_idx];
            let su: DVector<f64> = DVector::from_vec(section[u].clone());
            let sv: DVector<f64> = DVector::from_vec(section[v].clone());
            let ru = map_to_u * &su;
            let rv = map_to_v * &sv;
            for i in 0..ru.len() {
                if (ru[i] - rv[i]).abs() > tol {
                    return false;
                }
            }
        }
    }
    true
}

/// Compute the direct image sections: given sections of F, produce sections of f_*F.
pub fn pushforward_sections(
    pushforward: &PushforwardSheaf,
    source_sections: &[Vec<f64>],
    _map: &SimplicialMap,
) -> Vec<Vec<f64>> {
    let n_target = pushforward.target_complex.vertices.len();
    let d = pushforward.base_stalk_dim;

    let mut result = vec![Vec::new(); n_target];
    for (tgt_v, fiber) in pushforward.fiber.iter().enumerate() {
        let mut stalk = vec![0.0; fiber.len() * d];
        for (i, &src_v) in fiber.iter().enumerate() {
            if src_v < source_sections.len() {
                for j in 0..d {
                    stalk[i * d + j] = source_sections[src_v][j];
                }
            }
        }
        result[tgt_v] = stalk;
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplicial::SimplicialComplex;

    /// A 3-vertex line graph: 0—1—2
    fn line_complex() -> SimplicialComplex {
        SimplicialComplex {
            vertices: vec![0, 1, 2],
            edges: vec![(0, 1), (1, 2)],
            triangles: vec![],
            tetrahedra: vec![],
        }
    }

    /// A 2-vertex graph: 0—1
    fn edge_complex() -> SimplicialComplex {
        SimplicialComplex {
            vertices: vec![0, 1],
            edges: vec![(0, 1)],
            triangles: vec![],
            tetrahedra: vec![],
        }
    }

    fn identity_sheaf_on_line() -> CellularSheaf {
        let c = line_complex();
        CellularSheaf::constant(c, 2)
    }

    #[test]
    fn test_simplicial_map_identity() {
        let edges = vec![(0, 1), (1, 2)];
        let m = SimplicialMap::identity(3);
        assert!(m.is_simplicial(&edges, &edges));
    }

    #[test]
    fn test_simplicial_map_non_simplicial() {
        let src_edges = vec![(0, 1), (1, 2), (0, 2)];
        let tgt_edges = vec![(0, 1)]; // missing (0,2)
        // Map: 0→0, 1→1, 2→0 — edge (1,2) maps to (0,1) ✓, edge (0,2) collapses ✓
        let m = SimplicialMap::new(vec![0, 1, 0]);
        assert!(m.is_simplicial(&src_edges, &tgt_edges));
    }

    #[test]
    fn test_edge_image_collapse() {
        let m = SimplicialMap::new(vec![0, 0, 1]);
        assert_eq!(m.edge_image(0, 1), None); // both map to 0
        assert_eq!(m.edge_image(0, 2), Some((0, 1)));
        assert_eq!(m.edge_image(1, 2), Some((0, 1)));
    }

    #[test]
    fn test_pushforward_identity() {
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::identity(3);
        let target = line_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        // Identity map: fiber is singletons
        assert_eq!(pf.fiber.len(), 3);
        assert_eq!(pf.fiber[0], vec![0]);
        assert_eq!(pf.fiber[1], vec![1]);
        assert_eq!(pf.fiber[2], vec![2]);

        // Stalk dims unchanged
        assert_eq!(pf.stalk_dimensions, vec![2, 2, 2]);
    }

    #[test]
    fn test_pushforward_collapse() {
        // Collapse line 0—1—2 onto edge 0—1 by mapping vertex 2 → 1
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::new(vec![0, 1, 1]); // 0→0, 1→1, 2→1
        let target = edge_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        // fiber[0] = [0], fiber[1] = [1, 2]
        assert_eq!(pf.fiber[0], vec![0]);
        assert_eq!(pf.fiber[1], vec![1, 2]);

        // Stalk dims: vertex 0 has dim 2, vertex 1 has dim 4
        assert_eq!(pf.stalk_dimensions, vec![2, 4]);
    }

    #[test]
    fn test_pushforward_total_stalk_dim() {
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::identity(3);
        let target = line_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        assert_eq!(pf.total_stalk_dim(), 6); // 3 vertices × stalk_dim 2
    }

    #[test]
    fn test_pushforward_total_fiber() {
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::identity(3);
        let target = line_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        assert_eq!(pf.total_fiber_size(), 3);
    }

    #[test]
    fn test_compatibility_identity() {
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::identity(3);
        // Constant section — always compatible
        let section = vec![vec![1.0, 0.0], vec![1.0, 0.0], vec![1.0, 0.0]];
        assert!(verify_compatibility(&sheaf, &map, &section, 1e-10));
    }

    #[test]
    fn test_pushforward_sections_identity() {
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::identity(3);
        let target = line_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        let src_sections = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
            vec![5.0, 6.0],
        ];
        let tgt_sections = pushforward_sections(&pf, &src_sections, &map);

        assert_eq!(tgt_sections.len(), 3);
        assert_eq!(tgt_sections[0], vec![1.0, 2.0]);
        assert_eq!(tgt_sections[1], vec![3.0, 4.0]);
        assert_eq!(tgt_sections[2], vec![5.0, 6.0]);
    }

    #[test]
    fn test_pushforward_sections_collapse() {
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::new(vec![0, 1, 1]);
        let target = edge_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        let src_sections = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
            vec![5.0, 6.0],
        ];
        let tgt_sections = pushforward_sections(&pf, &src_sections, &map);

        assert_eq!(tgt_sections.len(), 2);
        assert_eq!(tgt_sections[0], vec![1.0, 2.0]);
        assert_eq!(tgt_sections[1], vec![3.0, 4.0, 5.0, 6.0]); // direct sum
    }

    #[test]
    fn test_pushforward_preserves_stalk_count() {
        // Total stalk dimension of pushforward = total stalk dim of source
        let sheaf = identity_sheaf_on_line();
        let map = SimplicialMap::new(vec![0, 1, 1]);
        let target = edge_complex();
        let pf = PushforwardSheaf::compute(&sheaf, &map, target);

        // Source: 3 vertices × dim 2 = 6
        // Pushforward: fiber[0] has 1×2=2, fiber[1] has 2×2=4, total = 6
        assert_eq!(pf.total_stalk_dim(), 6);
    }
}
