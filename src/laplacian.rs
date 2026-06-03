//! Sheaf Laplacian: generalizes graph Laplacian with sheaf-theoretic information.

use nalgebra::{DMatrix, DVector, SymmetricEigen};

use crate::sheaf::CellularSheaf;

/// The sheaf Laplacian L_F for a cellular sheaf F.
///
/// For a 0-cochain f (assignment of vectors to vertices):
/// (L_F f)(v) = Σ_{v~w} F_{v←w}^T (F_{v←w} f(v) - F_{w←v} f(w))
///
/// This generalizes the graph Laplacian: when all stalks are R and all maps
/// are 1, it reduces to the standard graph Laplacian.
pub struct SheafLaplacian {
    /// Dimension: n_vertices * stalk_dim.
    pub dimension: usize,
    /// The Laplacian matrix.
    pub matrix: DMatrix<f64>,
}

impl SheafLaplacian {
    /// Build the sheaf Laplacian from a cellular sheaf.
    pub fn from_sheaf(sheaf: &CellularSheaf) -> Self {
        let n_verts = sheaf.complex.vertices.len();
        let d = sheaf.stalk_dimension;
        let dim = n_verts * d;
        let mut mat = DMatrix::zeros(dim, dim);

        for (edge_idx, (a, b)) in sheaf.complex.edges.iter().enumerate() {
            if edge_idx >= sheaf.restriction_maps.len() {
                break;
            }
            let (map_to_a, map_to_b) = &sheaf.restriction_maps[edge_idx];

            // Diagonal block for vertex a: F_{a←b}^T * F_{a←b}
            let term_aa = map_to_a.transpose() * map_to_a;
            for i in 0..d {
                for j in 0..d {
                    mat[(a * d + i, a * d + j)] += term_aa[(i, j)];
                }
            }

            // Diagonal block for vertex b: F_{b←a}^T * F_{b←a}
            let term_bb = map_to_b.transpose() * map_to_b;
            for i in 0..d {
                for j in 0..d {
                    mat[(b * d + i, b * d + j)] += term_bb[(i, j)];
                }
            }

            // Off-diagonal: -F_{a←b}^T * F_{b←a}
            let off = map_to_a.transpose() * map_to_b;
            for i in 0..d {
                for j in 0..d {
                    mat[(a * d + i, b * d + j)] -= off[(i, j)];
                    mat[(b * d + j, a * d + i)] -= off[(j, i)];
                }
            }
        }

        Self {
            dimension: dim,
            matrix: mat,
        }
    }

    /// Build the standard graph Laplacian (trivial sheaf: all stalks R, maps 1).
    pub fn graph_laplacian(n_verts: usize, edges: &[(usize, usize)]) -> Self {
        let mut mat = DMatrix::zeros(n_verts, n_verts);
        for &(a, b) in edges {
            mat[(a, a)] += 1.0;
            mat[(b, b)] += 1.0;
            mat[(a, b)] -= 1.0;
            mat[(b, a)] -= 1.0;
        }
        Self {
            dimension: n_verts,
            matrix: mat,
        }
    }

    /// Build the per-edge sheaf Laplacian for a single edge.
    /// Returns the 2d × 2d Laplacian restricted to vertices of the edge.
    pub fn edge_laplacian(sheaf: &CellularSheaf, edge_idx: usize) -> Self {
        if edge_idx >= sheaf.complex.edges.len() {
            return Self {
                dimension: 0,
                matrix: DMatrix::zeros(0, 0),
            };
        }
        let (_a, _b) = sheaf.complex.edges[edge_idx];
        let d = sheaf.stalk_dimension;
        let (map_to_a, map_to_b) = &sheaf.restriction_maps[edge_idx];

        let dim = 2 * d;
        let mut mat = DMatrix::zeros(dim, dim);

        // Block (a,a): F_{a←b}^T * F_{a←b}
        let term_aa = map_to_a.transpose() * map_to_a;
        for i in 0..d {
            for j in 0..d {
                mat[(i, j)] += term_aa[(i, j)];
            }
        }

        // Block (b,b): F_{b←a}^T * F_{b←a}
        let term_bb = map_to_b.transpose() * map_to_b;
        for i in 0..d {
            for j in 0..d {
                mat[(d + i, d + j)] += term_bb[(i, j)];
            }
        }

        // Off-diagonal
        let off = map_to_a.transpose() * map_to_b;
        for i in 0..d {
            for j in 0..d {
                mat[(i, d + j)] -= off[(i, j)];
                mat[(d + j, i)] -= off[(j, i)];
            }
        }

        Self {
            dimension: dim,
            matrix: mat,
        }
    }

    /// Build the per-triangle sheaf Laplacian for a single triangle.
    /// Returns the 3d × 3d Laplacian restricted to vertices of the triangle.
    pub fn triangle_laplacian(sheaf: &CellularSheaf, tri_idx: usize) -> Self {
        if tri_idx >= sheaf.complex.triangles.len() {
            return Self {
                dimension: 0,
                matrix: DMatrix::zeros(0, 0),
            };
        }
        let (a, b, c) = sheaf.complex.triangles[tri_idx];
        let d = sheaf.stalk_dimension;
        let dim = 3 * d;

        // Find edge indices for the three edges of this triangle
        let edges = [(a.min(b), a.max(b)), (b.min(c), b.max(c)), (a.min(c), a.max(c))];
        let mut edge_indices = [usize::MAX, usize::MAX, usize::MAX];
        for (i, &(e0, e1)) in edges.iter().enumerate() {
            for (j, &(ee0, ee1)) in sheaf.complex.edges.iter().enumerate() {
                if e0 == ee0 && e1 == ee1 {
                    edge_indices[i] = j;
                    break;
                }
            }
        }

        let mut mat = DMatrix::zeros(dim, dim);

        for &ei in &edge_indices {
            if ei == usize::MAX || ei >= sheaf.restriction_maps.len() {
                continue;
            }
            let (ea, eb) = sheaf.complex.edges[ei];
            let (map_a, map_b) = &sheaf.restriction_maps[ei];

            // Map vertex indices to local positions (0=a, 1=b, 2=c)
            let local_a = if ea == a { 0 } else if ea == b { 1 } else { 2 };
            let local_b = if eb == a { 0 } else if eb == b { 1 } else { 2 };

            let term_aa = map_a.transpose() * map_a;
            let term_bb = map_b.transpose() * map_b;
            let off = map_a.transpose() * map_b;

            for i in 0..d {
                for j in 0..d {
                    mat[(local_a * d + i, local_a * d + j)] += term_aa[(i, j)];
                    mat[(local_b * d + i, local_b * d + j)] += term_bb[(i, j)];
                    mat[(local_a * d + i, local_b * d + j)] -= off[(i, j)];
                    mat[(local_b * d + j, local_a * d + i)] -= off[(j, i)];
                }
            }
        }

        Self {
            dimension: dim,
            matrix: mat,
        }
    }

    /// Multiply by a vector.
    pub fn mul_vec(&self, v: &DVector<f64>) -> DVector<f64> {
        &self.matrix * v
    }

    /// Compute eigenvalues via full symmetric eigen decomposition (small matrices).
    pub fn eigenvalues(&self) -> Vec<f64> {
        if self.dimension == 0 {
            return vec![];
        }
        let sym = SymmetricEigen::new(self.matrix.clone());
        let mut eigenvalues = sym.eigenvalues.iter().copied().collect::<Vec<_>>();
        eigenvalues.sort_by(|a, b| a.partial_cmp(b).unwrap());
        eigenvalues
    }

    /// Largest eigenvalue via power iteration.
    pub fn largest_eigenvalue(&self, iterations: usize) -> f64 {
        let n = self.dimension;
        if n == 0 {
            return 0.0;
        }
        let mut v = DVector::from_iterator(n, (0..n).map(|i| (i as f64 + 1.0) / n as f64));
        let norm = v.norm();
        if norm < 1e-15 {
            return 0.0;
        }
        v /= norm;

        for _ in 0..iterations {
            let mv = &self.matrix * &v;
            let norm = mv.norm();
            if norm < 1e-15 {
                return 0.0;
            }
            v = mv / norm;
        }

        let mv = &self.matrix * &v;
        v.dot(&mv)
    }

    /// Fiedler value (second smallest eigenvalue, algebraic connectivity).
    pub fn fiedler_value(&self) -> f64 {
        if self.dimension < 2 {
            return 0.0;
        }
        let ev = self.eigenvalues();
        if ev.len() >= 2 {
            ev[1].max(0.0)
        } else {
            0.0
        }
    }

    /// Trace of the Laplacian.
    pub fn trace(&self) -> f64 {
        (0..self.dimension).map(|i| self.matrix[(i, i)]).sum()
    }

    /// Frobenius norm of the Laplacian.
    pub fn frobenius_norm(&self) -> f64 {
        (0..self.dimension)
            .flat_map(|i| (0..self.dimension).map(move |j| self.matrix[(i, j)]))
            .map(|x| x * x)
            .sum::<f64>()
            .sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simplicial::SimplicialComplex;

    #[test]
    fn test_graph_laplacian() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        assert_eq!(l.dimension, 3);
        assert_eq!(l.matrix[(0, 0)], 1.0);
        assert_eq!(l.matrix[(1, 1)], 2.0);
        assert_eq!(l.matrix[(2, 2)], 1.0);
    }

    #[test]
    fn test_graph_laplacian_single_edge() {
        let edges = vec![(0, 1)];
        let l = SheafLaplacian::graph_laplacian(2, &edges);
        assert_eq!(l.matrix[(0, 0)], 1.0);
        assert_eq!(l.matrix[(1, 1)], 1.0);
        assert_eq!(l.matrix[(0, 1)], -1.0);
        assert_eq!(l.matrix[(1, 0)], -1.0);
    }

    #[test]
    fn test_laplacian_positive_semidefinite() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        let eig = l.largest_eigenvalue(100);
        assert!(eig > 0.0);
    }

    #[test]
    fn test_sheaf_laplacian() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 2);
        let l = SheafLaplacian::from_sheaf(&sheaf);
        assert_eq!(l.dimension, 4);
    }

    #[test]
    fn test_sheaf_laplacian_two_edges() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        let sheaf = CellularSheaf::constant(c, 1);
        let l = SheafLaplacian::from_sheaf(&sheaf);
        assert_eq!(l.dimension, 3);
    }

    #[test]
    fn test_mul_vec() {
        let edges = vec![(0, 1)];
        let l = SheafLaplacian::graph_laplacian(2, &edges);
        let v = DVector::from_vec(vec![1.0, -1.0]);
        let lv = l.mul_vec(&v);
        assert!((lv[0] - 2.0).abs() < 0.01);
        assert!((lv[1] + 2.0).abs() < 0.01);
    }

    #[test]
    fn test_eigenvalues() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        let ev = l.eigenvalues();
        assert_eq!(ev.len(), 3);
        assert!(ev[0].abs() < 1e-10); // smallest ≈ 0
    }

    #[test]
    fn test_fiedler_value() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        let f = l.fiedler_value();
        assert!(f > 0.0);
    }

    #[test]
    fn test_trace() {
        let edges = vec![(0, 1), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        let tr = l.trace();
        assert!((tr - 4.0).abs() < 1e-10); // 1+2+1
    }

    #[test]
    fn test_frobenius_norm() {
        let edges = vec![(0, 1)];
        let l = SheafLaplacian::graph_laplacian(2, &edges);
        let f = l.frobenius_norm();
        assert!(f > 0.0);
    }

    #[test]
    fn test_edge_laplacian() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 1);
        let el = SheafLaplacian::edge_laplacian(&sheaf, 0);
        assert_eq!(el.dimension, 2);
    }

    #[test]
    fn test_edge_laplacian_2d_stalk() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 2);
        let el = SheafLaplacian::edge_laplacian(&sheaf, 0);
        assert_eq!(el.dimension, 4);
    }

    #[test]
    fn test_triangle_laplacian() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        let sheaf = CellularSheaf::constant(c, 1);
        let tl = SheafLaplacian::triangle_laplacian(&sheaf, 0);
        assert_eq!(tl.dimension, 3);
    }

    #[test]
    fn test_triangle_laplacian_2d_stalk() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        let sheaf = CellularSheaf::constant(c, 2);
        let tl = SheafLaplacian::triangle_laplacian(&sheaf, 0);
        assert_eq!(tl.dimension, 6);
    }

    #[test]
    fn test_empty_laplacian() {
        let l = SheafLaplacian::graph_laplacian(0, &[]);
        assert_eq!(l.dimension, 0);
        let ev = l.eigenvalues();
        assert!(ev.is_empty());
    }

    #[test]
    fn test_symmetry() {
        let edges = vec![(0, 1), (0, 2), (1, 2)];
        let l = SheafLaplacian::graph_laplacian(3, &edges);
        for i in 0..l.dimension {
            for j in 0..l.dimension {
                assert!((l.matrix[(i, j)] - l.matrix[(j, i)]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_laplacian_kernel_dim() {
        // Two disconnected vertices → nullspace dimension 2 → both eigenvalues ≈ 0
        let edges = vec![];
        let l = SheafLaplacian::graph_laplacian(2, &edges);
        let ev = l.eigenvalues();
        assert_eq!(ev.len(), 2);
        assert!((ev[0]).abs() < 1e-10);
        assert!((ev[1]).abs() < 1e-10, "ev[1] = {} should be ~0 for disconnected graph", ev[1]);
    }
}
