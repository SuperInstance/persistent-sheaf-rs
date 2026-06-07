//! Sheaf Laplacian via Hodge theory: coboundary maps, Hodge Laplacian,
//! and harmonic sections.
//!
//! The sheaf Laplacian generalises the combinatorial Laplacian by encoding
//! restriction maps.  Given a cellular sheaf **F** on a simplicial complex,
//! we form the **coboundary operator** δ from cochain groups C⁰ → C¹,
//! and the Hodge Laplacian L₁ = δ*δ + δδ* acting on 1-cochains.

use nalgebra::{DMatrix, DVector, SymmetricEigen};

// ── Coboundary maps ─────────────────────────────────────────────────

/// A cochain of dimension-0: assigns a vector to each vertex.
pub type Cochain0 = DVector<f64>;

/// A cochain of dimension-1: assigns a vector to each edge.
pub type Cochain1 = DVector<f64>;

/// Build the **coboundary matrix** δ: C⁰(F) → C¹(F).
///
/// For edge e = (u, v) the row corresponding to e encodes:
///   (δs)(e) = F_{v←e} · s(v) − F_{u←e} · s(u)
///
/// `restriction_maps` is indexed by edge, each entry is
/// `(map_to_u, map_to_v)` where map_to_u : stalk(v) → stalk(u).
/// `n_vertices` and `stalk_dim` determine the cochain dimensions.
pub fn coboundary_matrix(
    n_vertices: usize,
    edges: &[(usize, usize)],
    stalk_dim: usize,
    restriction_maps: &[(DMatrix<f64>, DMatrix<f64>)],
) -> DMatrix<f64> {
    let n_edges = edges.len();
    let rows = n_edges * stalk_dim;
    let cols = n_vertices * stalk_dim;
    let mut delta = DMatrix::zeros(rows, cols);

    for (e_idx, &(u, v)) in edges.iter().enumerate() {
        let (map_to_u, map_to_v) = &restriction_maps[e_idx];
        // row block for edge e
        let row_start = e_idx * stalk_dim;
        // column block for vertex v: +F_{v←e}^T  (but map_to_v maps v→stalk at edge, so direct)
        // We place map_to_v into columns of v, and -map_to_u into columns of u.
        // Actually: δs(e) = F_{v←e} s(v) - F_{u←e} s(u)
        // map_to_u maps from v's stalk to u's stalk? Let's use the convention that
        // restriction_maps[e] = (map_to_a, map_to_b) where a = edges[e].0, b = edges[e].1
        // δs(e) = map_to_b * s(b) - map_to_a * s(a) ... 
        // We use a simpler convention: the row for edge e gets +R_e^b in column v, -R_e^a in column u.
        
        // Place +map_to_v into the block (row_start..row_start+stalk_dim, v*stalk_dim..(v+1)*stalk_dim)
        for i in 0..stalk_dim {
            for j in 0..stalk_dim {
                delta[(row_start + i, v * stalk_dim + j)] += map_to_v[(i, j)];
                delta[(row_start + i, u * stalk_dim + j)] -= map_to_u[(i, j)];
            }
        }
    }

    delta
}

// ── Hodge Laplacian ─────────────────────────────────────────────────

/// The Hodge Laplacian L₁ = δ*δ + δδ* acting on 1-cochains.
///
/// In practice we compute:
/// - L₀ = δ*δ (acting on 0-cochains) — this is the usual sheaf Laplacian.
/// - L₁ = δδ* + δ*δ (acting on 1-cochains).
#[derive(Debug, Clone)]
pub struct HodgeLaplacian {
    /// Coboundary δ: C⁰ → C¹.
    pub delta: DMatrix<f64>,
    /// Transpose / adjoint δ*: C¹ → C⁰.
    pub delta_star: DMatrix<f64>,
    /// L₀ = δ*δ (sheaf Laplacian on 0-cochains).
    pub l0: DMatrix<f64>,
    /// L₁ = δδ* + δ*δ₁ (higher Laplacian, if triangle maps provided).
    pub l1: DMatrix<f64>,
}

impl HodgeLaplacian {
    /// Build Hodge Laplacian from the coboundary matrix alone (no triangle data).
    ///
    /// L₀ = δᵀ δ   (dim = n_vertices × stalk_dim)
    /// L₁ = δ δᵀ    (dim = n_edges × stalk_dim)
    pub fn from_coboundary(delta: DMatrix<f64>) -> Self {
        let delta_star = delta.transpose();
        let l0 = &delta_star * &delta;
        let l1 = &delta * &delta_star;
        Self {
            delta,
            delta_star,
            l0,
            l1,
        }
    }

    /// Eigenvalues of L₀ (sheaf Laplacian on 0-cochains).
    pub fn l0_eigenvalues(&self) -> Vec<f64> {
        let eig = SymmetricEigen::new(self.l0.clone());
        let mut vals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        vals
    }

    /// Eigenvalues of L₁ (higher Laplacian on 1-cochains).
    pub fn l1_eigenvalues(&self) -> Vec<f64> {
        let eig = SymmetricEigen::new(self.l1.clone());
        let mut vals: Vec<f64> = eig.eigenvalues.iter().copied().collect();
        vals.sort_by(|a, b| a.partial_cmp(b).unwrap());
        vals
    }

    /// Dimension of the 0-cochain space.
    pub fn dim0(&self) -> usize {
        self.l0.nrows()
    }

    /// Dimension of the 1-cochain space.
    pub fn dim1(&self) -> usize {
        self.l1.nrows()
    }

    /// Compute the kernel dimension of L₀ = dim H⁰(F).
    ///
    /// Counts eigenvalues below `tol` as zero.
    pub fn h0_dimension(&self, tol: f64) -> usize {
        self.l0_eigenvalues()
            .iter()
            .filter(|&&v| v.abs() < tol)
            .count()
    }

    /// Compute the kernel dimension of L₁ = dim H¹(F).
    pub fn h1_dimension(&self, tol: f64) -> usize {
        self.l1_eigenvalues()
            .iter()
            .filter(|&&v| v.abs() < tol)
            .count()
    }
}

// ── Harmonic sections ────────────────────────────────────────────────

/// Find harmonic 0-sections: vectors s in C⁰ such that L₀ s = 0.
///
/// Returns an orthonormal basis for the kernel of L₀ (eigenvectors with
/// eigenvalue < tol).
pub fn harmonic_0_sections(hodge: &HodgeLaplacian, tol: f64) -> Vec<DVector<f64>> {
    let eig = SymmetricEigen::new(hodge.l0.clone());
    let mut basis = Vec::new();
    for (i, &lambda) in eig.eigenvalues.iter().enumerate() {
        if lambda.abs() < tol {
            let col = eig.eigenvectors.column(i);
            basis.push(DVector::from(col.iter().copied().collect::<Vec<f64>>()));
        }
    }
    basis
}

/// Find harmonic 1-sections: vectors s in C¹ such that L₁ s = 0.
pub fn harmonic_1_sections(hodge: &HodgeLaplacian, tol: f64) -> Vec<DVector<f64>> {
    let eig = SymmetricEigen::new(hodge.l1.clone());
    let mut basis = Vec::new();
    for (i, &lambda) in eig.eigenvalues.iter().enumerate() {
        if lambda.abs() < tol {
            let col = eig.eigenvectors.column(i);
            basis.push(DVector::from(col.iter().copied().collect::<Vec<f64>>()));
        }
    }
    basis
}

/// Verify the Hodge decomposition: dim(ker L₀) = dim H⁰(F)
/// and for a connected complex with trivial stalk maps, H⁰ ≅ R^k.
pub fn verify_hodge_decomposition(
    hodge: &HodgeLaplacian,
    expected_h0: usize,
    tol: f64,
) -> bool {
    hodge.h0_dimension(tol) == expected_h0
}

// ── Discrete vector bundles and connection Laplacian ────────────────

/// The **connection Laplacian** for a graph equipped with orthogonal
/// restriction maps (a discrete vector bundle with connection).
///
/// L_conn(v,w) = -F_{v←w}  if v~w
/// L_conn(v,v) = Σ_{v~w} I
///
/// This is always positive semi-definite; its kernel consists of
/// parallel sections.
pub fn connection_laplacian(
    n_vertices: usize,
    edges: &[(usize, usize)],
    stalk_dim: usize,
    restriction_maps: &[(DMatrix<f64>, DMatrix<f64>)],
) -> DMatrix<f64> {
    let dim = n_vertices * stalk_dim;
    let mut l = DMatrix::zeros(dim, dim);

    for (e_idx, &(u, v)) in edges.iter().enumerate() {
        let (map_to_u, map_to_v) = &restriction_maps[e_idx];

        // Off-diagonal block: L(u,v) = -map_to_u, L(v,u) = -map_to_v^T
        for i in 0..stalk_dim {
            for j in 0..stalk_dim {
                l[(u * stalk_dim + i, v * stalk_dim + j)] -= map_to_u[(i, j)];
                l[(v * stalk_dim + i, u * stalk_dim + j)] -= map_to_v[(i, j)];

                // Diagonal blocks accumulate identity
                l[(u * stalk_dim + i, u * stalk_dim + j)] += map_to_u[(i, j)];
                l[(v * stalk_dim + i, v * stalk_dim + j)] += map_to_v[(i, j)];
            }
        }
    }

    l
}

/// Compute the **sheaf Betti numbers** β₀ and β₁ from the Hodge Laplacian.
pub fn sheaf_betti_numbers(hodge: &HodgeLaplacian, tol: f64) -> (usize, usize) {
    (hodge.h0_dimension(tol), hodge.h1_dimension(tol))
}

// ── Spectral gap ─────────────────────────────────────────────────────

/// Compute the spectral gap of L₀ (smallest non-zero eigenvalue).
///
/// Returns `None` if all eigenvalues are zero (fully disconnected).
pub fn spectral_gap(hodge: &HodgeLaplacian, tol: f64) -> Option<f64> {
    let eigs = hodge.l0_eigenvalues();
    eigs.iter()
        .find(|&&v| v > tol)
        .copied()
}

/// Estimate the sheaf's connectivity from the spectral gap.
///
/// A larger spectral gap means stronger "agreement" between neighboring stalks.
pub fn sheaf_connectivity(hodge: &HodgeLaplacian, tol: f64) -> f64 {
    spectral_gap(hodge, tol).unwrap_or(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    /// Helper: build a triangle graph (3 vertices, 3 edges) with identity maps.
    fn triangle_hodge() -> HodgeLaplacian {
        let n = 3;
        let edges = vec![(0, 1), (1, 2), (0, 2)];
        let id = DMatrix::identity(1, 1);
        let maps = vec![(id.clone(), id.clone()); 3];
        let delta = coboundary_matrix(n, &edges, 1, &maps);
        HodgeLaplacian::from_coboundary(delta)
    }

    #[test]
    fn test_coboundary_matrix_shape() {
        let n = 4;
        let edges = vec![(0, 1), (1, 2), (2, 3)];
        let id = DMatrix::identity(2, 2);
        let maps = vec![(id.clone(), id.clone()); 3];
        let delta = coboundary_matrix(n, &edges, 2, &maps);
        assert_eq!(delta.nrows(), 6);  // 3 edges × stalk_dim 2
        assert_eq!(delta.ncols(), 8);  // 4 vertices × stalk_dim 2
    }

    #[test]
    fn test_coboundary_constant_section() {
        // Constant section s = [1, 1, 1] on a triangle with identity maps
        let n = 3;
        let edges = vec![(0, 1), (1, 2), (0, 2)];
        let id = DMatrix::identity(1, 1);
        let maps = vec![(id.clone(), id.clone()); 3];
        let delta = coboundary_matrix(n, &edges, 1, &maps);

        let s = DVector::from_vec(vec![1.0, 1.0, 1.0]);
        let ds = &delta * &s;
        // δ(constant) should be zero
        for i in 0..ds.len() {
            assert!(ds[i].abs() < 1e-10, "constant section should be in ker δ");
        }
    }

    #[test]
    fn test_hodge_l0_triangle() {
        let hodge = triangle_hodge();
        // L₀ for 3-vertex triangle with identity maps = graph Laplacian
        let eigs = hodge.l0_eigenvalues();
        assert_eq!(eigs.len(), 3);
        // Eigenvalues of triangle graph Laplacian: 0, 3, 3
        assert!(eigs[0].abs() < 1e-10, "smallest eigenvalue should be ~0");
        assert!((eigs[2] - 3.0).abs() < 1e-10, "largest eigenvalue should be ~3");
    }

    #[test]
    fn test_hodge_l1_triangle() {
        let hodge = triangle_hodge();
        let eigs = hodge.l1_eigenvalues();
        assert_eq!(eigs.len(), 3);
        // L₁ = δδᵀ for triangle: eigenvalues 0, 3, 3
        assert!(eigs[0].abs() < 1e-10);
    }

    #[test]
    fn test_h0_dimension_connected() {
        let hodge = triangle_hodge();
        // Connected graph with identity maps → H⁰ ≅ R
        assert_eq!(hodge.h0_dimension(1e-8), 1);
    }

    #[test]
    fn test_h1_dimension_triangle() {
        let hodge = triangle_hodge();
        // L₁ = δδᵀ (up-Laplacian only, no triangle coboundary).
        // For triangle graph: rank(δ) = 2, so ker(δδᵀ) has dim 1.
        // This kernel contains the cycle class not killed without the
        // second coboundary map δ₁.
        assert_eq!(hodge.h1_dimension(1e-8), 1);
    }

    #[test]
    fn test_harmonic_0_sections() {
        let hodge = triangle_hodge();
        let harm = harmonic_0_sections(&hodge, 1e-8);
        assert_eq!(harm.len(), 1, "constant sheaf on connected complex has 1 harmonic section");
        // Should be proportional to [1,1,1]
        let s = &harm[0];
        let norm = s.iter().map(|x| x * x).sum::<f64>().sqrt();
        assert!((norm - 1.0).abs() < 1e-10, "should be normalised");
    }

    #[test]
    fn test_harmonic_1_sections() {
        let hodge = triangle_hodge();
        let harm = harmonic_1_sections(&hodge, 1e-8);
        // Up-Laplacian only: kernel has dim 1 (cycle class)
        assert_eq!(harm.len(), 1);
    }

    #[test]
    fn test_spectral_gap() {
        let hodge = triangle_hodge();
        let gap = spectral_gap(&hodge, 1e-8).unwrap();
        assert!((gap - 3.0).abs() < 1e-10, "spectral gap of triangle graph should be 3");
    }

    #[test]
    fn test_connection_laplacian_psd() {
        let n = 3;
        let edges = vec![(0, 1), (1, 2), (0, 2)];
        let id = DMatrix::identity(1, 1);
        let maps = vec![(id.clone(), id.clone()); 3];
        let l_conn = connection_laplacian(n, &edges, 1, &maps);

        // Check PSD: all eigenvalues ≥ 0
        let eig = SymmetricEigen::new(l_conn);
        for &v in eig.eigenvalues.iter() {
            assert!(v >= -1e-10, "connection Laplacian should be PSD, got eigenvalue {v}");
        }
    }

    #[test]
    fn test_sheaf_betti_numbers() {
        let hodge = triangle_hodge();
        let (b0, b1) = sheaf_betti_numbers(&hodge, 1e-8);
        assert_eq!(b0, 1);
        assert_eq!(b1, 1); // up-Laplacian only; full H¹=0 requires triangle maps
    }

    #[test]
    fn test_two_component_graph() {
        let n = 4;
        let edges = vec![(0, 1), (2, 3)]; // two disconnected edges
        let id = DMatrix::identity(1, 1);
        let maps = vec![(id.clone(), id.clone()); 2];
        let delta = coboundary_matrix(n, &edges, 1, &maps);
        let hodge = HodgeLaplacian::from_coboundary(delta);

        // Two connected components → H⁰ ≅ R²
        assert_eq!(hodge.h0_dimension(1e-8), 2);
    }

    #[test]
    fn test_weighted_restriction_maps() {
        let n = 3;
        let edges = vec![(0, 1), (1, 2), (0, 2)];
        // Non-trivial restriction maps (scaling)
        let w = DMatrix::from_vec(1, 1, vec![2.0]);
        let maps = vec![
            (w.clone(), w.clone()),
            (w.clone(), w.clone()),
            (w.clone(), w.clone()),
        ];
        let delta = coboundary_matrix(n, &edges, 1, &maps);
        let hodge = HodgeLaplacian::from_coboundary(delta);

        // H⁰ should still be 1 (connected with consistent maps)
        assert_eq!(hodge.h0_dimension(1e-8), 1);

        // Spectral gap should be 4× the unweighted (since maps scale by 2)
        let gap = spectral_gap(&hodge, 1e-8).unwrap();
        assert!((gap - 12.0).abs() < 1e-10, "spectral gap should scale with map weights");
    }
}
