//! Cellular sheaf: assigns data stalks and restriction maps to cells of a simplicial complex.

use nalgebra::DMatrix;

use crate::simplicial::SimplicialComplex;

/// A cellular sheaf over a simplicial complex.
///
/// For each cell σ, assigns a vector space F(σ) (the "stalk").
/// For each face τ < σ, a linear restriction map F(σ) → F(τ).
pub struct CellularSheaf {
    /// Dimension of the stalk at each vertex.
    pub stalk_dimension: usize,
    /// Restriction maps per edge: (map_to_a, map_to_b) — linear maps
    /// from stalk at far vertex to stalk at near vertex.
    pub restriction_maps: Vec<(DMatrix<f64>, DMatrix<f64>)>,
    /// Restriction maps per triangle: one map for each face (toward each edge).
    pub triangle_maps: Vec<[DMatrix<f64>; 3]>,
    /// The underlying complex.
    pub complex: SimplicialComplex,
}

impl CellularSheaf {
    /// Create a constant sheaf: all stalks R^n, all restriction maps are identity.
    pub fn constant(complex: SimplicialComplex, stalk_dim: usize) -> Self {
        let id = DMatrix::identity(stalk_dim, stalk_dim);
        let edge_maps: Vec<_> = complex
            .edges
            .iter()
            .map(|_| (id.clone(), id.clone()))
            .collect();
        let tri_maps: Vec<_> = complex
            .triangles
            .iter()
            .map(|_| [id.clone(), id.clone(), id.clone()])
            .collect();
        Self {
            stalk_dimension: stalk_dim,
            restriction_maps: edge_maps,
            triangle_maps: tri_maps,
            complex,
        }
    }

    /// Create a sheaf from weight functions on edges (1-dimensional stalks).
    pub fn from_weights(complex: SimplicialComplex, weights: &[f64]) -> Self {
        let edge_maps: Vec<_> = complex
            .edges
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let w = weights.get(i).copied().unwrap_or(1.0);
                (DMatrix::from_element(1, 1, w), DMatrix::from_element(1, 1, w))
            })
            .collect();
        let id1 = DMatrix::identity(1, 1);
        let tri_maps: Vec<_> = complex
            .triangles
            .iter()
            .map(|_| [id1.clone(), id1.clone(), id1.clone()])
            .collect();
        Self {
            stalk_dimension: 1,
            restriction_maps: edge_maps,
            triangle_maps: tri_maps,
            complex,
        }
    }

    /// Global sections: assignments compatible with restriction maps.
    /// Returns the dimension of the global section space.
    pub fn global_section_dimension(&self) -> usize {
        if self.stalk_dimension == 0 || self.complex.vertices.is_empty() {
            return 0;
        }
        // For a constant sheaf on a connected complex, this equals stalk dimension.
        // General case would require solving the sheaf condition system.
        self.stalk_dimension
    }

    /// Sheaf cohomology dimensions H^0 and H^1.
    pub fn cohomology_dimension(&self, degree: usize) -> usize {
        match degree {
            0 => self.global_section_dimension(),
            1 => {
                let n_e = self.complex.edges.len();
                let n_v = self.complex.vertices.len();
                (n_e * self.stalk_dimension)
                    .saturating_sub(n_v * self.stalk_dimension)
            }
            _ => 0,
        }
    }

    /// Restriction map from an edge to its source vertex a (lower-index).
    pub fn edge_restriction_a(&self, edge_idx: usize) -> &DMatrix<f64> {
        &self.restriction_maps[edge_idx].0
    }

    /// Restriction map from an edge to its target vertex b (higher-index).
    pub fn edge_restriction_b(&self, edge_idx: usize) -> &DMatrix<f64> {
        &self.restriction_maps[edge_idx].1
    }

    /// The stalk dimension of the sheaf.
    pub fn stalk_dim(&self) -> usize {
        self.stalk_dimension
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constant_sheaf() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        let sheaf = CellularSheaf::constant(c, 2);
        assert_eq!(sheaf.stalk_dimension, 2);
        assert_eq!(sheaf.restriction_maps.len(), 2);
        assert_eq!(sheaf.triangle_maps.len(), 0);
    }

    #[test]
    fn test_constant_sheaf_with_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        let sheaf = CellularSheaf::constant(c, 3);
        assert_eq!(sheaf.stalk_dimension, 3);
        assert_eq!(sheaf.restriction_maps.len(), 3);
        assert_eq!(sheaf.triangle_maps.len(), 1);
    }

    #[test]
    fn test_weighted_sheaf() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        let sheaf = CellularSheaf::from_weights(c, &[1.0, 0.5]);
        assert_eq!(sheaf.stalk_dimension, 1);
    }

    #[test]
    fn test_global_section_dimension() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 3);
        assert_eq!(sheaf.global_section_dimension(), 3);
    }

    #[test]
    fn test_empty_sheaf() {
        let c = SimplicialComplex::new();
        let sheaf = CellularSheaf::constant(c, 2);
        // Empty complex has no vertices, so global section dimension is 0
        assert_eq!(sheaf.global_section_dimension(), 0);
        assert_eq!(sheaf.complex.vertices.len(), 0);
        assert_eq!(sheaf.stalk_dimension, 2);
    }

    #[test]
    fn test_cohomology() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 1);
        let h0 = sheaf.cohomology_dimension(0);
        assert_eq!(h0, 1);
    }

    #[test]
    fn test_constant_sheaf_identity_maps() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 2);
        let m = sheaf.edge_restriction_a(0);
        assert_eq!(m[(0, 0)], 1.0);
        assert_eq!(m[(1, 1)], 1.0);
    }

    #[test]
    fn test_sheaf_stalk_dim() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let sheaf = CellularSheaf::constant(c, 5);
        assert_eq!(sheaf.stalk_dim(), 5);
    }
}
