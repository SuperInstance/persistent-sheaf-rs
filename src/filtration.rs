//! Filtration: a sequence of nested simplicial complexes.

use crate::persistence::PersistenceDiagram;
use crate::simplicial::SimplicialComplex;

/// A filtration: an increasing sequence of simplicial complexes indexed by a scale parameter.
pub struct Filtration {
    /// The complexes at each threshold.
    pub complexes: Vec<(f64, SimplicialComplex)>,
}

impl Filtration {
    /// Create an empty filtration.
    pub fn new() -> Self {
        Self { complexes: vec![] }
    }

    /// Add a complex at a given threshold.
    pub fn add(&mut self, threshold: f64, complex: SimplicialComplex) {
        self.complexes.push((threshold, complex));
    }

    /// Build a Vietoris-Rips filtration from a point cloud (N × D).
    pub fn from_point_cloud(points: &[Vec<f64>], num_steps: usize) -> Self {
        if points.is_empty() {
            return Self::new();
        }
        let n = points.len();
        let mut max_dist = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                let d = crate::simplicial::euclidean_dist(&points[i], &points[j]);
                if d > max_dist {
                    max_dist = d;
                }
            }
        }
        let mut filt = Self::new();
        for step in 0..=num_steps {
            let epsilon = if num_steps > 0 {
                max_dist * step as f64 / num_steps as f64
            } else {
                0.0
            };
            let complex = SimplicialComplex::vietoris_rips(points, epsilon);
            filt.add(epsilon, complex);
        }
        filt
    }

    /// Build a filtration from a distance matrix by varying epsilon.
    pub fn from_distance_matrix(distances: &[Vec<f64>], num_steps: usize) -> Self {
        if distances.is_empty() {
            return Self::new();
        }
        let n = distances.len();
        let mut max_dist = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                if distances[i][j] > max_dist {
                    max_dist = distances[i][j];
                }
            }
        }
        let mut filt = Self::new();
        for step in 0..=num_steps {
            let epsilon = if num_steps > 0 {
                max_dist * step as f64 / num_steps as f64
            } else {
                0.0
            };
            let complex = SimplicialComplex::vietoris_rips_from_distances(distances, epsilon);
            filt.add(epsilon, complex);
        }
        filt
    }

    /// Compute persistent homology from this filtration.
    /// Returns a PersistenceDiagram tracking birth/death of H0 and H1 features.
    pub fn compute_persistence(&self) -> PersistenceDiagram {
        let mut diagram = PersistenceDiagram::new();

        // Track H0 (connected components) and H1 (cycles) across the filtration
        let mut prev_complex: Option<&SimplicialComplex> = None;

        for (threshold, complex) in &self.complexes {
            let betti = complex.betti_numbers();
            let beta0 = betti.first().copied().unwrap_or(0);
            let beta1 = betti.get(1).copied().unwrap_or(0);

            if let Some(prev) = prev_complex {
                let prev_betti = prev.betti_numbers();
                let prev_beta0 = prev_betti.first().copied().unwrap_or(0);

                // If components merged, record death at this threshold
                if prev_beta0 > beta0 {
                    let n_merged = prev_beta0 - beta0;
                    for _ in 0..n_merged {
                        diagram.add(0.0, *threshold, 0);
                    }
                }

                // Check for new H1 features (cycles appearing)
                let prev_beta1 = prev_betti.get(1).copied().unwrap_or(0);
                if beta1 > prev_beta1 {
                    let n_new = beta1 - prev_beta1;
                    for _ in 0..n_new {
                        diagram.add(*threshold, f64::INFINITY, 1);
                    }
                }

                // Check for H1 features dying (being filled by triangles)
                if beta1 < prev_beta1 {
                    let n_dying = prev_beta1 - beta1;
                    // Feature died when a triangle filled it
                    for _ in 0..n_dying {
                        diagram.add(0.0, *threshold, 1);
                    }
                }
            } else {
                // At the first step, all vertices are births for H0
                for _ in 0..beta0 {
                    diagram.add(*threshold, f64::INFINITY, 0);
                }
            }

            prev_complex = Some(complex);
        }

        diagram
    }

    /// Threshold values.
    pub fn thresholds(&self) -> Vec<f64> {
        self.complexes.iter().map(|(t, _)| *t).collect()
    }

    /// Number of steps in the filtration.
    pub fn len(&self) -> usize {
        self.complexes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.complexes.is_empty()
    }

    /// Get the complex at a specific step.
    pub fn get(&self, index: usize) -> Option<&(f64, SimplicialComplex)> {
        self.complexes.get(index)
    }

    /// Get the maximum threshold.
    pub fn max_threshold(&self) -> f64 {
        self.complexes
            .last()
            .map(|(t, _)| *t)
            .unwrap_or(0.0)
    }

    /// Make a constant sheaf on every complex in the filtration.
    pub fn make_sheaf_filtration(&self, stalk_dim: usize) -> Vec<crate::sheaf::CellularSheaf> {
        self.complexes
            .iter()
            .map(|(_, c)| crate::sheaf::CellularSheaf::constant(c.clone(), stalk_dim))
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_filtration() {
        let f = Filtration::new();
        assert!(f.is_empty());
        assert_eq!(f.len(), 0);
    }

    #[test]
    fn test_from_distance_matrix() {
        let distances = vec![
            vec![0.0, 1.0, 3.0],
            vec![1.0, 0.0, 1.0],
            vec![3.0, 1.0, 0.0],
        ];
        let f = Filtration::from_distance_matrix(&distances, 5);
        assert!(f.len() > 0);
    }

    #[test]
    fn test_from_point_cloud() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![2.0, 0.0],
        ];
        let f = Filtration::from_point_cloud(&points, 5);
        assert!(f.len() > 0);
    }

    #[test]
    fn test_compute_persistence() {
        let distances = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ];
        let f = Filtration::from_distance_matrix(&distances, 10);
        let diagram = f.compute_persistence();
        assert!(diagram.len() > 0);
    }

    #[test]
    fn test_filtration_grows() {
        let distances = vec![
            vec![0.0, 1.0, 2.0],
            vec![1.0, 0.0, 1.0],
            vec![2.0, 1.0, 0.0],
        ];
        let f = Filtration::from_distance_matrix(&distances, 5);
        for w in f.complexes.windows(2) {
            assert!(w[0].0 <= w[1].0);
        }
    }

    #[test]
    fn test_thresholds() {
        let distances = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let f = Filtration::from_distance_matrix(&distances, 4);
        let thresholds = f.thresholds();
        assert_eq!(thresholds.len(), 5);
        assert!((thresholds[0] - 0.0).abs() < 1e-10);
        assert!((thresholds[4] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_get_complex() {
        let points = vec![vec![0.0], vec![1.0]];
        let f = Filtration::from_point_cloud(&points, 2);
        let entry = f.get(0);
        assert!(entry.is_some());
        assert_eq!(entry.unwrap().1.num_simplices(0), 2);
    }

    #[test]
    fn test_max_threshold() {
        let distances = vec![vec![0.0, 2.0], vec![2.0, 0.0]];
        let f = Filtration::from_distance_matrix(&distances, 10);
        assert!((f.max_threshold() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn test_make_sheaf_filtration() {
        let points = vec![vec![0.0], vec![1.0], vec![2.0]];
        let f = Filtration::from_point_cloud(&points, 3);
        let sheaves = f.make_sheaf_filtration(2);
        assert_eq!(sheaves.len(), f.len());
        for s in &sheaves {
            assert_eq!(s.stalk_dimension, 2);
        }
    }

    #[test]
    fn test_empty_point_cloud() {
        let points: Vec<Vec<f64>> = vec![];
        let f = Filtration::from_point_cloud(&points, 5);
        assert!(f.is_empty());
    }

    #[test]
    fn test_empty_distance_matrix() {
        let d: Vec<Vec<f64>> = vec![];
        let f = Filtration::from_distance_matrix(&d, 5);
        assert!(f.is_empty());
    }

    #[test]
    fn test_compute_persistence_single_point() {
        let points = vec![vec![0.0]];
        let f = Filtration::from_point_cloud(&points, 5);
        let diagram = f.compute_persistence();
        assert!(diagram.len() > 0);
    }

    #[test]
    fn test_persistence_two_points() {
        let distances = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let f = Filtration::from_distance_matrix(&distances, 10);
        let diagram = f.compute_persistence();
        // Two points start separate, merge → should have deaths
        assert!(diagram.len() >= 1);
    }
}
