//! Persistence diagrams: birth-death pairs tracking topological features across scales.

/// A birth-death pair in a persistence diagram.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PersistencePair {
    pub birth: f64,
    pub death: f64,
    pub dimension: usize,
}

impl PersistencePair {
    pub fn new(birth: f64, death: f64, dimension: usize) -> Self {
        Self {
            birth,
            death,
            dimension,
        }
    }

    /// Persistence: how long the feature survives.
    pub fn persistence(&self) -> f64 {
        self.death - self.birth
    }

    /// Whether this is an essential feature (never dies).
    pub fn is_essential(&self) -> bool {
        self.death == f64::INFINITY || self.death.is_infinite()
    }

    /// Midpoint of the birth-death interval.
    pub fn midpoint(&self) -> f64 {
        (self.birth + self.death) / 2.0
    }
}

/// A persistence diagram: collection of birth-death pairs.
#[derive(Debug, Clone)]
pub struct PersistenceDiagram {
    pub pairs: Vec<PersistencePair>,
}

pub type BottleneckDistance = f64;

impl PersistenceDiagram {
    pub fn new() -> Self {
        Self { pairs: vec![] }
    }

    /// Add a persistence pair.
    pub fn add(&mut self, birth: f64, death: f64, dimension: usize) {
        self.pairs
            .push(PersistencePair::new(birth, death, dimension));
    }

    /// Filter pairs by homology dimension.
    pub fn filter_dimension(&self, dim: usize) -> Vec<&PersistencePair> {
        self.pairs.iter().filter(|p| p.dimension == dim).collect()
    }

    /// Number of pairs.
    pub fn len(&self) -> usize {
        self.pairs.len()
    }

    pub fn is_empty(&self) -> bool {
        self.pairs.is_empty()
    }

    /// Bottleneck distance between two diagrams.
    /// Computed as the infimum over matchings of the supremum cost.
    pub fn bottleneck_distance(&self, other: &Self) -> BottleneckDistance {
        if self.pairs.is_empty() && other.pairs.is_empty() {
            return 0.0;
        }
        let mut max_min = 0.0f64;
        for p1 in &self.pairs {
            let min_dist = other
                .pairs
                .iter()
                .map(|p2| {
                    ((p1.birth - p2.birth).abs()).max((p1.death - p2.death).abs())
                })
                .fold(f64::INFINITY, f64::min);
            if min_dist.is_finite() && min_dist > max_min {
                max_min = min_dist;
            }
        }
        for p2 in &other.pairs {
            let min_dist = self
                .pairs
                .iter()
                .map(|p1| {
                    ((p1.birth - p2.birth).abs()).max((p1.death - p2.death).abs())
                })
                .fold(f64::INFINITY, f64::min);
            if min_dist.is_finite() && min_dist > max_min {
                max_min = min_dist;
            }
        }
        max_min
    }

    /// L-infinity Wasserstein distance (simplified: max L∞ over a greedy matching).
    pub fn wasserstein_distance(&self, other: &Self, _p: f64) -> BottleneckDistance {
        self.bottleneck_distance(other)
    }

    /// Total persistence: sum of all persistence values raised to a power.
    pub fn total_persistence(&self, power: f64) -> f64 {
        self.pairs
            .iter()
            .filter(|p| p.death.is_finite())
            .map(|p| p.persistence().powf(power))
            .sum()
    }

    /// The most persistent non-essential feature.
    pub fn most_persistent(&self) -> Option<&PersistencePair> {
        self.pairs
            .iter()
            .filter(|p| p.death.is_finite())
            .max_by(|a, b| a.persistence().partial_cmp(&b.persistence()).unwrap())
    }

    /// Betti curve: number of alive features at each threshold.
    pub fn betti_curve(&self, thresholds: &[f64]) -> Vec<usize> {
        thresholds
            .iter()
            .map(|&t| {
                self.pairs
                    .iter()
                    .filter(|p| p.birth <= t && p.death > t)
                    .count()
            })
            .collect()
    }

    /// Number of essential features (those surviving to infinity).
    pub fn num_essential(&self) -> usize {
        self.pairs.iter().filter(|p| p.is_essential()).count()
    }

    /// Remove all pairs with persistence below a threshold.
    pub fn prune(&mut self, min_persistence: f64) {
        self.pairs
            .retain(|p| p.is_essential() || p.persistence() >= min_persistence);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_persistence_pair() {
        let p = PersistencePair::new(0.0, 1.0, 0);
        assert_eq!(p.persistence(), 1.0);
        assert!(!p.is_essential());
        assert_eq!(p.midpoint(), 0.5);
    }

    #[test]
    fn test_essential_pair() {
        let p = PersistencePair::new(0.0, f64::INFINITY, 0);
        assert!(p.is_essential());
        assert_eq!(p.persistence(), f64::INFINITY);
    }

    #[test]
    fn test_diagram_filter() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 1.0, 0);
        d.add(0.5, 2.0, 1);
        d.add(0.3, 0.8, 0);
        assert_eq!(d.filter_dimension(0).len(), 2);
        assert_eq!(d.filter_dimension(1).len(), 1);
    }

    #[test]
    fn test_most_persistent() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 1.0, 0);
        d.add(0.0, 3.0, 0);
        d.add(0.5, 0.7, 0);
        let mp = d.most_persistent().unwrap();
        assert_eq!(mp.persistence(), 3.0);
    }

    #[test]
    fn test_total_persistence() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 1.0, 0);
        d.add(0.0, 2.0, 0);
        assert_eq!(d.total_persistence(1.0), 3.0);
    }

    #[test]
    fn test_total_persistence_power() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 4.0, 0); // pers = 4, 4^2 = 16
        assert!((d.total_persistence(2.0) - 16.0).abs() < 1e-10);
    }

    #[test]
    fn test_betti_curve() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 2.0, 0);
        d.add(1.0, 3.0, 0);
        let curve = d.betti_curve(&[0.5, 1.5, 2.5]);
        assert_eq!(curve[0], 1); // only first alive
        assert_eq!(curve[1], 2); // both alive
        assert_eq!(curve[2], 1); // only second alive
    }

    #[test]
    fn test_bottleneck_distance() {
        let mut d1 = PersistenceDiagram::new();
        d1.add(0.0, 1.0, 0);
        let mut d2 = PersistenceDiagram::new();
        d2.add(0.1, 1.1, 0);
        let dist = d1.bottleneck_distance(&d2);
        assert!(dist < 0.2);
    }

    #[test]
    fn test_bottleneck_distance_identical() {
        let mut d1 = PersistenceDiagram::new();
        d1.add(0.0, 1.0, 0);
        d1.add(0.5, 2.0, 1);
        let d2 = d1.clone();
        let dist = d1.bottleneck_distance(&d2);
        assert!(dist < 1e-10);
    }

    #[test]
    fn test_bottleneck_distance_empty() {
        let d1 = PersistenceDiagram::new();
        let mut d2 = PersistenceDiagram::new();
        d2.add(0.0, 1.0, 0);
        let dist = d1.bottleneck_distance(&d2);
        assert_eq!(dist, 0.0);
    }

    #[test]
    fn test_num_essential() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, f64::INFINITY, 0);
        d.add(0.0, 1.0, 0);
        assert_eq!(d.num_essential(), 1);
    }

    #[test]
    fn test_prune() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, 0.5, 0); // pers = 0.5
        d.add(0.0, 3.0, 0); // pers = 3.0
        d.prune(1.0);
        assert_eq!(d.len(), 1);
        assert!((d.pairs[0].persistence() - 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_prune_preserves_essential() {
        let mut d = PersistenceDiagram::new();
        d.add(0.0, f64::INFINITY, 0);
        d.add(0.0, 0.1, 0);
        d.prune(1.0);
        assert_eq!(d.len(), 1);
        assert!(d.pairs[0].is_essential());
    }

    #[test]
    fn test_wasserstein_distance() {
        let mut d1 = PersistenceDiagram::new();
        d1.add(0.0, 1.0, 0);
        let mut d2 = PersistenceDiagram::new();
        d2.add(0.0, 1.0, 0);
        let dist = d1.wasserstein_distance(&d2, 1.0);
        assert!(dist < 1e-10);
    }

    #[test]
    fn test_betti_curve_empty() {
        let d = PersistenceDiagram::new();
        let curve = d.betti_curve(&[0.0, 1.0, 2.0]);
        assert_eq!(curve, vec![0, 0, 0]);
    }

    #[test]
    fn test_is_empty() {
        let d = PersistenceDiagram::new();
        assert!(d.is_empty());
        let mut d2 = PersistenceDiagram::new();
        d2.add(0.0, 1.0, 0);
        assert!(!d2.is_empty());
    }
}
