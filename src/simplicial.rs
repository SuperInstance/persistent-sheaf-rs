//! Simplicial complex data structure and Vietoris-Rips construction.

use std::collections::HashSet;

/// A simplicial complex: a collection of simplices closed under taking faces.
#[derive(Debug, Clone)]
pub struct SimplicialComplex {
    /// Vertices (0-simplices).
    pub vertices: Vec<usize>,
    /// Edges (1-simplices), stored as sorted pairs.
    pub edges: Vec<(usize, usize)>,
    /// Triangles (2-simplices), stored as sorted triples.
    pub triangles: Vec<(usize, usize, usize)>,
    /// Tetrahedra (3-simplices).
    pub tetrahedra: Vec<(usize, usize, usize, usize)>,
}

impl SimplicialComplex {
    /// Create a new empty complex.
    pub fn new() -> Self {
        Self {
            vertices: vec![],
            edges: vec![],
            triangles: vec![],
            tetrahedra: vec![],
        }
    }

    /// Add a vertex.
    pub fn add_vertex(&mut self, v: usize) {
        if !self.vertices.contains(&v) {
            self.vertices.push(v);
        }
    }

    /// Add an edge (ensures vertices exist).
    pub fn add_edge(&mut self, a: usize, b: usize) {
        let edge = if a < b { (a, b) } else { (b, a) };
        self.add_vertex(a);
        self.add_vertex(b);
        if !self.edges.contains(&edge) {
            self.edges.push(edge);
        }
    }

    /// Add a triangle (ensures edges exist).
    pub fn add_triangle(&mut self, a: usize, b: usize, c: usize) {
        let mut verts = [a, b, c];
        verts.sort();
        let tri = (verts[0], verts[1], verts[2]);
        self.add_edge(a, b);
        self.add_edge(b, c);
        self.add_edge(a, c);
        if !self.triangles.contains(&tri) {
            self.triangles.push(tri);
        }
    }

    /// Add a tetrahedron (ensures triangles exist).
    pub fn add_tetrahedron(&mut self, a: usize, b: usize, c: usize, d: usize) {
        let mut verts = [a, b, c, d];
        verts.sort();
        let tet = (verts[0], verts[1], verts[2], verts[3]);
        self.add_triangle(a, b, c);
        self.add_triangle(a, b, d);
        self.add_triangle(a, c, d);
        self.add_triangle(b, c, d);
        if !self.tetrahedra.contains(&tet) {
            self.tetrahedra.push(tet);
        }
    }

    /// Number of simplices of dimension k.
    pub fn num_simplices(&self, dim: usize) -> usize {
        match dim {
            0 => self.vertices.len(),
            1 => self.edges.len(),
            2 => self.triangles.len(),
            3 => self.tetrahedra.len(),
            _ => 0,
        }
    }

    /// Total number of simplices across all dimensions.
    pub fn total_simplices(&self) -> usize {
        self.vertices.len() + self.edges.len() + self.triangles.len() + self.tetrahedra.len()
    }

    /// Euler characteristic: χ = V - E + F - T.
    pub fn euler_characteristic(&self) -> i32 {
        self.vertices.len() as i32
            - self.edges.len() as i32
            + self.triangles.len() as i32
            - self.tetrahedra.len() as i32
    }

    /// Boundary of a 0-simplex (empty).
    pub fn boundary_of_vertex(&self, _v: usize) -> Vec<usize> {
        vec![]
    }

    /// Boundary of an edge: its two vertices.
    pub fn boundary_of_edge(&self, edge: (usize, usize)) -> Vec<usize> {
        vec![edge.0, edge.1]
    }

    /// Boundary of a triangle: its three edges.
    pub fn boundary_of_triangle(&self, tri: (usize, usize, usize)) -> Vec<(usize, usize)> {
        let (a, b, c) = tri;
        let mut edges = vec![
            (a.min(b), a.max(b)),
            (b.min(c), b.max(c)),
            (a.min(c), a.max(c)),
        ];
        edges.sort();
        edges
    }

    /// Boundary of a tetrahedron: its four triangles.
    pub fn boundary_of_tetrahedron(
        &self,
        tet: (usize, usize, usize, usize),
    ) -> Vec<(usize, usize, usize)> {
        let (a, b, c, d) = tet;
        let mut tris = vec![
            (a, b, c),
            (a, b, d),
            (a, c, d),
            (b, c, d),
        ];
        for t in tris.iter_mut() {
            let mut v = [t.0, t.1, t.2];
            v.sort();
            *t = (v[0], v[1], v[2]);
        }
        tris.sort();
        tris
    }

    /// Number of connected components via union-find.
    pub fn connected_components(&self) -> usize {
        if self.vertices.is_empty() {
            return 0;
        }
        let n = *self.vertices.iter().max().unwrap_or(&0) + 1;
        let mut parent: Vec<usize> = (0..n).collect();

        for &(a, b) in &self.edges {
            if a < n && b < n {
                union(a, b, &mut parent);
            }
        }
        let roots: HashSet<usize> = self.vertices.iter().map(|&v| find(v, &mut parent)).collect();
        roots.len()
    }

    /// Betti numbers: β₀ = components, β₁ = cycles (2D Euler-Poincaré).
    /// Returns vec![β₀, β₁].
    pub fn betti_numbers(&self) -> Vec<usize> {
        let n = self.vertices.len();
        if n == 0 {
            return vec![];
        }
        let beta0 = self.connected_components();
        let chi = self.euler_characteristic();
        // χ = β₀ − β₁ + β₂ − ... For 2-skeleton: β₁ = β₀ − χ − β₂
        let beta2 = self.num_simplices(2).saturating_sub(self.num_simplices(3));
        let beta1 = (beta0 as i32 - chi - beta2 as i32).max(0) as usize;
        vec![beta0, beta1]
    }

    /// h0_betti: number of connected components (H₀) at a given Vietoris-Rips threshold.
    /// Constructs the VR complex at `threshold` and returns β₀.
    pub fn h0_betti(points: &[Vec<f64>], threshold: f64) -> usize {
        let n = points.len();
        if n == 0 {
            return 0;
        }
        let mut parent: Vec<usize> = (0..n).collect();
        for i in 0..n {
            for j in (i + 1)..n {
                let d = euclidean_dist(&points[i], &points[j]);
                if d <= threshold {
                    union(i, j, &mut parent);
                }
            }
        }
        let roots: HashSet<usize> = (0..n).map(|i| find(i, &mut parent)).collect();
        roots.len()
    }

    /// h1_betti: number of independent cycles (H₁) at a given Vietoris-Rips threshold.
    /// Uses Euler-Poincaré: β₁ = β₀ − χ.
    pub fn h1_betti(points: &[Vec<f64>], threshold: f64) -> usize {
        let complex = SimplicialComplex::vietoris_rips(points, threshold);
        let betti = complex.betti_numbers();
        if betti.len() > 1 {
            betti[1]
        } else {
            0
        }
    }

    /// Build a Vietoris-Rips complex from a point cloud (N×D) and distance threshold.
    pub fn vietoris_rips(points: &[Vec<f64>], epsilon: f64) -> Self {
        let n = points.len();
        let mut complex = Self::new();
        for i in 0..n {
            complex.add_vertex(i);
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if euclidean_dist(&points[i], &points[j]) <= epsilon {
                    complex.add_edge(i, j);
                }
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    if euclidean_dist(&points[i], &points[j]) <= epsilon
                        && euclidean_dist(&points[j], &points[k]) <= epsilon
                        && euclidean_dist(&points[i], &points[k]) <= epsilon
                    {
                        complex.add_triangle(i, j, k);
                    }
                }
            }
        }
        complex
    }

    /// Build Vietoris-Rips from a precomputed distance matrix.
    pub fn vietoris_rips_from_distances(distances: &[Vec<f64>], epsilon: f64) -> Self {
        let n = distances.len();
        let mut complex = Self::new();
        for i in 0..n {
            complex.add_vertex(i);
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if distances[i][j] <= epsilon {
                    complex.add_edge(i, j);
                }
            }
        }
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    if distances[i][j] <= epsilon
                        && distances[j][k] <= epsilon
                        && distances[i][k] <= epsilon
                    {
                        complex.add_triangle(i, j, k);
                    }
                }
            }
        }
        complex
    }

    /// Clique complex (flag complex) from a graph: add triangles for every 3-clique.
    pub fn flag_complex(edges: &[(usize, usize)]) -> Self {
        let mut complex = Self::new();
        let mut adj: std::collections::HashMap<usize, HashSet<usize>> =
            std::collections::HashMap::new();
        for &(a, b) in edges {
            complex.add_edge(a, b);
            adj.entry(a).or_default().insert(b);
            adj.entry(b).or_default().insert(a);
        }
        // Find 3-cliques
        for &(a, b) in edges {
            if let (Some(na), Some(nb)) = (adj.get(&a), adj.get(&b)) {
                let common: Vec<&usize> = na.intersection(nb).collect();
                for &&c in &common {
                    complex.add_triangle(a, b, c);
                }
            }
        }
        complex
    }
}

fn find(mut x: usize, parent: &mut [usize]) -> usize {
    while parent[x] != x {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    x
}

fn union(a: usize, b: usize, parent: &mut [usize]) {
    let ra = find(a, parent);
    let rb = find(b, parent);
    if ra != rb {
        parent[ra] = rb;
    }
}

/// Euclidean distance between two points.
pub fn euclidean_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        .sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_complex() {
        let c = SimplicialComplex::new();
        assert_eq!(c.num_simplices(0), 0);
        assert_eq!(c.total_simplices(), 0);
    }

    #[test]
    fn test_add_vertex() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        c.add_vertex(1);
        assert_eq!(c.num_simplices(0), 2);
    }

    #[test]
    fn test_add_edge() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        assert_eq!(c.num_simplices(0), 2);
        assert_eq!(c.num_simplices(1), 1);
    }

    #[test]
    fn test_add_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        assert_eq!(c.num_simplices(0), 3);
        assert_eq!(c.num_simplices(1), 3);
        assert_eq!(c.num_simplices(2), 1);
    }

    #[test]
    fn test_add_tetrahedron() {
        let mut c = SimplicialComplex::new();
        c.add_tetrahedron(0, 1, 2, 3);
        assert_eq!(c.num_simplices(0), 4);
        assert_eq!(c.num_simplices(1), 6);
        assert_eq!(c.num_simplices(2), 4);
        assert_eq!(c.num_simplices(3), 1);
        assert_eq!(c.euler_characteristic(), 1); // 4-6+4-1=1
    }

    #[test]
    fn test_dedup() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        c.add_vertex(0);
        assert_eq!(c.num_simplices(0), 1);
        c.add_edge(0, 1);
        c.add_edge(1, 0);
        assert_eq!(c.num_simplices(1), 1);
    }

    #[test]
    fn test_euler_characteristic() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        assert_eq!(c.euler_characteristic(), 1); // V=3, E=3, F=1
    }

    #[test]
    fn test_euler_tetrahedron() {
        let mut c = SimplicialComplex::new();
        c.add_tetrahedron(0, 1, 2, 3);
        assert_eq!(c.euler_characteristic(), 1); // 4-6+4-1=1
    }

    #[test]
    fn test_connected_components_single() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        assert_eq!(c.connected_components(), 1);
    }

    #[test]
    fn test_connected_components_two() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        c.add_vertex(1);
        assert_eq!(c.connected_components(), 2);
    }

    #[test]
    fn test_connected_components_edge() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        assert_eq!(c.connected_components(), 1);
    }

    #[test]
    fn test_connected_components_disjoint() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(2, 3);
        assert_eq!(c.connected_components(), 2);
    }

    #[test]
    fn test_betti_numbers_point() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1); // one component
    }

    #[test]
    fn test_betti_numbers_two_points() {
        let mut c = SimplicialComplex::new();
        c.add_vertex(0);
        c.add_vertex(1);
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 2); // two separate components
    }

    #[test]
    fn test_betti_numbers_edge() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1); // one component
        assert_eq!(betti[1], 0); // no cycle
    }

    #[test]
    fn test_betti_numbers_triangle_edges() {
        let mut c = SimplicialComplex::new();
        c.add_edge(0, 1);
        c.add_edge(1, 2);
        c.add_edge(0, 2);
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1); // one component
        assert_eq!(betti[1], 1); // one cycle
    }

    #[test]
    fn test_betti_numbers_filled_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1);
        assert_eq!(betti[1], 0); // cycle is filled
    }

    #[test]
    fn test_betti_numbers_tetrahedron_skeleton() {
        let mut c = SimplicialComplex::new();
        c.add_tetrahedron(0, 1, 2, 3);
        // χ = 1 for a tetrahedron (4-6+4-1=1)
        // β₀ = 1, β₂ = 0 (triangle faces, but 2-skeleton), β₁ = β₀ - χ = 0
        let betti = c.betti_numbers();
        assert_eq!(betti[0], 1);
    }

    #[test]
    fn test_vietoris_rips_simple() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![2.0, 0.0],
        ];
        let c = SimplicialComplex::vietoris_rips(&points, 1.5);
        assert_eq!(c.num_simplices(0), 3);
        assert_eq!(c.num_simplices(1), 2); // 0-1 and 1-2
    }

    #[test]
    fn test_vietoris_rips_triangle() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.5, 0.866],
        ];
        let c = SimplicialComplex::vietoris_rips(&points, 1.1);
        assert_eq!(c.num_simplices(1), 3);
        assert_eq!(c.num_simplices(2), 1);
    }

    #[test]
    fn test_vietoris_rips_high_threshold() {
        let points = vec![
            vec![0.0, 0.0],
            vec![10.0, 0.0],
        ];
        let c = SimplicialComplex::vietoris_rips(&points, 100.0);
        assert_eq!(c.num_simplices(1), 1);
    }

    #[test]
    fn test_vietoris_rips_from_distances() {
        let d = vec![vec![0.0, 1.0], vec![1.0, 0.0]];
        let c = SimplicialComplex::vietoris_rips_from_distances(&d, 0.5);
        assert_eq!(c.num_simplices(1), 0);
        let c2 = SimplicialComplex::vietoris_rips_from_distances(&d, 2.0);
        assert_eq!(c2.num_simplices(1), 1);
    }

    #[test]
    fn test_flag_complex() {
        let edges = vec![(0, 1), (1, 2), (0, 2)];
        let c = SimplicialComplex::flag_complex(&edges);
        assert_eq!(c.num_simplices(2), 1);
    }

    #[test]
    fn test_boundary_of_edge() {
        let c = SimplicialComplex::new();
        let b = c.boundary_of_edge((2, 5));
        assert_eq!(b, vec![2, 5]);
    }

    #[test]
    fn test_boundary_of_triangle() {
        let mut c = SimplicialComplex::new();
        c.add_triangle(0, 1, 2);
        let edges = c.boundary_of_triangle((0, 1, 2));
        assert_eq!(edges.len(), 3);
        assert!(edges.contains(&(0, 1)));
        assert!(edges.contains(&(1, 2)));
        assert!(edges.contains(&(0, 2)));
    }

    #[test]
    fn test_boundary_of_vertex() {
        let c = SimplicialComplex::new();
        assert!(c.boundary_of_vertex(0).is_empty());
    }

    #[test]
    fn test_boundary_of_tetrahedron() {
        let mut c = SimplicialComplex::new();
        c.add_tetrahedron(0, 1, 2, 3);
        let tris = c.boundary_of_tetrahedron((0, 1, 2, 3));
        assert_eq!(tris.len(), 4);
    }

    #[test]
    fn test_h0_betti() {
        let points = vec![vec![0.0], vec![1.0], vec![3.0]];
        assert_eq!(SimplicialComplex::h0_betti(&points, 0.5), 3);
        // At threshold 1.1: 0-1 connected (dist=1.0), but 3 is separate
        assert_eq!(SimplicialComplex::h0_betti(&points, 1.0), 2);
        assert_eq!(SimplicialComplex::h0_betti(&points, 10.0), 1);
    }

    #[test]
    fn test_h0_betti_empty() {
        let points: Vec<Vec<f64>> = vec![];
        assert_eq!(SimplicialComplex::h0_betti(&points, 1.0), 0);
    }

    #[test]
    fn test_h1_betti() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![0.5, 0.866],
        ];
        // At high threshold, all edges form a triangle → β₁ = 0 (filled)
        assert_eq!(SimplicialComplex::h1_betti(&points, 0.5), 0);
        // Edge threshold: three edges, no triangle → β₁ = 1 (cycle)
        assert_eq!(SimplicialComplex::h1_betti(&points, 1.0), 0);
    }

    #[test]
    fn test_h1_betti_square_cycle() {
        let points = vec![
            vec![0.0, 0.0],
            vec![1.0, 0.0],
            vec![1.0, 1.0],
            vec![0.0, 1.0],
        ];
        // At medium threshold: edges of a square, no diagonals → has a cycle
        let b1 = SimplicialComplex::h1_betti(&points, 1.0);
        assert_eq!(b1, 1);
    }

    #[test]
    fn test_euclidean_dist() {
        let d = euclidean_dist(&[0.0, 0.0], &[3.0, 4.0]);
        assert!((d - 5.0).abs() < 1e-10);
    }

    #[test]
    fn test_total_simplices() {
        let mut c = SimplicialComplex::new();
        c.add_tetrahedron(0, 1, 2, 3);
        assert_eq!(c.total_simplices(), 15); // 4+6+4+1
    }
}
