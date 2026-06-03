//! Integration tests for persistent-sheaf

use persistent_sheaf::*;

// ---- Simplicial Complex ----

#[test]
fn test_simplicial_complex_empty() {
    let sc = SimplicialComplex::new();
    assert_eq!(sc.num_simplices(0), 0);
    assert_eq!(sc.euler_characteristic(), 0);
}

#[test]
fn test_simplicial_complex_triangle() {
    let mut sc = SimplicialComplex::new();
    sc.add_triangle(0, 1, 2);
    assert_eq!(sc.num_simplices(0), 3);
    assert_eq!(sc.num_simplices(1), 3);
    assert_eq!(sc.num_simplices(2), 1);
    assert_eq!(sc.euler_characteristic(), 1);
}

#[test]
fn test_simplicial_complex_two_triangles() {
    let mut sc = SimplicialComplex::new();
    sc.add_triangle(0, 1, 2);
    sc.add_triangle(1, 2, 3);
    assert_eq!(sc.num_simplices(0), 4);
    assert_eq!(sc.num_simplices(1), 5);
    assert_eq!(sc.num_simplices(2), 2);
    assert_eq!(sc.euler_characteristic(), 1);
}

#[test]
fn test_simplicial_complex_dedup() {
    let mut sc = SimplicialComplex::new();
    sc.add_vertex(0);
    sc.add_vertex(0);
    assert_eq!(sc.num_simplices(0), 1);
    sc.add_edge(0, 1);
    sc.add_edge(1, 0);
    assert_eq!(sc.num_simplices(1), 1);
}

#[test]
fn test_tetrahedron() {
    let mut sc = SimplicialComplex::new();
    sc.add_tetrahedron(0, 1, 2, 3);
    assert_eq!(sc.num_simplices(0), 4);
    assert_eq!(sc.num_simplices(1), 6);
    assert_eq!(sc.num_simplices(2), 4);
    assert_eq!(sc.num_simplices(3), 1);
}

#[test]
fn test_h0_betti_integration() {
    let points = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![3.0, 0.0],
    ];
    assert_eq!(SimplicialComplex::h0_betti(&points, 0.5), 3);
    assert_eq!(SimplicialComplex::h0_betti(&points, 1.1), 2);
    assert_eq!(SimplicialComplex::h0_betti(&points, 2.1), 1);
}

#[test]
fn test_h1_betti_integration() {
    let points = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![1.0, 1.0],
        vec![0.0, 1.0],
    ];
    // Near threshold: square with no diagonals = cycle
    assert_eq!(SimplicialComplex::h1_betti(&points, 1.0), 1);
}

// ---- Persistence Diagram ----

#[test]
fn test_persistence_diagram_basic() {
    let mut pd = PersistenceDiagram::new();
    assert!(pd.is_empty());
    pd.add(0.0, 1.0, 0);
    pd.add(0.5, 2.0, 1);
    assert_eq!(pd.len(), 2);

    let dim0 = pd.filter_dimension(0);
    assert_eq!(dim0.len(), 1);
    assert!((dim0[0].persistence() - 1.0).abs() < 1e-10);
}

#[test]
fn test_persistence_pair_essential() {
    let mut pd = PersistenceDiagram::new();
    pd.add(0.0, f64::INFINITY, 0);
    assert!(pd.pairs[0].is_essential());
    assert!(pd.pairs[0].persistence().is_infinite());
}

#[test]
fn test_persistence_pair_midpoint() {
    let mut pd = PersistenceDiagram::new();
    pd.add(1.0, 5.0, 0);
    assert!((pd.pairs[0].midpoint() - 3.0).abs() < 1e-10);
}

#[test]
fn test_total_persistence() {
    let mut pd = PersistenceDiagram::new();
    pd.add(0.0, 2.0, 0);
    pd.add(1.0, 4.0, 1);
    let tp = pd.total_persistence(1.0);
    assert!((tp - 5.0).abs() < 1e-10);
}

#[test]
fn test_most_persistent() {
    let mut pd = PersistenceDiagram::new();
    pd.add(0.0, 1.0, 0);
    pd.add(0.0, 5.0, 1);
    let mp = pd.most_persistent().unwrap();
    assert!((mp.persistence() - 5.0).abs() < 1e-10);
}

#[test]
fn test_bottleneck_distance_integration() {
    let mut d1 = PersistenceDiagram::new();
    d1.add(0.0, 1.0, 0);
    let mut d2 = PersistenceDiagram::new();
    d2.add(0.1, 1.1, 0);
    let dist = d1.bottleneck_distance(&d2);
    assert!(dist < 0.2);
}

// ---- Sheaf ----

#[test]
fn test_constant_sheaf() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    sc.add_edge(1, 2);
    let sheaf = CellularSheaf::constant(sc, 2);
    assert_eq!(sheaf.stalk_dimension, 2);
    assert_eq!(sheaf.restriction_maps.len(), 2);
}

#[test]
fn test_sheaf_from_weights() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    let sheaf = CellularSheaf::from_weights(sc, &[0.5]);
    assert_eq!(sheaf.stalk_dimension, 1);
}

#[test]
fn test_global_sections_integration() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    let sheaf = CellularSheaf::constant(sc, 3);
    assert_eq!(sheaf.global_section_dimension(), 3);
}

// ---- Sheaf Laplacian ----

#[test]
fn test_sheaf_laplacian_integration() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    let sheaf = CellularSheaf::constant(sc, 2);
    let lap = SheafLaplacian::from_sheaf(&sheaf);
    assert_eq!(lap.dimension, 4);
}

#[test]
fn test_edge_laplacian_integration() {
    let mut sc = SimplicialComplex::new();
    sc.add_edge(0, 1);
    let sheaf = CellularSheaf::constant(sc, 1);
    let el = SheafLaplacian::edge_laplacian(&sheaf, 0);
    assert_eq!(el.dimension, 2);
}

#[test]
fn test_triangle_laplacian_integration() {
    let mut sc = SimplicialComplex::new();
    sc.add_triangle(0, 1, 2);
    let sheaf = CellularSheaf::constant(sc, 1);
    let tl = SheafLaplacian::triangle_laplacian(&sheaf, 0);
    assert_eq!(tl.dimension, 3);
}

// ---- Filtration ----

#[test]
fn test_filtration_from_point_cloud() {
    let points = vec![
        vec![0.0, 0.0],
        vec![1.0, 0.0],
        vec![2.0, 0.0],
    ];
    let f = Filtration::from_point_cloud(&points, 5);
    assert_eq!(f.len(), 6);
}

#[test]
fn test_filtration_thresholds() {
    let points = vec![vec![0.0], vec![1.0]];
    let f = Filtration::from_point_cloud(&points, 10);
    let thresh = f.thresholds();
    assert_eq!(thresh.len(), 11);
    assert!((thresh[0] - 0.0).abs() < 1e-10);
    assert!((thresh[10] - 1.0).abs() < 1e-10);
}
