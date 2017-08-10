use super::*;

#[test]
fn test_get_inside_vertices() {
    let positions = vec![Vector3D(0.0, 0.0, 0.0); 11];
    let faces = vec![(0, 1, 3), (1, 2, 4), (1, 4, 3), (3, 7, 6), (3, 4, 7), (4, 8, 7),
                     (4, 5, 8), (6, 7, 9), (9, 7, 10), (7, 8, 10)];
    let mesh = Mesh::new(positions, faces);
    let ring = Ring::new(&[0, 1, 2, 4, 5, 8, 10, 9, 6, 3]);
    let mut result = mesh.get_inside_vertices(&ring);
    result.sort();
    
    assert_eq!(result, vec![7]);
    
    let positions = vec![Vector3D(0.0, 0.0, 0.0); 16];
    let faces = vec![(0, 1, 4),
                     (1, 5, 4),
                     (1, 2, 5),
                     (2, 6, 5),
                     (2, 3, 7),
                     (2, 7, 6),
                     (4, 5, 9),
                     (4, 9, 8),
                     (5, 6, 10),
                     (5, 10, 9),
                     (6, 7, 10),
                     (7, 11, 10),
                     (8, 9, 13),
                     (8, 13, 12),
                     (9, 10, 13),
                     (10, 14, 13),
                     (10, 11, 14),
                     (11, 15, 14)];
    let mesh = Mesh::new(positions, faces);
    let ring = Ring::new(&[0, 1, 2, 3, 7, 11, 15, 14, 13, 12, 8, 4]);
    let mut result = mesh.get_inside_vertices(&ring);
    result.sort();
    
    assert_eq!(result, vec![5, 6, 9, 10]);
    
    let positions = vec![Vector3D(0.0, 0.0, 0.0); 17];
    let faces = vec![(0, 1, 4),
                     (1, 5, 4),
                     (1, 2, 5),
                     (2, 6, 5),
                     (2, 3, 7),
                     (2, 7, 6),
                     (4, 5, 9),
                     (4, 9, 8),
                     (5, 6, 16),
                     (6, 10, 16),
                     (10, 9, 16),
                     (9, 5, 16),
                     (6, 7, 10),
                     (7, 11, 10),
                     (8, 9, 13),
                     (8, 13, 12),
                     (9, 10, 13),
                     (10, 14, 13),
                     (10, 11, 14),
                     (11, 15, 14)];
    let mesh = Mesh::new(positions, faces);
    let ring = Ring::new(&[0, 1, 2, 3, 7, 11, 15, 14, 13, 12, 8, 4]);
    let mut result = mesh.get_inside_vertices(&ring);
    result.sort();
    
    assert_eq!(result, vec![5, 6, 9, 10, 16]);
}

#[test]
fn test_extrude() {
    let positions = vec![Vector3D(0.0, 0.0, 0.0),
                         Vector3D(1.0, 0.0, 0.0),
                         Vector3D(3.0, 0.0, 0.0),
                         Vector3D(4.0, 0.0, 0.0),
                         Vector3D(0.0, 1.0, 0.0),
                         Vector3D(3.0, 1.0, 0.0),
                         Vector3D(2.0, 2.0, 0.0),
                         Vector3D(0.0, 3.0, 0.0),
                         Vector3D(1.0, 3.0, 0.0),
                         Vector3D(3.0, 3.0, 0.0),
                         Vector3D(4.0, 3.0, 0.0),
                         Vector3D(0.0, 0.0, -9.0)];
    let faces = vec![(0, 1, 4),
                     (1, 2, 5),
                     (2, 3, 5),
                     (1, 6, 4),
                     (1, 5, 6),
                     (5, 3, 10),
                     (6, 5, 10),
                     (4, 6, 7),
                     (7, 6, 8),
                     (8, 6, 9),
                     (9, 6, 10),
                     (0, 11, 1),
                     (1, 11, 2),
                     (2, 11, 3),
                     (3, 11, 10),
                     (10, 11, 9),
                     (9, 11, 8),
                     (8, 11, 7),
                     (7, 11, 4),
                     (4, 11, 0)];
    let mut mesh = Mesh::new(positions, faces);
    let ring = Ring::new(&[0, 1, 2, 3, 10, 9, 8, 7, 4]);
    
    mesh.extrude(&ring, Vector3D(0.0, 0.0, 1.0));
}
