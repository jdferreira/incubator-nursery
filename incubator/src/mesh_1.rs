#![allow(dead_code)]

use std::collections::HashMap;
use vector::Vector3D;

/// Internal (private) index for vertices in the mesh data structure.
/// We use this so that we do not accidentally pass an edge index or face index where a vertex
/// index is needed.
struct VertexIdx(usize);

/// We allow direct conversion from usize to VertexIdx
impl From<usize> for VertexIdx {
    fn from(idx: usize) -> VertexIdx {
        VertexIdx(idx)
    }
}

/// Internal (private) index for edges in the mesh data structure.
/// We use this so that we do not accidentally pass a vertex index or face index where an edge
/// index is needed.
struct EdgeIdx(usize);

/// Internal (private) index for faces in the mesh data structure.
/// We use this so that we do not accidentally pass a vertex index or an edge index where a face
/// index is needed.
struct FaceIdx(usize);


/// Represents a vertex in the mesh data structure.
/// The vertex contains a 3D vector defining its geometric position, as well as a list of edges
/// that start or end in this vertex.
struct Vertex {
    
    /// The spatial position of the vertex
    coordinates: Vector3D,
    
    /// The edges that start or end in this vertex.
    /// The edges are stored in counter-clockwise order and any mutation to this field must
    /// keep that invariant
    edges: Vec<EdgeIdx>,
    
}

/// Represents an edge in the mesh data structure.
/// We represent the mesh using the winged-edge data structure, which means that the edge contains
/// pointers to:
///   - the start and end vertices of the edge;
///   - the left and right faces connected to this edge (left and right are related to the direction
///     defined by the start and end vertices);
///   - the four edges that make up the "butterfly wings" of the edge: the two left ones and the two
///     right ones. Given that faces are always triangles, the two left (resp. right) faces share
///     a vertex.
struct Edge {
    
    /// The vertex where this edge starts. The direction from start to end defines the orientation
    /// of the edge and, therefore, the meaning of the rest of the fields in this struct.
    start: VertexIdx,
    
    /// The vertex where this edge ends. The direction from start to end defines the orientation
    /// of the edge and, therefore, the meaning of the rest of the fields in this struct.
    end: VertexIdx,
    
    /// The face to the left of the the edge.
    left_face: FaceIdx,
    
    /// The face to the right of the the edge.
    right_face: FaceIdx,
    
    /// The edge that:
    ///   - starts in the `start` vertex,
    ///   - goes to the left of the edge (i.e. it is adjacent to the left face), and
    ///   - is closest (in counter-clockwise rotation) to this edge.
    start_left_vertex: EdgeIdx,
    
    /// The edge that:
    ///   - starts in the `start` vertex,
    ///   - goes to the right of the edge (i.e. it is adjacent to the right face), and
    ///   - is closest (in clockwise rotation) to this edge.
    start_right_vertex: EdgeIdx,
    
    /// The edge that:
    ///   - starts in the `end` vertex,
    ///   - goes to the left of the edge (i.e. it is adjacent to the left face), and
    ///   - is closest (in clockwise rotation) to this edge.
    end_left_vertex: EdgeIdx,
    
    /// The edge that:
    ///   - starts in the `end` vertex,
    ///   - goes to the right of the edge (i.e. it is adjacent to the right face), and
    ///   - is closest (in counter-clockwise rotation) to this edge.
    end_right_vertex: EdgeIdx,
    
}

/// Represents a face in the mesh data structure
/// Each face is represented by the three vertices that compose it. The vertices are given in
/// counter-clockwise order, with any of the three the first one. This means that two `Face`s
/// represent the same face if the three vertices are equal up to rotation:
///
/// ```
/// assert!(Face(1, 2, 3) == Face(2, 3, 1));
/// assert!(Face(1, 2, 3) != Face(2, 1, 3));
/// ```
struct Face(VertexIdx, VertexIdx, VertexIdx);

/// Two `Face`s represent the same face if the three vertices are equal up to rotation:
///
/// ```
/// assert!(Face(1, 2, 3) == Face(2, 3, 1));
/// assert!(Face(1, 2, 3) != Face(2, 1, 3));
/// ```
impl PartialEq for Face {
    fn eq(&self, other: &Face) -> bool {
        let &Face(VertexIdx(self_a ), VertexIdx(self_b ), VertexIdx(self_c )) = self ;
        let &Face(VertexIdx(other_a), VertexIdx(other_b), VertexIdx(other_c)) = other;
        
        if self_a == other_a {
            self_b == other_b && self_c == other_c
        }
        else if self_a == other_b {
            self_b == other_c && self_c == other_a
        }
        else if self_a == other_c {
            self_b == other_a && self_c == other_b
        }
        else {
            false
        }
    }
}

/// Two `Face`s represent the same face if the three vertices are equal up to rotation:
///
/// ```
/// assert!(Face(1, 2, 3) == Face(2, 3, 1));
/// assert!(Face(1, 2, 3) != Face(2, 1, 3));
/// ```
impl Eq for Face {}

/// Represents the entire mesh.
struct Mesh {
    
    /// The vertices of the mesh. Vertices are stored with an index, which is used throughout the
    /// vertices, edges and faces to identify the vertex. This are publicly accessed through a
    /// `VertexIdx` struct.
    vertices: Vec<Vertex>,
    
    /// The edges of the mesh. Edges are stored with an index, which is used throughout the
    /// vertices, edges and faces to identify the edge. This are publicly accessed through a
    /// `EdgeIdx` struct.
    edges: HashMap<(VertexIdx, VertexIdx), Edge>,
    
    /// The faces of the mesh. Faces are stored with an index, which is used throughout the
    /// vertices, edges and faces to identify the face. This are publicly accessed through a
    /// `FaceIdx` struct.
    faces: Vec<Face>,
    
}

/// Initializes a list of vertices based on their position. The vertices will be empty, in the sense
/// that no edge information will be included in the returned vector.
fn init_vertices(positions: &[Vector3D]) -> Vec<Vertex> {
    positions.iter().map(|&v| Vertex { coordinates: v, edges: vec![] }).collect()
}

/// Initializes a list of edges based on the faces. To avoid duplicates, we will create only the
/// edges where the start vertex is smaller than the end vertex. The function requires a second
/// argument that specifies the number of edges to create. This is used to directly allocate the
/// space needed for the vector of edges and must be calculated upstream.
fn init_edges(faces: &[(usize, usize, usize)], n_edges: usize) -> Vec<Edge> {
    let mut edges = Vec::with_capacity(n_edges);
    
    edges
}

/// Initializes a list of faces based on the index of its vertices.
fn init_faces(faces: &[(usize, usize, usize)]) -> Vec<Face> {
    faces.iter().map(|&(a, b, c)| Face(a.into(), b.into(), c.into())).collect()
}

impl Mesh {
    
    /// Creates a new mesh with the given vertices and faces. The method assumes correctness of the
    /// arguments. In particular, all indices in `faces` must be between 0 and
    /// `positions.len() - 1`, and they must be given in counter-clockwise order.
    pub fn new(positions: Vec<Vector3D>, faces: Vec<(usize, usize, usize)>) -> Mesh {
        
        let n_edges = faces.len() + positions.len() - 2;
        
        let vertices = init_vertices(&positions);
        let edges = init_edges(&faces, n_edges);
        let faces = init_faces(&faces);
        
        
        
        return Mesh {
            edges: edges,
            faces: faces,
            vertices: vertices,
        };
        
        // let mut neighbours = vec![HashMap::new(); positions.len()];
        
        // for &(i1, i2, i3) in faces.iter() {
        //     for &(a, b, c) in [(i1, i2, i3), (i2, i3, i1), (i3, i1, i2)].iter() {
        //         let mut map = neighbours.get_mut(a).unwrap();
        //         match map.entry(b) {
        //             Entry::Vacant(entry)   => { entry.insert((c, 0)); },
        //             Entry::Occupied(entry) => { entry.into_mut().0 = c; },
        //         }
        //         match map.entry(c) {
        //             Entry::Vacant(entry)   => { entry.insert((0, b)); },
        //             Entry::Occupied(entry) => { entry.into_mut().1 = b; },
        //         }
        //     }
        // }
        
        // let mut internal_faces = HashMap::new();
        // for &(i1, i2, i3) in faces.iter() {
        //     // Sort the indices in ascending order
        //     let vertices = sort_vertices(i1, i2, i3);
        //     let face = Face::new(&positions, vertices);
        //     internal_faces.insert(vertices, face);
        // }
        
        // Mesh {
        //     positions: positions,
        //     faces: internal_faces,
        //     neighbours: neighbours,
        // }
    }
    
}
