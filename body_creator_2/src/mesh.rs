use std::collections::{HashSet,HashMap};
use std::collections::hash_map::Entry;

use vector::Vector3D;

#[inline(always)]
fn sort_vertices(a: usize, b: usize, c:usize) -> (usize, usize, usize) {
    if      a < b && a < c    { (a, b, c) }
    else if b < a && b < c    { (b, c, a) }
    else /* c < a && c < b */ { (c, a, b) }
}

pub struct Face {
    vertices: (usize, usize, usize),
    min_w: f64,
    min_h: f64,
    history: Vec<Option<FaceHistoryStep>>,
}

impl Face {
    fn new(positions: &[Vector3D], vertices: (usize, usize, usize)) -> Face {
        let step = FaceHistoryStep::new(positions, vertices);
        Face {
            vertices: vertices,
            min_w: step.width,
            min_h: step.height,
            history: vec![Some(step)]
        }
    }
    
    fn update_present(&mut self, positions: &[Vector3D]) {
        let step = FaceHistoryStep::new(positions, self.vertices);
        self.min_w = self.min_w.max(step.width);
        self.min_h = self.min_h.max(step.height);
        self.history.push(Some(step));
    }
    
    fn update_absent(&mut self) {
        self.history.push(None);
    }
}

struct FaceHistoryStep {
    /// The index of the first vertex in the mesh
    i1: usize,
    /// The index of the second vertex in the mesh
    i2: usize,
    /// The index of the third vertex in the mesh
    i3: usize,
    /// The x-coordinate of the second vertex of the triangle in UV space (width of the triangle)
    width: f64,
    /// The y-coordinate of the third vertex of the triangle in UV space (height of the triangle)
    height: f64,
    /// The x-coordinate of the third vertex of the triangle in UV space
    mid: f64,
}

impl FaceHistoryStep {
    fn new(positions: &[Vector3D], (i1, i2, i3): (usize, usize, usize)) -> FaceHistoryStep {
        // Get the position of the three vertices
        let v1 = positions[i1];
        let v2 = positions[i2];
        let v3 = positions[i3];
        
        // Get the sides of the trinagle
        let a = (v1 - v2).norm();
        let b = (v2 - v3).norm();
        let c = (v3 - v1).norm();
        
        // Calculate the area of the triangle
        let p = (a + b + c ) / 2.0;
        let area = (p * (p - a) * (p - b) * (p - c)).sqrt();
        
        // Depeding on which is the largest side, return something a bit different
        if a > b && a > c {
            let h = 2.0 * area / a;
            let m = (c * c - h * h).sqrt();
            FaceHistoryStep { i1: i1, i2: i2, i3: i3, width: a, height: h, mid: m }
        }
        else if b > a && b > c {
            let h = 2.0 * area / b;
            let m = (a * a - h * h).sqrt();
            FaceHistoryStep { i1: i2, i2: i3, i3: i1, width: b, height: h, mid: m }
        }
        else {
            let h = 2.0 * area / c;
            let m = (b * b - h * h).sqrt();
            FaceHistoryStep { i1: i3, i2: i1, i3: i2, width: c, height: h, mid: m }
        }
    }
}

/// A Mesh contains vertices, edges and faces.
pub struct Mesh {
    /// The positions of each vertex. Vertices are identified by their index in this vector.
    pub positions: Vec<Vector3D>,
    
    /// The faces in the mesh
    pub faces: HashMap<(usize, usize, usize), Face>,
    
    /// Each element contains the neighbours of the vertex. The key in the map is the neighbour
    /// (allows for quick checking of neighbours) and the value is a tuple, where the first value
    /// is the left vertex of this edge and the second value is the right vertex
    pub neighbours: Vec<HashMap<usize, (usize, usize)>>,
}


impl Mesh {
    /// Creates a new mesh with the given vertices and faces. The method assumes correctness of the
    /// arguments. In particular, all indices in `faces` must be between 0 and
    /// `positions.len() - 1`.
    pub fn new(positions: Vec<Vector3D>, faces: HashSet<(usize, usize, usize)>) -> Mesh {
        let mut neighbours = vec![HashMap::new(); positions.len()];
        
        for &(i1, i2, i3) in faces.iter() {
            for &(a, b, c) in [(i1, i2, i3), (i2, i3, i1), (i3, i1, i2)].iter() {
                let mut map = neighbours.get_mut(a).unwrap();
                match map.entry(b) {
                    Entry::Vacant(entry)   => { entry.insert((c, 0)); },
                    Entry::Occupied(entry) => { entry.into_mut().0 = c; },
                }
                match map.entry(c) {
                    Entry::Vacant(entry)   => { entry.insert((0, b)); },
                    Entry::Occupied(entry) => { entry.into_mut().1 = b; },
                }
            }
        }
        
        let mut internal_faces = HashMap::new();
        for &(i1, i2, i3) in faces.iter() {
            // Sort the indices in ascending order
            let vertices = sort_vertices(i1, i2, i3);
            let face = Face::new(&positions, vertices);
            internal_faces.insert(vertices, face);
        }
        
        Mesh {
            positions: positions,
            faces: internal_faces,
            neighbours: neighbours,
        }
    }
    
    // // Function to make a face
    // #[inline(always)]
    // fn make_face(&mut self, a: usize, b: usize, c: usize) {
    //     let (a, b, c) = sort_vertices(a, b, c);
    //     let uv_info = UVInfo::new(self.positions[a], self.positions[b], self.positions[c]);
    //     self.faces.insert(Face::new((a, b, c), uv_info));
    // }
    
    // // Function to make an edge
    // #[inline(always)]
    // fn make_edge(&mut self, start: usize, end: usize, left: usize, right: usize) {
    //     self.neighbours[start].insert(end, (left, right));
    //     self.neighbours[end].insert(start, (right, left));
    // }
    
    /// The ring must be given such that the up neighbour (FIXME I don't remember which one) is
    /// inside the region defined by the ring.
    #[allow(dead_code)]
    pub fn extrude(&mut self, ring: &HashMap<usize, usize>, region: &[usize], dir: Vector3D) {
        // For each vertex in the ring, make a new one just beside it (shifted by `dir`). Also, for
        // each of these vertices, make a new one between it and its next one on the ring, also
        // shifted with `dir`
        let mut new_ring: HashMap<usize, usize> = HashMap::new();
        let mut wcn_indices: HashMap<usize, (usize, usize)> = HashMap::new();
        
        for (&curr, &next) in ring.iter() {
            let new_pos = self.positions[curr] + dir;
            self.positions.push(new_pos);
            new_ring.insert(curr, self.positions.len() - 1);
            
            // Make the vertex between curr and next
            let wcn_pos = (self.positions[curr] + self.positions[next] + dir) / 2.0;
            self.positions.push(wcn_pos);
            let wcn = self.positions.len() - 1;
            
            wcn_indices.entry(curr).or_insert((0, wcn)).1 = wcn;
            wcn_indices.entry(next).or_insert((wcn, 0)).0 = wcn;
        }
        
        // Also move the rest of the region in that direction
        for &index in region.iter() {
            if ring.contains_key(&index) { continue; }
            self.positions[index] = self.positions[index] + dir;
        }
        
        // Go around the ring and cut and glue the vertices with new edges
        for (&curr, &next) in ring.iter() {
            let up = self.neighbours[curr][&next].0;
            let down1 = {
                let d = self.neighbours[next][&up].1;
                *new_ring.get(&d).unwrap_or(&d)
            };
            let down2 = {
                let d = self.neighbours[up][&curr].1;
                *new_ring.get(&d).unwrap_or(&d)
            };
            
            let wcn = wcn_indices[&curr].1;
            let w_c = wcn_indices[&curr].0;
            let wn_ = wcn_indices[&next].1;
            
            // Remove the face (curr, next, up) and the edges (next, up) and (up, curr)
            self.faces.remove(&sort_vertices(curr, next, up));
            self.neighbours[next].remove(&up);
            self.neighbours[up].remove(&next);
            self.neighbours[up].remove(&curr);
            self.neighbours[curr].remove(&up);
            
            // Make face (curr', next', up) (where A' means new_ring[A])
            self.make_face(new_ring[&curr], new_ring[&next], up);
            self.make_edge(new_ring[&next], up, new_ring[&curr], down1);
            self.make_edge(up, new_ring[&curr], new_ring[&next], down2);
            self.make_edge(new_ring[&curr], new_ring[&next], up, wcn);
            
            // Make the new four faces ...
            self.make_face(curr, next, wcn);
            self.make_face(new_ring[&next], new_ring[&curr], wcn);
            self.make_face(next, wn_, wcn);
            self.make_face(wcn, wn_, new_ring[&next]);
            
            // ... and the new five edges
            self.make_edge(next, wcn, curr, wn_);
            self.make_edge(wcn, curr, next, w_c);
            self.make_edge(wcn, new_ring[&next], new_ring[&curr], wn_);
            self.make_edge(new_ring[&curr], wcn, new_ring[&next], w_c);
            self.make_edge(wcn, wn_, new_ring[&next], next);
            
            // Update the adjacent vertices on the (curr, next) edge
            self.neighbours[curr].get_mut(&next).unwrap().0 = wcn;
            self.neighbours[next].get_mut(&curr).unwrap().1 = wcn;
        }
    }
}

