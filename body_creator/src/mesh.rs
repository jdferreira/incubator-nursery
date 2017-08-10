use std::collections::{HashSet,HashMap};
use std::collections::hash_map::Entry;
use vector::Vector3D;

#[inline(always)]
fn new_face(a: usize, b: usize, c:usize) -> (usize, usize, usize) {
    if      a < b && a < c    { (a, b, c) }
    else if b < a && b < c    { (b, c, a) }
    else /* c is the least */ { (c, a, b) }
}

/// A Mesh contains vertices, edges and faces.
pub struct Mesh {
    /// The positions of each vertex. Vertices are identified by their index in this vector.
    pub positions: Vec<Vector3D>,
    
    /// The faces in the mesh
    pub faces: HashSet<(usize, usize, usize)>,
    
    /// Each element contains the neighbours of the vertex. The key in the map is the neighbour
    /// (allows for quick checking of neighbours) and the value is a tuple, where the first value
    /// is the left vertex of this edge and the second value is the right vertex
    pub neighbours: Vec<HashMap<usize, (usize, usize)>>
}


impl Mesh {
    
    /// Creates a new mesh with the given vertices and faces. The method assumes correctness of the
    /// arguments. In particular, all indices in `faces` must be between 0 and
    /// `positions.len() - 1`. Otherwise, a panic occurs.
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
        
        Mesh {
            positions: positions,
            faces: faces,
            neighbours: neighbours,
        }
    }
    
    // Function to make a face
    #[inline(always)]
    fn make_face(&mut self, a: usize, b: usize, c: usize) {
        self.faces.insert(new_face(a, b, c));
    }
    
    // Function to make an edge
    #[inline(always)]
    fn make_edge(&mut self, start: usize, end: usize, left: usize, right: usize) {
        self.neighbours[start].insert(end, (left, right));
        self.neighbours[end].insert(start, (right, left));
    }
    
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
            self.faces.remove(&new_face(curr, next, up));
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


// const MAX_IDX: usize = ::std::usize::MAX;

// /// A SymmetricMesh is a mesh that is symmetric around the plane `x = 0`. The main functional
// /// difference between a Mesh and a SymmetricMesh is in the `extrude` method, which must take care
// /// of extrusions around the symmetry plane in the SymmetricMesh
// pub struct SymmetricMesh {
//     /// The positions of each vertex. Vertices are identified by their index in this vector.
//     pub positions: Vec<Vector3D>,
    
//     /// The set of vertices that are locked in the symmetry plane
//     pub locked: HashSet<usize>,
    
//     /// The faces in the mesh
//     pub faces: HashSet<(usize, usize, usize)>,
    
//     /// Each element contains the neighbours of the vertex. The key in the map is the neighbour
//     /// (allows for quick checking of neighbours) and the value is a tuple, where the first value
//     /// is the left vertex of this edge and the second value is the right vertex. The edges that
//     /// are between locked vertices do not have either a left or a right adjacent vertex, which is
//     /// marked with the `std::usize::MAX` value
//     pub neighbours: Vec<HashMap<usize, (usize, usize)>>
// }


// impl SymmetricMesh {
    
//     // Function to make a face
//     #[inline(always)]
//     fn make_face(&mut self, a: usize, b: usize, c: usize) {
//         self.faces.insert(new_face(a, b, c));
//     }
    
//     // Function to make an edge
//     #[inline(always)]
//     fn make_edge(&mut self, start: usize, end: usize, left: usize, right: usize) {
//         self.neighbours[start].insert(end, (left, right));
//         self.neighbours[end].insert(start, (right, left));
//     }
    
//     /// There can be two types of rings: one where the ring crosses the symmetry plane and one
//     /// that does not. A ring that crosses the symmetry plane must have a key that is
//     /// `std::usize::MAX`, which points to one of the vertices that is locked to the symmetry
//     /// plane; the last vertex on the ring points again to a `std::usize::MAX` value, signalling
//     /// the end of the ring. This method does not check that the ring is valid and assumes it is.
//     /// The ring must be given such that the up neighbour (FIXME I don't remember which one) is
//     /// inside the region defined by the ring.
//     pub fn extrude(&mut self, ring: &HashMap<usize, usize>, region: &[usize], dir: Vector3D) {
//         // For each vertex in the ring, make a new one just beside it (shifted by `dir`). Also, for
//         // each of these vertices, make a new one between it and its next one on the ring, also
//         // shifted with `dir`
//         let mut new_ring: HashMap<usize, usize> = HashMap::new();
//         let mut wcn_indices: HashMap<usize, (usize, usize)> = HashMap::new();
        
//         for (&curr, &next) in ring.iter() {
//             if curr == MAX_IDX || next == MAX_IDX { continue; }
            
//             let new_pos = self.positions[curr] + dir;
//             self.positions.push(new_pos);
//             let new_idx = self.positions.len() - 1;
//             new_ring.insert(curr, new_idx);
//             if self.locked.contains(&curr) {
//                 self.locked.insert(new_idx);
//             }
            
//             // Make the vertex between curr and next, but only if the curr or next are not locked
//             if !self.locked.contains(&curr) && !self.locked.contains(&next) {
//                 let wcn_pos = (self.positions[curr] + self.positions[next] + dir) / 2.0;
//                 self.positions.push(wcn_pos);
//                 let wcn = self.positions.len() - 1;
                
//                 wcn_indices.entry(curr).or_insert((0, wcn)).1 = wcn;
//                 wcn_indices.entry(next).or_insert((wcn, 0)).0 = wcn;
//             }
//         }
        
//         // Also move the rest of the region in that direction
//         for &index in region.iter() {
//             if ring.contains_key(&index) { continue; }
//             self.positions[index] = self.positions[index] + dir;
//         }
        
//         // Go around the ring and cut and glue the vertices with new edges
//         for (&curr, &next) in ring.iter() {
//             let up = self.neighbours[curr][&next].0;
//             let down1 = {
//                 let d = self.neighbours[next][&up].1;
//                 *new_ring.get(&d).unwrap_or(&d)
//             };
//             let down2 = {
//                 let d = self.neighbours[up][&curr].1;
//                 *new_ring.get(&d).unwrap_or(&d)
//             };
            
//             let wcn = wcn_indices[&curr].1;
//             let w_c = wcn_indices[&curr].0;
//             let wn_ = wcn_indices[&next].1;
            
//             // Remove the face (curr, next, up) and the edges (next, up) and (up, curr)
//             self.faces.remove(&new_face(curr, next, up));
//             self.neighbours[next].remove(&up);
//             self.neighbours[up].remove(&next);
//             self.neighbours[up].remove(&curr);
//             self.neighbours[curr].remove(&up);
            
//             // Make face (curr', next', up) (where A' means new_ring[A])
//             self.make_face(new_ring[&curr], new_ring[&next], up);
//             self.make_edge(new_ring[&next], up, new_ring[&curr], down1);
//             self.make_edge(up, new_ring[&curr], new_ring[&next], down2);
//             self.make_edge(new_ring[&curr], new_ring[&next], up, wcn);
            
//             // Make the new four faces ...
//             self.make_face(curr, next, wcn);
//             self.make_face(new_ring[&next], new_ring[&curr], wcn);
//             self.make_face(next, wn_, wcn);
//             self.make_face(wcn, wn_, new_ring[&next]);
            
//             // ... and the new five edges
//             self.make_edge(next, wcn, curr, wn_);
//             self.make_edge(wcn, curr, next, w_c);
//             self.make_edge(wcn, new_ring[&next], new_ring[&curr], wn_);
//             self.make_edge(new_ring[&curr], wcn, new_ring[&next], w_c);
//             self.make_edge(wcn, wn_, new_ring[&next], next);
            
//             // Update the adjacent vertices on the (curr, next) edge
//             self.neighbours[curr].get_mut(&next).unwrap().0 = wcn;
//             self.neighbours[next].get_mut(&curr).unwrap().1 = wcn;
//         }
//     }
// }
