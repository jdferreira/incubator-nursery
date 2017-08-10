#![allow(dead_code)]

mod ring;
mod extrusion;

#[cfg(test)] mod tests;

use std::collections::HashMap;
use std::collections::hash_map::Entry;

use vector::Vector3D;
use self::ring::Ring;

/// Represents the entire mesh.
pub struct Mesh {
    /// The 3-dimensional position of each vertex.
    /// Vertices are identified by their index in this vector.
    pub positions: Vec<Vector3D>,
    
    /// The faces in the mesh
    pub faces: Vec<(usize, usize, usize)>,
    
    /// Each element contains the neighbours of the vertex. The key in the map is the neighbour
    /// (allows for quick checking of neighbours) and the value is a tuple, where the first value
    /// is the left vertex of this edge and the second value is the right vertex
    pub neighbours: Vec<HashMap<usize, EdgeData>>,
}

/// Each edge has a left vertex and a right vertex, which correspond to the third vertex in the face
/// to the left of the edge and to the third vertex in the face to the right of the edge
#[derive(Copy, Clone)]
pub struct EdgeData {
    pub left: usize,
    pub right: usize,
}

#[inline(always)]
fn sort_vertices(&(a, b, c): &(usize, usize, usize)) -> (usize, usize, usize) {
         if a < b && a < c    { (a, b, c) }
    else if b < a && b < c    { (b, c, a) }
    else /* c < a && c < b */ { (c, a, b) }
}

impl Mesh {
    /// Creates a new mesh with the given vertices and faces. The method assumes correctness of the
    /// arguments. In particular, all indices in `faces` must be between 0 and
    /// `positions.len() - 1`, and they must be given in counter-clockwise order.
    pub fn new(positions: Vec<Vector3D>, faces: Vec<(usize, usize, usize)>) -> Mesh {
        
        let mut neighbours = vec![HashMap::new(); positions.len()];
        
        for &(i1, i2, i3) in faces.iter() {
            for &(a, b, c) in [(i1, i2, i3), (i2, i3, i1), (i3, i1, i2)].iter() {
                let mut map = neighbours.get_mut(a).unwrap();
                match map.entry(b) {
                    Entry::Vacant(entry)   => { entry.insert(EdgeData { left: c, right: !0 });   }
                    Entry::Occupied(entry) => { entry.into_mut().left = c; }
                }
                match map.entry(c) {
                    Entry::Vacant(entry)   => { entry.insert(EdgeData { left: !0, right: b });   }
                    Entry::Occupied(entry) => { entry.into_mut().right = b; }
                }
            }
        }
        
        let faces = faces.iter().map(sort_vertices).collect();
        
        Mesh {
            positions: positions,
            faces: faces,
            neighbours: neighbours,
        }
    }
    
    pub fn extrude(&mut self, ring: &Ring, dir: Vector3D) {
        extrusion::extrude(self, ring, dir);
    }
    
    pub fn get_inside_vertices(&self, ring: &Ring) -> Vec<usize> {
        let mut result = vec![];
        
        let start = ring.get_random_vertex();
        let mut current_edge = (start, ring.next(start));
        let mut ring = ring.clone();
        
        while ring.len() > 2 {
            let left = self.neighbours[current_edge.0][&current_edge.1].left;
            
            if !ring.contains_vertex(left) {
                result.push(left);
                ring.insert_edge(current_edge.0, left);
                current_edge = (left, current_edge.1);
            }
            else if ring.prev(left) == current_edge.1 {
                ring.remove_vertex(current_edge.1);
                current_edge = (current_edge.0, left)
            }
            else if ring.next(left) == current_edge.0 {
                ring.remove_vertex(current_edge.0);
                current_edge = (left, current_edge.1)
            }
            else {
                // The ring DOES contain this left vertex, but not nearby. We need to split
                // the ring in two and determine the inside of each component separately.
                
                let ring1 = {
                    let mut current_vertex = current_edge.1;
                    let mut aux_vec = vec![current_vertex];
                    while current_vertex != left {
                        current_vertex = ring.next(current_vertex);
                        aux_vec.push(current_vertex)
                    }
                    Ring::new(&aux_vec)
                };
                
                let ring2 = {
                    let mut current_vertex = left;
                    let mut aux_vec = vec![current_vertex];
                    while current_vertex != current_edge.0 {
                        current_vertex = ring.next(current_vertex);
                        aux_vec.push(current_vertex)
                    }
                    Ring::new(&aux_vec)
                };
                
                result.extend(self.get_inside_vertices(&ring1));
                result.extend(self.get_inside_vertices(&ring2));
                
                break;
            }
            
        }
        
        result
    }
}

// /// Keeps the necessary data to perform extrusion on a mesh with a given ring
// struct ExtrusionData<'a> {
//     mesh: &'a mut Mesh,
//     ring: Ring,
//     dir: Vector3D,
//     inside_vertices: Vec<usize>,
//     extrud_vertices: HashMap<usize, usize>,
//     middle_vertices: HashMap<usize, usize>,
// }


// impl<'a> ExtrusionData<'a> {
    
//     fn perform(mesh: &mut Mesh, ring: Ring, dir: Vector3D) {
//         let data = ExtrusionData {
//             inside_vertices: mesh.get_inside_vertices(&ring),
//             extrud_vertices: mesh.create_extrud_vertices(&ring, dir),
//             middle_vertices: mesh.create_middle_vertices(&ring, dir),
//             ring: ring,
//             mesh: mesh,
//             dir: dir,
//         };
        
//         data.extrude();
//     }
    
//     fn new(mesh: &mut Mesh, ring: Ring, dir: Vector3D) -> Self {
//         let inside_vertices = mesh.get_inside_vertices(ring);
        
//         let extrud_vertices = Self::create_extrud_vertices(mesh, &ring, dir);
//         let middle_vertices = Self::create_middle_vertices(mesh, &ring, dir);
        
//         ExtrusionData {
            
//         }
//     }
    
//     fn extrude(self) {
        
//         // Translate inside vertices in the direction of the extrusion
//         self.mesh.translate_inside(&self.ring, self.dir);
        
//         // Create new vertices
//         let extrud_vertices = self.mesh.create_extrud_vertices(&self.ring, self.dir);
//         let middle_vertices = self.mesh.create_middle_vertices(&self.ring, self.dir);
        
//         // Change the faces on the inside of the ring that contain at least one vertex on the ring
//         // so that they instead use the new vertices on top the old ones
//         self.mesh.rewire_faces_inside_ring(&self.ring, &extrud_vertices);
//         self.mesh.rewire_ring_edges(&self.ring, &extrud_vertices, &middle_vertices);
        
//     }
    
//     fn translate_inside(&mut self, ring: &Ring, dir: Vector3D) {
//         for idx in self.get_inside_vertices(ring) {
//             self.positions[idx] += dir;
//         }
//     }
    
//     fn create_extrud_vertices(mesh: &mut Mesh, ring: &Ring, dir: Vector3D) -> HashMap<usize, usize> {
//         let mut new_vertices = HashMap::new();
        
//         for v in ring.vertices() {
//             let new_pos = self.positions[v] + dir;
//             let new_vertex_idx = self.positions.len();
            
//             self.positions.push(new_pos);
//             new_vertices.insert(v, new_vertex_idx);
//         }
        
//         new_vertices
//     }
    
//     fn create_middle_vertices(mesh: &mut Mesh, ring: &Ring, dir: Vector3D) -> HashMap<usize, usize> {
//         let mut new_vertices = HashMap::new();
        
//         for (v1, v2) in ring.edges() {
//             let new_pos = (self.positions[v1] + self.positions[v2]) / 2.0 + dir / 2.0;
//             let new_vertex_idx = self.positions.len();
            
//             self.positions.push(new_pos);
//             new_vertices.insert(v1, new_vertex_idx);
//         }
        
//         new_vertices
//     }
    
//     fn rewire_faces_inside_ring(&mut self, ring: &Ring, extrud_vertices: &HashMap<usize, usize>) {
        
//         for &mut (ref mut a, ref mut b, ref mut c) in self.faces.iter_mut() {
            
//             let raised_a = extrud_vertices.contains_key(a);
//             let raised_b = extrud_vertices.contains_key(b);
//             let raised_c = extrud_vertices.contains_key(c);
            
//             let raised_ab = ring.has_edge(*a, *b);
//             let raised_bc = ring.has_edge(*b, *c);
//             let raised_ca = ring.has_edge(*c, *a);
            
//             if raised_a {
//                 *a = extrud_vertices[&*a];
//             }
            
//             if raised_b {
//                 *b = extrud_vertices[&*b];
//             }
            
//             if raised_c {
//                 *c = extrud_vertices[&*c];
//             }
//         }
        
        
//     }
    
//     fn rewire_ring_edges(&mut self,
//                          ring: &Ring,
//                          extrud_vertices: &HashMap<usize, usize>,
//                          middle_vertices: &HashMap<usize, usize>) {
        
//         for (v1, v2) in ring.edges() {
//             let old_left = self.neighbours[v1][&v2].left;
//             let middle = middle_vertices[&v1];
            
//             self.neighbours[v1].get_mut(&v2).unwrap().left = middle;
//             self.neighbours[v2].get_mut(&v1).unwrap().right = middle;
            
//             let extruded_v1 = extrud_vertices[&v1];
//             let extruded_v2 = extrud_vertices[&v2];
            
//             self.neighbours[extruded_v1].insert(extruded_v2, EdgeData { left: old_left, right: middle });
//             self.neighbours[extruded_v2].insert(extruded_v1, EdgeData { left: middle, right: old_left });
            
//             self.faces.push(sort_vertices(&(v1, v2, middle)));
//             self.faces.push(sort_vertices(&(v2, extruded_v2, middle)));
//             self.faces.push(sort_vertices(&(extruded_v2, extruded_v1, middle)));
//             self.faces.push(sort_vertices(&(extruded_v1, v1, middle)));
            
//             let middle_prev = middle_vertices[&ring.prev(v1)];
//             let middle_next = middle_vertices[&ring.next(v2)];
            
//             self.neighbours[v1].insert(middle, EdgeData { left: extruded_v1, right: v2} );
//             self.neighbours[middle].insert(v1, EdgeData { left: v2, right: extruded_v1} );
            
//             self.neighbours[v2].insert(middle, EdgeData { left: v1, right: extruded_v2} );
//             self.neighbours[middle].insert(v2, EdgeData { left: extruded_v2, right: v1} );
            
//             self.neighbours[extruded_v2].insert(middle, EdgeData { left: v2, right: extruded_v1} );
//             self.neighbours[middle].insert(extruded_v2, EdgeData { left: extruded_v1, right: v2} );
            
//             self.neighbours[extruded_v1].insert(middle, EdgeData { left: extruded_v2, right: v1} );
//             self.neighbours[middle].insert(extruded_v1, EdgeData { left: v1, right: extruded_v2} );
            
//             self.neighbours[v1].insert(extruded_v1, EdgeData { left: middle_prev, right: middle} );
//             self.neighbours[extruded_v1].insert(v1, EdgeData { left: middle, right: middle_prev} );
            
//             self.neighbours[v2].insert(extruded_v2, EdgeData { left: middle, right: middle_next} );
//             self.neighbours[extruded_v2].insert(v2, EdgeData { left: middle_next, right: middle} );
//         }
        
//     }
    
// }


