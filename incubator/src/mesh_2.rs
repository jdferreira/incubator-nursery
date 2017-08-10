#![allow(dead_code)]

use std::collections::HashMap;
use std::collections::hash_map::{Entry};

use vector::Vector3D;

/// Represents a closed ring of vertices in a `Mesh`
#[derive(Clone,PartialEq,Eq,Debug)]
pub struct Ring {
    inner: HashMap<usize, RingData>,
}

#[derive(Clone,PartialEq,Eq,Debug)]
struct RingData {
    prev: usize,
    next: usize,
}

impl Ring {
    
    /// Creates a `Ring` from a slice.
    pub fn new(slice: &[usize]) -> Ring {
        slice.into()
    }
    
    /// Determines the size of this ring
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    
    /// Gets a random vertex from the ring
    pub fn get_random_vertex(&self) -> usize {
        *self.inner.keys().take(1).next().unwrap()
    }
    
    /// Given a vertex in the ring, returns the next vertex
    pub fn next(&self, vertex: usize) -> usize {
        self.inner[&vertex].next
    }
    
    /// Given a vertex in the ring, returns the next vertex
    pub fn prev(&self, vertex: usize) -> usize {
        self.inner[&vertex].prev
    }
    
    pub fn contains_vertex(&self, v: usize) -> bool {
        self.inner.contains_key(&v)
    }
    
    pub fn insert_edge(&mut self, src: usize, dest: usize) {
        assert!(self.inner.contains_key(&src));
        assert!(!self.inner.contains_key(&dest));
        
        let next = self.next(src);
        
        self.inner.get_mut(&src).unwrap().next = dest;
        self.inner.get_mut(&next).unwrap().prev = dest;
        
        self.inner.insert(dest, RingData { prev: src, next: next });
    }
    
    pub fn remove_vertex(&mut self, vertex: usize) {
        assert!(self.len() >= 3);
        assert!(self.inner.contains_key(&vertex));
        
        let prev = self.prev(vertex);
        let next = self.next(vertex);
        
        self.inner.get_mut(&prev).unwrap().next = next;
        self.inner.get_mut(&next).unwrap().prev = prev;
        
        self.inner.remove(&vertex);
    }
    
    pub fn vertices<'a>(&'a self) -> impl Iterator<Item=usize> + 'a {
        self.inner.keys().map(|&x| x)
    }

    pub fn edges<'a>(&'a self) -> impl Iterator<Item=(usize, usize)> + 'a {
        self.inner.iter().map(|(&x, &RingData { prev: _, next })| (x, next))
    }
    
    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        match self.inner.get(&a) {
            None => false,
            Some(&ref data) => data.next == b || data.prev == b
        }
    }
}

/// Creates a vector containing the vertices of the ring, starting with a random vertex. The
/// vertex next to the first is placed in the second position in the vector, the vertex next to
/// the second in the third position, and so on.
impl <'a> From<&'a Ring> for Vec<usize> {
    fn from(ring: &Ring) -> Vec<usize> {
        let sentinel = ring.get_random_vertex();
        
        let mut result = Vec::with_capacity(ring.len());
        let mut current = sentinel;
        while ring.next(current) != sentinel {
            result.push(current);
            current = ring.next(current);
        }
        
        result.push(current);
        result
    }
}

/// Creates a Ring from a slice of `usize` integers, where each integer in the slice
/// is mapped into the next one, and the last integer is mapped to the first one.
//
// Assumes the given slice has at least two vertices. Otherwise, this conversion panics.
impl <'a> From<&'a [usize]> for Ring {
    fn from(slice: &[usize]) -> Ring {
        assert!(slice.len() >= 2);
        
        // Make a vector that contains all the 
        let n = slice.len();
        let window = {
            let mut window = vec![];
            window.push(slice[n - 1]);
            window.extend(slice);
            window.push(slice[0]);
            window
        };
        let prev = window[0 .. n  ].iter();
        let this = window[1 .. n+1].iter();
        let next = window[2 .. n+2].iter();
        
        // Zip through the given slice by forming pairs of adjacent vertices in the slice
        Ring {
            inner:
                izip!(prev, this, next)
                .map(|(&prev, &this, &next)| (this, RingData { prev: prev, next: next }))
                .collect()
        }
    }
}

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
    pub neighbours: Vec<HashMap<usize, (usize, usize)>>,
    
}

#[inline(always)]
fn sort_vertices(&(a, b, c): &(usize, usize, usize)) -> (usize, usize, usize) {
    if      a < b && a < c    { (a, b, c) }
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
                    Entry::Vacant(entry)   => { entry.insert((c, 0)); },
                    Entry::Occupied(entry) => { entry.into_mut().0 = c; },
                }
                match map.entry(c) {
                    Entry::Vacant(entry)   => { entry.insert((0, b)); },
                    Entry::Occupied(entry) => { entry.into_mut().1 = b; },
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
    
    #[allow(unused_variables)]
    pub fn extrude(&mut self, ring: &Ring, dir: Vector3D) {
        // Translate inside vertices in the direction of the extrusion
        self.translate_inside(ring, dir);
        
        // Create new vertices
        let extrud_vertices = self.create_extrud_vertices(ring, dir);
        let middle_vertices = self.create_middle_vertices(ring, dir);
        
        // Change the faces on the inside of the ring that contain at least one vertex on the ring
        // so that they instead use the new vertices on top the old ones
        self.rewire_faces_inside_ring(ring, &extrud_vertices);
        self.rewire_ring_edges(ring, &extrud_vertices, &middle_vertices);
        
    }
    
    fn translate_inside(&mut self, ring: &Ring, dir: Vector3D) {
        for idx in self.get_inside_vertices(ring) {
            self.positions[idx] += dir;
        }
    }
    
    fn create_extrud_vertices(&mut self, ring: &Ring, dir: Vector3D) -> HashMap<usize, usize> {
        let mut new_vertices = HashMap::new();
        
        for v in ring.vertices() {
            let new_pos = self.positions[v] + dir;
            let new_vertex_idx = self.positions.len();
            
            self.positions.push(new_pos);
            new_vertices.insert(v, new_vertex_idx);
        }
        
        new_vertices
    }
    
    fn create_middle_vertices(&mut self, ring: &Ring, dir: Vector3D) -> HashMap<usize, usize> {
        let mut new_vertices = HashMap::new();
        
        for (v1, v2) in ring.edges() {
            let new_pos = (self.positions[v1] + self.positions[v2]) / 2.0 + dir / 2.0;
            let new_vertex_idx = self.positions.len();
            
            self.positions.push(new_pos);
            new_vertices.insert(v1, new_vertex_idx);
        }
        
        new_vertices
    }
    
    fn rewire_faces_inside_ring(&mut self, ring: &Ring, extrud_vertices: &HashMap<usize, usize>) {
        
        for &mut (ref mut a, ref mut b, ref mut c) in self.faces.iter_mut() {
            
            let raised_a = extrud_vertices.contains_key(a);
            let raised_b = extrud_vertices.contains_key(b);
            let raised_c = extrud_vertices.contains_key(c);
            
            let raised_ab = ring.has_edge(*a, *b);
            let raised_bc = ring.has_edge(*b, *c);
            let raised_ca = ring.has_edge(*c, *a);
            
            if raised_a {
                *a = extrud_vertices[&*a];
            }
            if raised_b {
                *b = extrud_vertices[&*b];
            }
            if raised_c {
                *c = extrud_vertices[&*c];
            }
        }
        
        
    }
    
    fn rewire_ring_edges(&mut self, ring: &Ring,
                         extrud_vertices: &HashMap<usize, usize>,
                         middle_vertices: &HashMap<usize, usize>) {
        
        for (v1, v2) in ring.edges() {
            let old_left = self.neighbours[v1][&v2].0;
            let middle = middle_vertices[&v1];
            
            self.neighbours[v1].get_mut(&v2).unwrap().0 = middle;
            self.neighbours[v2].get_mut(&v1).unwrap().1 = middle;
            
            let extruded_v1 = extrud_vertices[&v1];
            let extruded_v2 = extrud_vertices[&v2];
            
            self.neighbours[extruded_v1].insert(extruded_v2, (old_left, middle));
            self.neighbours[extruded_v2].insert(extruded_v1, (middle, old_left));
            
            self.faces.push(sort_vertices(&(v1, v2, middle)));
            self.faces.push(sort_vertices(&(v2, extruded_v2, middle)));
            self.faces.push(sort_vertices(&(extruded_v2, extruded_v1, middle)));
            self.faces.push(sort_vertices(&(extruded_v1, v1, middle)));
            
            let middle_prev = middle_vertices[&ring.prev(v1)];
            let middle_next = middle_vertices[&ring.next(v2)];
            
            self.neighbours[v1].insert(middle, (extruded_v1, v2));
            self.neighbours[middle].insert(v1, (v2, extruded_v1));
            
            self.neighbours[v2].insert(middle, (v1, extruded_v2));
            self.neighbours[middle].insert(v2, (extruded_v2, v1));
            
            self.neighbours[extruded_v2].insert(middle, (v2, extruded_v1));
            self.neighbours[middle].insert(extruded_v2, (extruded_v1, v2));
            
            self.neighbours[extruded_v1].insert(middle, (extruded_v2, v1));
            self.neighbours[middle].insert(extruded_v1, (v1, extruded_v2));
            
            self.neighbours[v1].insert(extruded_v1, (middle_prev, middle));
            self.neighbours[extruded_v1].insert(v1, (middle, middle_prev));
            
            self.neighbours[v2].insert(extruded_v2, (middle, middle_next));
            self.neighbours[extruded_v2].insert(v2, (middle_next, middle));
        }
        
    }
    
    fn get_inside_vertices(&self, ring: &Ring) -> Vec<usize> {
        let mut result = vec![];
        
        let start = ring.vertices().take(1).next().unwrap();
        let mut current_edge = (start, ring.next(start));
        let mut ring = ring.clone();
        
        while ring.len() > 2 {
            let left = self.neighbours[current_edge.0][&current_edge.1].0;
            
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

#[cfg(test)]
mod tests {
    
    use super::*;
    
    #[test]
    #[should_panic]
    fn test_new_ring_empty() {
        Ring::new(&[]);
    }
    
    #[test]
    #[should_panic]
    fn test_new_ring_size_1() {
        Ring::new(&[1]);
    }

    #[test]
    fn test_slice_to_ring_small() {
        let ring = Ring::new(&[1, 2]);

        assert_eq!(ring.next(1), 2);
        assert_eq!(ring.next(2), 1);

        assert_eq!(ring.prev(1), 2);
        assert_eq!(ring.prev(2), 1);
    }
    
    #[test]
    fn test_slice_to_ring() {
        let ring = Ring::new(&[1, 2, 3, 4, 5]);

        assert_eq!(ring.next(1), 2);
        assert_eq!(ring.next(2), 3);
        assert_eq!(ring.next(3), 4);
        assert_eq!(ring.next(4), 5);
        assert_eq!(ring.next(5), 1);

        assert_eq!(ring.prev(1), 5);
        assert_eq!(ring.prev(2), 1);
        assert_eq!(ring.prev(3), 2);
        assert_eq!(ring.prev(4), 3);
        assert_eq!(ring.prev(5), 4);
    }
    
    #[test]
    #[should_panic]
    fn test_ring_next_on_non_existing_vertex() {
        let ring = Ring::new(&[1, 2, 3, 4, 5]);
        ring.next(10);
    }
    
    fn shift_slice(slice: &mut [usize], start: usize) {
        
        let shift = match slice.iter().position(|&el| el == start) {
            Some(idx) => idx,
            None => return,
        };
        
        let new_vec = {
            let iter = slice.iter().skip(shift).chain(slice.iter().take(shift));
            iter.cloned().collect::<Vec<_>>()
        };
        
        for (idx, &i) in new_vec.iter().enumerate() {
            slice[idx] = i;
        }
    }
    
    fn test_ring_to_slice_unit(slice: &[usize]) {
        let ring = Ring::new(slice);
        let mut vec = Vec::from(&ring);
        shift_slice(&mut vec, slice[0]);
        assert_eq!(vec, slice);
    }
    
    #[test]
    fn test_ring_to_slice() {
        test_ring_to_slice_unit(&[1, 2]);
        test_ring_to_slice_unit(&[1, 2, 3]);
        test_ring_to_slice_unit(&[1, 2, 3, 4]);
        test_ring_to_slice_unit(&[1, 2, 3, 4, 5]);
    }
    
    #[test]
    fn test_get_inside_vertices() {
        let positions = vec![Vector3D(0.0, 0.0, 0.0); 11];
        let faces = vec![
            (0,1,3),
            (1,2,4),
            (1,4,3),
            (3,7,6),
            (3,4,7),
            (4,8,7),
            (4,5,8),
            (6,7,9),
            (9,7,10),
            (7,8,10),
        ];
        let mesh = Mesh::new(positions, faces);
        let ring = Ring::new(&[0, 1, 2, 4, 5, 8, 10, 9, 6, 3]);
        let mut result = mesh.get_inside_vertices(&ring);
        result.sort();

        assert_eq!(result, vec![7]);
        
        let positions = vec![Vector3D(0.0, 0.0, 0.0); 16];
        let faces = vec![
            (0,1,4),
            (1,5,4),
            (1,2,5),
            (2,6,5),
            (2,3,7),
            (2,7,6),
            (4,5,9),
            (4,9,8),
            (5,6,10),
            (5,10,9),
            (6,7,10),
            (7,11,10),
            (8,9,13),
            (8,13,12),
            (9,10,13),
            (10,14,13),
            (10,11,14),
            (11,15,14),
        ];
        let mesh = Mesh::new(positions, faces);
        let ring = Ring::new(&[0,1,2,3,7,11,15,14,13,12,8,4]);
        let mut result = mesh.get_inside_vertices(&ring);
        result.sort();

        assert_eq!(result, vec![5,6,9,10]);
        
        let positions = vec![Vector3D(0.0, 0.0, 0.0); 17];
        let faces = vec![
            (0,1,4),
            (1,5,4),
            (1,2,5),
            (2,6,5),
            (2,3,7),
            (2,7,6),
            (4,5,9),
            (4,9,8),
            (5,6,16),
            (6,10,16),
            (10,9,16),
            (9,5,16),
            (6,7,10),
            (7,11,10),
            (8,9,13),
            (8,13,12),
            (9,10,13),
            (10,14,13),
            (10,11,14),
            (11,15,14),
        ];
        let mesh = Mesh::new(positions, faces);
        let ring = Ring::new(&[0,1,2,3,7,11,15,14,13,12,8,4]);
        let mut result = mesh.get_inside_vertices(&ring);
        result.sort();

        assert_eq!(result, vec![5,6,9,10,16]);
    }
    
    #[test]
    fn test_extrude() {
        let positions = vec![
            Vector3D(0.0, 0.0, 0.0),
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
            Vector3D(0.0, 0.0, -9.0),
        ];
        let faces = vec![
            ( 0,  1,  4),
            ( 1,  2,  5),
            ( 2,  3,  5),
            ( 1,  6,  4),
            ( 1,  5,  6),
            ( 5,  3, 10),
            ( 6,  5, 10),
            ( 4,  6,  7),
            ( 7,  6,  8),
            ( 8,  6,  9),
            ( 9,  6, 10),
            ( 0, 11,  1),
            ( 1, 11,  2),
            ( 2, 11,  3),
            ( 3, 11, 10),
            (10, 11,  9),
            ( 9, 11,  8),
            ( 8, 11,  7),
            ( 7, 11,  4),
            ( 4, 11,  0),
        ];
        let mut mesh = Mesh::new(positions, faces);
        let ring = Ring::new(&[0,1,2,3,10, 9, 8, 7, 4]);
        
        mesh.extrude(&ring, Vector3D(0.0, 0.0, 1.0));
    }
}
