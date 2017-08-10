pub mod consts;

use std::f64;
use std::collections::{HashSet,HashMap};
use std::collections::hash_map::Entry;

use vector::Vector3D;
use field::Field;
use genome::Genome;
use mesh::Mesh;
use face_history::{FaceHistory,FaceStep};

use utils::math::convert_to_gaussian;


#[inline(always)]
fn sort_vertices(a: usize, b: usize, c:usize) -> (usize, usize, usize) {
    if      a < b && a < c    { (a, b, c) }
    else if b < a && b < c    { (b, c, a) }
    else /* c < a && c < b */ { (c, a, b) }
}

pub struct Body {
    /// The positions of each vertex. Vertices are identified by their index in this vector.
    /// Vertices that disappear (after extrusion) are relocated to Vector3D(0.0, 0.0, 0.0).
    pub positions: Vec<Vector3D>,
    
    /// The faces in the mesh
    pub faces: HashMap<(usize, usize, usize), FaceHistory>,
    
    /// Each element contains the neighbours of the vertex. The key in the map is the neighbour
    /// (allows for quick checking of neighbours) and the value is a tuple, where the first value
    /// is the left vertex of this edge and the second value is the right vertex
    pub neighbours: Vec<HashMap<usize, (usize, usize)>>,
    
    
    // Private fields
    
    /// The field that defines the force moving the vertices
    field: Field,
    
    /// Since the body is mirrored in the x-plane, all vertices have a mirror image. This map
    /// associates each vertex with its mirror image. Verttices on the x-plane are their own mirror
    /// images.
    mirror: HashMap<usize, usize>,
    
    /// The velocities that each vertex has in a given time instance. The velocities are updated
    /// based on the environment (field, nieghbour vertices, etc.) and they, in turn, are used
    /// to update the positions of the next iteration
    velocities: Vec<Vector3D>,
    
    /// The current iteration
    iteration: usize,
}

// TODO: Each vertex will also have to have some color information to paint the faces. That
// information **may** change with time (for extra flexibility) but I still need to think about
// that. In that case, color may be based on a field of its own.

impl Body {
    
    pub fn new() -> Body {
        let mut bb = BodyBuilder::new(consts::RADIUS);
        bb.divide(consts::EMBRYO_DIVISIONS);
        
        // Finally, return a body that contains the mesh corresponding to the vertices and faces
        // calculated so far and the mirror information. The body needs also the field and the
        // velocities associated with each vertex
        let faces = bb.faces();
        let neighbours = bb.neighbours();
        let n = bb.positions.len();
        
        Body {
            positions: bb.positions,
            faces: faces,
            neighbours: neighbours,
            field: Field::new(consts::FIELD_DIVISION),
            mirror: bb.mirror,
            velocities: vec![Vector3D(0.0, 0.0, 0.0); n],
            iteration: 0
        }
    }
    
    /// Uses a genome to create the field of values
    pub fn attach_genome(&mut self, genome: &mut Genome) {
        
        for _ in 0..consts::NUM_PARTICLES {
            // Create a new particle and add it to the field. The parameters of the particle
            // are retrieved from the genome.
            let x = genome.get_float();
            let y = genome.get_float();
            let z = genome.get_float();
            
            // attraction is a "corrected" gaussian value with mean 0 and standard deviation 1.
            // The correction is made to make sure the value is between -1 and 1
            let attraction = convert_to_gaussian(genome.get_float(), 0.0, 1.0, -1.0, 1.0);
            
            // The spread is also a "corrected" gaussian with mean 0 and standard deviation 1,
            // with a minimum value of 1. We scale the value so that we have a value not
            // restricted to be greater than 1 but greater than the minimum spread.
            let spread = convert_to_gaussian(genome.get_float(), 0.0, 1.0, 1.0, f64::INFINITY) *
                         consts::MIN_SPREAD;
            
            self.field.add_particle(Vector3D(x, y, z), attraction, spread);
        }
        
        // // We also add a central repulsive particle so that the body contains a
        // // central mass. This seems to be good to create nice looking bodies
        // let attraction = convert_to_gaussian(genome.get_float(), 1.0, 1.0, -1.0, 0.0);
        // let spread = convert_to_gaussian(genome.get_float(), 0.0, 1.0, 1.0, f64::INFINITY) *
        //              consts::CENTRAL_SPREAD;
                     
        // self.add_particle(Vector3D(0.0, 0.0, 0.0), attraction, spread);
    }
    
    pub fn step(&mut self) {
        for index in 0..self.positions.len() {
            self.update_velocities(index);
        }
        
        self.move_vertices();
    }
    
    fn move_vertices(&mut self) {
        for index in 0..self.positions.len() {
            let pos = self.positions[index];
            let vel = self.velocities[index];
            self.positions[index] = pos + vel * consts::TIME_STEP;
        }
        
        // We now need to go into each face and update their history
        for (&(i1, i2, i3), history) in self.faces.iter_mut() {
            // Get the position of the three vertices
            let v1 = self.positions[i1];
            let v2 = self.positions[i2];
            let v3 = self.positions[i3];
            
            // Get the sides of the trinagle
            let a = (v1 - v2).norm();
            let b = (v2 - v3).norm();
            let c = (v3 - v1).norm();
            
            // Calculate the area of the triangle
            let p = (a + b + c ) / 2.0;
            let area = (p * (p - a) * (p - b) * (p - c)).sqrt();
            
            // Calculate the information of the triangle
            let (width, height, mid, corner);
            if a > b && a > c {
                corner = i1;
                width = a;
                height = 2.0 * area / width;
                mid = (c * c - height * height).sqrt();
            }
            else if b > a && b > c {
                corner = i2;
                width = b;
                height = 2.0 * area / width;
                mid = (a * a - height * height).sqrt();
            }
            else {
                corner = i3;
                width = c;
                height = 2.0 * area / width;
                mid = (b * b - height * height).sqrt();
            }
            
            // Add a new step to the history
            history.add_step(self.iteration, corner, width, height, mid);
        }
    }
    
    fn update_velocities(&mut self, index: usize) {
        let pos = self.positions[index];
        
        // Where is the field sending me to?
        let field_acc = consts::FIELD_STRENGTH * self.field.get_nabla(pos);
        
        // This vertex wants to run away from other vertices.
        let mut proximity_acc = Vector3D(0.0, 0.0, 0.0);
        for other in 0..self.positions.len() {
            if other != index {
                // This is a vertex other than the one whose velocity is
                // being computed, so we must react to its presence
                let delta = pos - self.positions[other];
                let dist = delta.norm();
                let factor = consts::PROXIMITY_STRENGTH / (dist * dist * dist);
                
                proximity_acc = proximity_acc + factor * delta;
            }
        }
        
        // Another force at play is the attraction to the neighbours.
        let mut pivot = Vector3D(0.0, 0.0, 0.0);
        for (&n, _) in self.neighbours[index].iter() {
            pivot = pivot + self.positions[n];
        }
        let n_neighbours = self.neighbours[index].len() as f64;
        let neighbour_acc = (pivot / n_neighbours - pos) * consts::NEIGHBOURS_STRENGTH / n_neighbours;
        
        // We also want to be stay clear of symmetry plane
        let mut symmetry_acc = Vector3D(0.0, 0.0, 0.0);
        
        if self.mirror[&index] != index {
            if 0.0 < pos.0 && pos.0 < consts::SYMMETRY_EDGE_THICKNESS {
                symmetry_acc.0 = consts::SYMMETRY_EDGE_THICKNESS - pos.0;
            }
            else if consts::SYMMETRY_EDGE_THICKNESS < pos.0 && pos.0 < 0.0 {
                symmetry_acc.0 = -consts::SYMMETRY_EDGE_THICKNESS - pos.0;
            }
        }
        
        // Now let's merge together all the effects into a single acceleration
        // and update the velocity of the vertex
        let vel = self.velocities[index];
        let acc = field_acc + proximity_acc + neighbour_acc + symmetry_acc;
        
        // We take into account both the drag and the time step
        let mut new_vel = consts::DRAG * (vel + acc * consts::TIME_STEP);
        if self.mirror[&index] == index {
            new_vel.0 = 0.0;
        }
        
        let speed = new_vel.norm();
        if speed > consts::MAX_VELOCITY {
            new_vel = new_vel / speed * consts::MAX_VELOCITY;
        }
        
        self.velocities[index] = new_vel;
    }
    
    fn compute_region(&self, ring: &HashMap<usize, usize>) -> Vec<usize> {
        // Make a copy of the ring so that we can change it
        let mut ring: HashMap<usize, usize> = ring.iter().cloned().collect();
        
        // Get one of the vertices of the ring (anyone can do), along with its previous and next
        // vertices in the ring
        let mut prev = *ring.keys().next().unwrap();
        let mut curr = ring[prev];
        let mut next = ring[curr];
        let result = vec![];
        
        // Process the ring as long as it has more than three vertices As soon as three vertices
        // are reached, no more internal vertices can exist and we are done
        while ring.len() > 3 {
            // Get the left vertice to the edge between `curr` and `next`
            let left = self.neighbours[curr][next].0;
            if !ring.contains_key(left) {
                // The left vertex is inside the ring, which means that it will be part of the
                // result of this function
                result.push(left);
                
                // If this left vertex is the same as the left vertex between `prev` and `curr`, we
                // can route the ring through `prev` -- `left` -- `next`
                if left == self.neighbours[prev][curr].0 &&  {
                    ring[prev] = left;
                    ring[left] = next;
                    curr = left;
                }
                // Otherwise, gobble this face and update the `curr` and `prev`
                else {
                    ring[curr] = left;
                    ring[left] = next;
                    curr = left;
                    prev = curr;
                }
            }
            else {
                // If the left vertex is the ring vertex after `next`, we can gobble on face
                if ring[next] == left {
                    ring[curr] = left;
                    next = left;
                }
                // The same can be done if the left vertex is the `prev`
                else if left == prev {
                    ring[prev] = next;
                    curr = next;
                    next = ring[next];
                }
                // Otherwise, just move along the ring until we find a place that we can process
                else {
                    prev = curr;
                    curr = next;
                    next = ring[next];
                }
            }
        }
        
        result
    }
    
    /// Given a ring of vertices, extrude them in a certain direction so that new faces emerge.
    ///
    /// The ring must be given such that the left neighbour is inside the region defined by the
    /// ring.                                ^^^^ TODO: Fix me if not correct
    pub fn extrude(&mut self, ring: &HashMap<usize, usize>, dir: Vector3D) {
        
        // For each edge on the ring, we place a vertex on its middle, shifted by `dir/ 2.0`
        let mut middle_vertices = HashMap::new();
        for (&curr, &next) in ring.iter() {
            self.positions.push((self.positions[curr] + self.positions[next] + dir) / 2.0);
            
            // The velocity of the new vertex will be the average of the velocities of the two
            // old vertices
            self.velocities.push((self.velocities[curr] + self.velocities[next]) / 2.0);
            
            // Associate this current vertex to the new vertex
            middle_vertices.insert(curr, self.positions.len() - 1);
        }
        
        // For each of the new edges created above (the edges have not been created yet, but
        // will later), create yet another vertex in its middle point, shifted by `dir / 2.0`
        let mut top_vertices = HashMap::new();
        for (&curr, &next) in ring.iter() {
            let v1 = middle_vertices[curr];
            let v2 = middle_vertices[next];
            
            self.positions.push((self.positions[v1] + self.positions[v2] + dir) / 2.0);
            self.velocities.push((self.velocities[v1] + self.velocities[v2]) / 2.0);
            top_vertices.insert(next, self.positions.len() - 1);
        }
        
        // We also need to shift the inner points of the region by `dir`
        let region = self.compute_region(ring);
        for &index in region.iter() {
            self.positions[index] = self.positions[index] + dir;
        }
        
        // All the faces that have at least one vertex in the ring and all the other vertices
        // in the ring or within the region must be updated
        for &(v1, v2, v3) in self.faces.keys() {
            let in_ring_1 = ring.contains_key(v1);
            let in_ring_2 = ring.contains_key(v2);
            let in_ring_3 = ring.contains_key(v3);
            
            if !in_ring_1 && !in_ring_2 && !in_ring_3 {
                // No vertex on the ring
                continue;
            }
            
            let in_region_1 = if in_ring_1 { false } else { region.iter().position() }
        }
        
        // Go around the ring and cut and glue the vertices with new edges.
        // This requires that we take care of the faces, including their history.
        for (&curr, &next) in ring.iter() {
            // The faces that go away are the ones around the vertex `next`, according to this:
            // 1. let (prev, top) = (curr, left of the edge (curr, next))
            // 2. remove face (prev, next, top)
            // 3. let (prev, top) = (top, left of the edge (top, next))
            // 4.
            //   a. if top == next of `next`, break;
            //   b. otherwise repeat from 2
            // For each face removed, we need to add the corresponding one in the extruded region,
            // where vertices in the ring are replaced by their counterpart (with the
            // `top_vertices` map). Notice that it can happen that the same face is processed twice
            // with this; in that case, do not add the face again.
            
            
            let (mut prev, mut top) = (curr, self.neighbours[curr][next].0);
            loop {
                self.faces[sort_vertices(prev, next, top)].remove();
                
                // TODO: What about the left and right of the edges that are changed?????
                
                if prev == curr { prev = top_vertices[&curr]; }
                let face = sort_vertices(prev, ring[&next], top);
                if !self.faces.contains_key(&face) {
                    self.faces[face] = FaceHistory::new_on_iteration(self.iteration);
                }
                
                prev = top;
                top = self.neighbours[top][next].0
                if ring[&next] == top { break; }
            }
            
            // We have now handled the faces that will be removed and their extruded counterparts.
            // We still need to take care of the faces on the edge of the extrusion
            
        }
        
        
        
        
        
        // For each vertex in the ring, make a new one just beside it (shifted by `dir`). Also, for
        // each of these vertices, make a new one between it and its next one on the ring, also
        // shifted, but now by `dir / 2.0`
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

/// This struct is used to initialize the body mesh, and takes care of dividing the faces of the
/// initial icosahedron a certain amount of times.
struct BodyBuilder {
    positions: Vec<Vector3D>,
    faces: HashSet<(usize, usize, usize)>,
    mirror: HashMap<usize, usize>,
    radius: f64,
    mid_cache: HashMap<(usize, usize), usize>,
}

impl BodyBuilder {
    
    fn new(radius: f64) -> BodyBuilder {
        // Start by defining the position of the first 12 vertices of the icosahedron
        let phi = (1.0 + 5f64.sqrt()) / 2.0; // TODO: constexpr; alas!, Rust doesn't have them yet
        
        let Vector3D(x0, y0, z0) = Vector3D(1.0, phi, 0.0).normalize() * radius;
        let positions = vec![
            Vector3D( x0,  y0,  z0), // 0
            Vector3D( x0, -y0,  z0), // 1
            Vector3D(-x0,  y0,  z0), // 2
            Vector3D(-x0, -y0,  z0), // 3
            Vector3D( z0,  x0,  y0), // 4
            Vector3D( z0,  x0, -y0), // 5
            Vector3D( z0, -x0,  y0), // 6
            Vector3D( z0, -x0, -y0), // 7
            Vector3D( y0,  z0,  x0), // 8
            Vector3D(-y0,  z0,  x0), // 9
            Vector3D( y0,  z0, -x0), // 10
            Vector3D(-y0,  z0, -x0)  // 11
        ];
        
        // Create the faces that go with them
        let faces = [
            (0, 2, 4), (0, 4, 8), (0, 8, 10), (0, 10, 5), (0, 5, 2),  // 5 faces around point 0
            (2, 9, 4), (4, 6, 8), (8, 1, 10), (10, 7, 5), (5, 11, 2), // 5 adjacent faces
            (3, 1, 6), (3, 6, 9), (3, 9, 11), (3, 11, 7), (3, 7, 1),  // 5 faces around point 3
            (1, 8, 6), (6, 4, 9), (9, 2, 11), (11, 5, 7), (7, 10, 1)  // 5 adjacent faces
        ].iter().cloned().collect();
        
        // Calculate mirror vertices
        let mut mirror = HashMap::new();
        
        // Points that are mirrored through the plane `x = 0`
        for &(i1, i2) in [(0, 2), (1, 3), (8, 9), (10, 11)].iter() {
            mirror.insert(i1, i2);
            mirror.insert(i2, i1);
        }
        // Points that are on the plane `x = 0`
        for &i in [4, 5, 6, 7].iter() {
            mirror.insert(i, i);
        }
        
        BodyBuilder {
            positions: positions,
            faces: faces,
            mirror: mirror,
            radius: radius,
            mid_cache: HashMap::new(),
        }
    }
    
    fn divide(&mut self, divisions: u32) {
        for _ in 0..divisions {
            // Copy the set of the current faces to a temporary location and clear it
            let old_faces = self.faces.clone();
            self.faces.clear();
            
            // Divide each of the old faces into four new ones
            for &(i1, i2, i3) in old_faces.iter() {
                let i4 = self.get_middle_vertex(i1, i2);
                let i5 = self.get_middle_vertex(i2, i3);
                let i6 = self.get_middle_vertex(i3, i1);
                self.faces.insert((i1, i4, i6));
                self.faces.insert((i2, i5, i4));
                self.faces.insert((i3, i6, i5));
                self.faces.insert((i4, i5, i6));
            }
        }
    }

    fn get_middle_vertex(&mut self, i: usize, j: usize) -> usize {
        let key = if i < j { (i, j) } else { (j, i) };
        
        // If the pair has already been used to calculate a new vertex, return
        // that vertex
        if let Some(&result) = self.mid_cache.get(&key) {
            return result;
        }
        
        // Get the middle point between the two vertices
        let v1 = self.positions[i];
        let v2 = self.positions[j];
        let v3 = ((v1 + v2) / 2.0).normalize() * self.radius;
        // We use the normalize method to extend the vertex so that it lies on
        // the surface of the sphere.
        
        // Add the new vertex to the list of vertices and hold its index
        self.positions.push(v3);
        let index = self.positions.len() - 1;
        
        let mi = self.mirror[&i];
        let mj = self.mirror[&j];
        
        if mi == i && mj == j {
            // The two end points are on the symmetry plane, so this one should be as well
            self.mirror.insert(index, index);
        }
        else if mi == j && mj == i {
            // The end points of the edge being splitted are their own mirror images. This means
            // that the resulting vertex must be on the symmetry plane
            self.mirror.insert(index, index);
        }
        else {
            // There is a mirror image. It may have already been created or it will be created
            // later. Since we need the index of the mirror vertex, we update `self.mirror` only
            // when the mirror image already exists. We update both ways. This ensures consistency.
            let key = if mi < mj { (mi, mj) } else { (mj, mi) };
            if let Some(&result) = self.mid_cache.get(&key) {
                self.mirror.insert(index, result);
                self.mirror.insert(result, index);
            }
        }
        
        // Cache the value
        self.mid_cache.insert(key, index);
        
        // And finally return the index of the middle vertex
        index
    }
    
    fn faces(&self) -> HashMap<(usize, usize, usize), FaceHistory> {
        let mut faces = HashMap::new();
        for &(i1, i2, i3) in self.faces.iter() {
            faces.insert(sort_vertices(i1, i2, i3), FaceHistory::new());
        }
        faces
    }
    
    fn neighbours(&self) -> Vec<HashMap<usize, (usize, usize)>> {
        let mut neighbours = vec![HashMap::new(); self.positions.len()];
        
        for &(i1, i2, i3) in self.faces.iter() {
            for &(a, b, c) in [(i1, i2, i3), (i2, i3, i1), (i3, i1, i2)].iter() {
                let map = neighbours.get_mut(a).unwrap();
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
        
        neighbours
    }
}
