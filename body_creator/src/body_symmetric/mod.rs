pub mod consts;

use std::f64;
use std::collections::{HashSet,HashMap};

use vector::Vector3D;
use field::Field;
use genome::Genome;
use mesh::Mesh;

use utils::math::convert_to_gaussian;


pub struct Body {
    /// The mesh that defines the body
    pub mesh: Mesh,
    
    // Private fields
    field: Field,
    mirror: HashMap<usize, usize>,
    velocities: Vec<Vector3D>,
}


impl Body {
    
    pub fn new() -> Body {
        // Start by defining the position of the first 12 vertices
        let phi = (1.0 + 5f64.sqrt()) / 2.0; // TODO: constexpr; alas!, Rust doesn't have them yet
        let Vector3D(x0, y0, z0) = Vector3D(1.0, phi, 0.0).normalize() * consts::RADIUS;
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
        
        // Find mirror vertices
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
        
        // Now let's divide the faces and create new vertices. We will need a BodyBuilder to pass
        // around the information (it was either that or a big function signature with several
        // arguments)
        let mut bb = BodyBuilder {
            positions: positions,
            faces: faces,
            mirror: mirror,
            radius: consts::RADIUS,
            divisions: consts::EMBRYO_DIVISIONS,
            mid_cache: HashMap::new(),
        };
        
        // Divide this bad boy
        bb.divide();
        
        // Finally, return a body that contains the mesh corresponding to the vertices and faces
        // calculated so far and the mirror information. The body needs also the field and the
        // velocities associated with each vertex
        Body {
            field: Field::new(consts::FIELD_DIVISION),
            velocities: vec![Vector3D(0.0, 0.0, 0.0); bb.positions.len()],
            mesh: Mesh::new(bb.positions, bb.faces),
            mirror: bb.mirror,
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
        for index in 0..self.mesh.positions.len() {
            self.update_velocities(index);
        }
        
        self.move_vertices();
    }
    
    fn move_vertices(&mut self) {
        for index in 0..self.mesh.positions.len() {
            let pos = self.mesh.positions[index];
            let vel = self.velocities[index];
            self.mesh.positions[index] = pos + vel * consts::TIME_STEP;
        }
    }
    
    fn update_velocities(&mut self, index: usize) {
        let pos = self.mesh.positions[index];
        
        // Where is the field sending me to?
        let field_acc = consts::FIELD_STRENGTH * self.field.get_nabla(pos);
        
        // This vertex wants to run away from other vertices.
        let mut proximity_acc = Vector3D(0.0, 0.0, 0.0);
        for other in 0..self.mesh.positions.len() {
            if other != index {
                // This is a vertex other than the one whose velocity is
                // being computed, so we must react to its presence
                let delta = pos - self.mesh.positions[other];
                let dist = delta.norm();
                let factor = consts::PROXIMITY_STRENGTH / (dist * dist * dist);
                
                proximity_acc = proximity_acc + factor * delta;
            }
        }
        
        // Another force at play is the attraction to the neighbours.
        let mut pivot = Vector3D(0.0, 0.0, 0.0);
        for (&n, _) in self.mesh.neighbours[index].iter() {
            pivot = pivot + self.mesh.positions[n];
        }
        let n_neighbours = self.mesh.neighbours[index].len() as f64;
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
    
}


struct BodyBuilder {
    positions: Vec<Vector3D>,
    faces: HashSet<(usize, usize, usize)>,
    mirror: HashMap<usize, usize>,
    radius: f64,
    divisions: u32,
    mid_cache: HashMap<(usize, usize), usize>,
}

impl BodyBuilder {
    fn divide(&mut self) {
        for _ in 0..self.divisions {
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
}
