use std::collections::HashSet;
use std::collections::HashMap;
use std::rc::Rc;
use std::cell::RefCell;
use std::f64;

use vector::Vector3D;
use field::{Field, Particle};
use genome::Genome;
use utils::math::interpolate;

#[derive(Clone,Debug)]
struct Edge {
    start: usize,
    end: usize,
    up: usize,
    down: usize
}

impl Edge {
    fn new_ref(start: usize, end: usize) -> Rc<RefCell<Edge>> {
        Rc::new(RefCell::new(
            Edge { start: start, end: end, up: 0, down: 0 }
        ))
    }
}

/// Represents half an icosphere, which will serve as the starting point
/// tto construct teh symmetric body
#[derive(Debug)]
struct HalfIcosphere {
    /// The initial radius of the icosphere that will be evolved into the body
    radius: f64,
    
    /// The vertices. usize integers identify each vertex
    vertices: Vec<Vector3D>,
    
    /// A set of vertices that are locked in the symmetry axis x = 0
    locked: HashSet<usize>,
    
    /// The faces in the mesh
    faces: HashSet<(usize, usize, usize)>,
    
    /// A map from pairs of vertices to the edge between them (if any exists)
    pair_to_edge: HashMap<(usize, usize), Rc<RefCell<Edge>>>,
    
    /// A map from vertices to the edges that they are part of
    vertices_to_edges: Vec<Vec<Rc<RefCell<Edge>>>>
}


impl HalfIcosphere {
    
    fn new(divisions: u32, radius: f64) -> HalfIcosphere {
        
        // Start by defining the position of the first 10 vertices
        let phi = (1.0 + 5f64.sqrt()) / 2.0; // FIXME: Should be a constexpr,
                                             // but Rust doesn't have them
        let Vector3D(x0, y0, _) = Vector3D(1.0, phi, 0.0).normalize() * radius;
        let y1 = radius;
        
        let vertices = vec![
            Vector3D(x0,  y0,  0.),
            Vector3D(x0, -y0,  0.),
            Vector3D(0.,  y1,  0.),
            Vector3D(0., -y1,  0.),
            Vector3D(0.,  x0,  y0),
            Vector3D(0.,  x0, -y0),
            Vector3D(0., -x0,  y0),
            Vector3D(0., -x0, -y0),
            Vector3D(y0,  0.,  x0),
            Vector3D(y0,  0., -x0)
        ];
        
        // The locked vertices are also fully specified here. These are the
        // ones that have x = 0
        let locked = [2, 3, 4, 5, 6, 7].iter().cloned().collect();
        let faces = [
            (0, 2, 4), (0, 4, 8), (0, 8, 9), (0, 9, 5), (0, 5, 2), (3, 1, 6),
            (4, 6, 8), (8, 1, 9), (9, 7, 5), (3, 7, 1), (1, 8, 6), (7, 9, 1)
        ].iter().cloned().collect();
        
        // We create the result in a halfstate. We do this becuase it is easier
        // to refer to self.radius ... etc than to pass the values and other
        // objects as arguments to functions
        let mut result = HalfIcosphere {
            radius: radius,
            vertices: vertices,
            locked: locked,
            faces: faces,
            pair_to_edge: HashMap::new(),
            vertices_to_edges: vec![]
        };
        
        // Refine the shape using the number of requested divisions
        result.divide(divisions);
        result
    }
    
    fn divide(&mut self, divisions: u32) {
        let mut vertex_cache = HashMap::new();
        for _ in 0..divisions {
            let old_faces = self.faces.clone();
            self.faces.clear();
            
            for &(i1, i2, i3) in old_faces.iter() {
                let i4 = self.get_middle_vertex(i1, i2, &mut vertex_cache);
                let i5 = self.get_middle_vertex(i2, i3, &mut vertex_cache);
                let i6 = self.get_middle_vertex(i3, i1, &mut vertex_cache);
                self.faces.insert((i1, i4, i6));
                self.faces.insert((i2, i5, i4));
                self.faces.insert((i3, i6, i5));
                self.faces.insert((i4, i5, i6));
            }
        }
        
        self.vertices_to_edges = vec![vec![]; self.vertices.len()];
        
        // We need to clone the faces in order to avoid borrowing mut self
        // while still making changes inside the iterator
        let cloned_faces = self.faces.clone();
        for &(v1, v2, v3) in cloned_faces.iter() {
            for &(i1, i2, i3) in [(v1, v2, v3),
                                  (v2, v3, v1),
                                  (v3, v1, v2)].iter() {
                
                // Note: I need to comment throroughly this section of code.
                // Rc and RefCell and heap and stuff...
                let edge;
                if self.pair_to_edge.contains_key(&(i1, i2)) {
                    edge = self.pair_to_edge.get(&(i1, i2)).unwrap().clone();
                    edge.borrow_mut().down = i3;
                }
                else {
                    edge = Edge::new_ref(i1, i2);
                    edge.borrow_mut().up = i3;
                    self.pair_to_edge.insert((i1, i2), edge.clone());
                    self.pair_to_edge.insert((i2, i1), edge.clone());
                }
                
                self.vertices_to_edges[i1].push(edge.clone());
                self.vertices_to_edges[i2].push(edge.clone());
            }
        }
    }
    
    fn get_middle_vertex(&mut self,
                         i: usize, j: usize,
                         cache: &mut HashMap<(usize, usize), usize>) -> usize {
        let key = if i < j { (i, j) } else { (j, i) };
        
        // If the pair has already been used to calculate a new vertex, return
        // that vertex
        match cache.get(&key) {
            Some(&result) => return result,
            None => ()
        }
        
        // Get the middle point between the two vertices
        let v1 = self.vertices[i];
        let v2 = self.vertices[j];
        let v3 = (0.5 * (v1 + v2)).normalize() * self.radius;
        // We use the normalize method to extend the vertex so that it lies on
        // the surface of the sphere.
        
        // Add the new vertex to the list of vertices and hold its index
        self.vertices.push(v3);
        let index = self.vertices.len() - 1;
        cache.insert(key, index);
        
        // If this vertex is in between locked vertices, it is also locked.
        if self.locked.contains(&i) && self.locked.contains(&j) {
            self.locked.insert(index);
        }
        
        index
    }
}

/// This module contains only constants that are used by the Outline struct
mod consts {
    use vector::Vector3D;
    
    /// We want the field to be calculated with a certain number of divisions.
    /// If more accuracy is necessary, this can be increased.
    pub const FIELD_DIVISION: u32 = 20;
    
    /// We want the field to be have certain bounds. We also want the bounds in
    /// all three dimensions to be exactly the same. Thus, we will have a cube
    /// outside the body.
    pub const FIELD_BOUND: f64 = 1.0;
    
    /// These number of attraction/repulsion particles to create in the field
    pub const NUM_PARTICLES: usize = 150;
    
    /// We want a margin (a place near the border of the cube of the field
    /// where no particles are ever generated), which reduces the amount of
    /// space the particles can be located in.
    pub const FIELD_MARGIN: f64 = 0.05;
    
    /// Given the margin's size, the actual useful size of the fiels is reduced
    pub const FIELD_USEFUL_SIZE: f64 = 2.0 * (FIELD_BOUND - FIELD_MARGIN);
    
    /// We also want each particle to have at least some spread factor
    pub const MIN_SPREAD: f64 = 0.1;
    
    /// This decides how many divisions on the starting embryo should be
    /// calculated. 0 means use the simple icosahedron; higher values lead to
    /// further subdivision steps
    pub const EMBRYO_DIVISIONS: u32 = 3;
    
    /// The distance from the center of the field, which is where the vertices
    /// are initialized
    pub const VERTEX_START_DISTANCE: f64 = 0.05;
    
    /// The vertices respond to the field. This constant defines how strong
    /// this responsiveness is
    pub const FIELD_STRENGTH: f64 = 0.125;
    
    /// The vertices also respond to the presence of other vertices nearby.
    /// This constant defines how strong that responsiveness is
    pub const PROXIMITY_STRENGTH: f64 = 0.00026875;
    
    /// The maximum distance between two vertices for them to have an effect on
    /// each other
    pub const PROXIMITY_MAX_DISTANCE: f64 = 0.3;
    
    /// The proximity contribution to the force sensed by each vertex is not
    /// continuous for vertices that cross the threshold of
    /// PROXIMITY_MAX_DISTANCE; closer vertices have a positive effect, further
    /// vertices have a null effect. To smooth this, we make a small correction
    /// to the proximity contribution by removing just enough so that vertices
    /// exactly on the threshold also have a null effect
    /// FIXME: This is a function atm, since powi is not valid in a const
    #[inline(always)]
    pub fn PROXIMITY_CORRECTION() -> f64 {
        PROXIMITY_STRENGTH / PROXIMITY_MAX_DISTANCE.powi(3)
    }
    
    /// This constant defines the strength governing the movement of a vertex
    /// towards the point in between its neighbours
    pub const NEIGHBORS_STRENGTH: f64 = 30.0;
    
    /// The preferred distance between the vertices that share an edge
    pub const NEIGHBOR_DISTANCE: f64 = 0.11;
        
    /// The maximum distance between a vertex and the field's edge to activate
    /// the edge repulsion
    pub const EDGE_THICKNESS: f64 = 0.025;
    
    /// Vertices must stay away from the edges of the world. These constants
    /// define how strong that repulsive force is
    pub const EDGE_STRENGTH: f64 = 0.2;
    
    /// ???
    pub const SYMMETRY_EDGE_STRENGTH: f64 = 20.0;
    
    /// Vertices experience drag, which reduces their velocity
    pub const DRAG: f64 = 0.8;
    
    /// The maximum allowed velocity in each time frame
    pub const MAX_VELOCITY: f64 = 0.05;
    
    /// Time step
    pub const TIME_STEP: f64 = 0.1;
    
    /// Whether the vertices are to be centered at the end of each iteration
    pub const CENTER_SERVERS: bool = true;
    
    /// THe center of the field
    pub const FIELD_CENTER: Vector3D = Vector3D(0.0, 0.0, 0.0);
    
    /// The minimum distance that leads an edge to want to split in two
    pub const SPLIT_DISTANCE: f64 = 0.1; // Note: Not used???
    
    /// The minimum value for the norm of the vectorial difference between the
    /// velocities of the two vertices of an edge that must be observed before
    /// the edge splits
    pub const SPLIT_VELOCITY_DIFFERENCE: f64 = 0.00; // Note: Not used???
    
    /// The minimum length of an edge; an edge smaller than this is collapsed
    pub const MERGE_DISTANCE: f64 = 0.001; // Note: Not used???
    
    /// The number of divisions made on the grid of vertices. Notice that the
    /// number of cells in the grid is this number cubed
    pub const GRID_DIVISIONS: u32 = 5;
}

pub struct Outline {
    /// The field where the vertices exist
    field: Field,
    
    /// The position of each vertex
    pub positions: Vec<Vector3D>,
    
    /// The current velocity of each vertex
    velocities: Vec<Vector3D>,
    
    /// The neighbours of each vertex
    neighbours: Vec<Vec<usize>>,
    
    /// The set of vertices (identified by thie index in the positions and
    /// velocities vectors) that are locked in the symmetry axis x = 0
    pub locked: HashSet<usize>,
    
    /// The faces of the model, for later reconstruction
    pub faces: Vec<(usize, usize, usize)>,
}

impl Outline {
    pub fn new() -> Outline {
        let field = Field::new(consts::FIELD_DIVISION, consts::FIELD_BOUND);
        
        let mesh = HalfIcosphere::new(
            consts::EMBRYO_DIVISIONS, consts::VERTEX_START_DISTANCE);
        
        let positions = mesh.vertices;
        let velocities = vec![Vector3D(0.0, 0.0, 0.0); positions.len()];
        let locked = mesh.locked;
        let faces = mesh.faces.iter().cloned().collect();
        
        let mut neighbours = vec![vec![]; positions.len()];
        for &(i, j) in mesh.pair_to_edge.keys() {
            neighbours[i].push(j);
            neighbours[j].push(i);
        }
        
        Outline {
            field: field,
            positions: positions,
            velocities: velocities,
            neighbours: neighbours,
            locked: locked,
            faces: faces
        }
    }
    
    pub fn init(&mut self, mut genome: Genome) {
        for _ in 0..consts::NUM_PARTICLES {
            // Get coordinates from the genome
            let x = -consts::FIELD_BOUND + consts::FIELD_MARGIN +
                    genome.get_float() * consts::FIELD_USEFUL_SIZE;
            let y = -consts::FIELD_BOUND + consts::FIELD_MARGIN +
                    genome.get_float() * consts::FIELD_USEFUL_SIZE;
            let z = -consts::FIELD_BOUND + consts::FIELD_MARGIN +
                    genome.get_float() * consts::FIELD_USEFUL_SIZE;
            
            // attraction is a "corrected" gaussian value with mean 0 and
            // standard deviation 1. The correction is made to make sure the
            // value is between -1 and 1
            let attraction = convert_to_gaussian(genome.get_float(),
                                                 0.0, 1.0, -1.0, 1.0);
            
            // The spread is also a "corrected" gaussian with mean 0 and
            // standard deviation 1, with a minimum value of 1
            // We scale the value so that we have a value not restricted to be
            // greater than 1 but greater than the minimum spread
            let spread = convert_to_gaussian(genome.get_float(),
                                             0.0, 1.0, 1.0, f64::INFINITY)
                         * consts::MIN_SPREAD;
            
            let particle = Particle::new(Vector3D(x, y, z), attraction, spread);
            self.field.add_particle(&particle);
            
            println!("# DEBUG Particle");
            println!("#   x, y, z = {:.8}, {:.8}, {:.8}", x, y, z);
            println!("#   attraction = {:.8}", attraction);
            println!("#   spread = {:.8}", spread);
        }
        
        // We also add a central repulsive particle so that the body contains a
        // central mass. This seems to be good to create nice looking bodies
        let attraction = convert_to_gaussian(genome.get_float(),
                                             1.0, 1.0, -1.0, 0.0);
        let spread = convert_to_gaussian(genome.get_float(),
                                         0.0, 1.0, 1.0, f64::INFINITY)
                     * consts::MIN_SPREAD * 2.0;
        let particle = Particle::new(Vector3D(0.0, 0.0, 0.0),
                                     attraction, spread);
        self.field.add_particle(&particle);
    }
    
    pub fn step(&mut self) {
        for index in 0..self.positions.len() {
            self.update_velocities(index);
        }
        
        self.move_vertices();
        self.correct_astray_vertices();
    }
    
    
    fn update_velocities(&mut self, index: usize) {
        let pos = self.positions[index];
        
        // Where is the field sending me to?
        let field_acc = consts::FIELD_STRENGTH * self.field.get_nabla(pos);
        
        // This vertex wants to run away from other nearby vertices.
        // Since we are dealling with a symmetric body by ignoring half of the
        // outline, the points with x < 0 are not actually represented in
        // memory and so we must emulate their presence
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
            
            // We also need to react to the mirror image of this vertex
            // We only react to the mirror image if the vertex is not locked,
            // since being locked means that the vertex does not have a mirror
            // image.
            if !self.locked.contains(&other) {
                let Vector3D(new_x, new_y, new_z) = self.positions[other];
                let mirrored = Vector3D(-new_x, new_y, new_z);
                let delta = pos - mirrored;
                let dist = delta.norm();
                let factor = consts::PROXIMITY_STRENGTH / (dist * dist * dist);
                
                proximity_acc = proximity_acc + factor * delta;
            }
        }
        
        // Another force at play is the attraction to the neighbours.
        let mut pivot = Vector3D(0.0, 0.0, 0.0);
        for &neighbour in self.neighbours[index].iter() {
            pivot = pivot + self.positions[neighbour];
        }
        let n_neighbours = self.neighbours[index].len() as f64;
        let neighbour_acc = (pivot / n_neighbours - pos)
                           * consts::NEIGHBORS_STRENGTH / n_neighbours;
        
        // We also want to be stay clear of the edge of the field
        let mut edge_acc = Vector3D(0.0, 0.0, 0.0);
        if pos.0 > consts::FIELD_BOUND - consts::EDGE_THICKNESS {
            edge_acc.0 -= consts::EDGE_STRENGTH;
        }
        
        if pos.1 < -consts::FIELD_BOUND + consts::EDGE_THICKNESS {
            edge_acc.1 += consts::EDGE_STRENGTH;
        }
        else if pos.1 > consts::FIELD_BOUND - consts::EDGE_THICKNESS {
            edge_acc.1 -= consts::EDGE_STRENGTH;
        }
        
        if pos.2 < -consts::FIELD_BOUND + consts::EDGE_THICKNESS {
            edge_acc.2 += consts::EDGE_STRENGTH;
        }
        else if pos.2 > consts::FIELD_BOUND - consts::EDGE_THICKNESS {
            edge_acc.2 -= consts::EDGE_STRENGTH;
        }
        
        
        // Now let's merge together all the effects into a single acceleration
        // and update the velocity of the vertex
        let vel = self.velocities[index];
        let acc = field_acc + proximity_acc + neighbour_acc + edge_acc;
        
        if index == 0 {
            println!("# DEBUG index=0;");
            println!("#   field_acc = {:?}", field_acc);
            println!("#   proximity_acc = {:?}", proximity_acc);
            println!("#   neighbour_acc = {:?}", neighbour_acc);
            println!("#   edge_acc = {:?}", edge_acc);
        }
        
        // We take into account both the drag and the time step
        let mut new_vel = consts::DRAG * (vel + acc * consts::TIME_STEP);
        if self.locked.contains(&index) {
            new_vel.0 = 0.0;
        }
        
        let speed = new_vel.norm();
        if speed > consts::MAX_VELOCITY {
            new_vel = new_vel / speed * consts::MAX_VELOCITY;
        }
        
        self.velocities[index] = new_vel;
    }
    
    fn move_vertices(&mut self) {
        for index in 0..self.positions.len() {
            let pos = self.positions[index];
            let vel = self.velocities[index];
            self.positions[index] = pos + vel * consts::TIME_STEP;
        }
    }
    
    fn correct_astray_vertices(&mut self) {
        // Make sure that no vertex crosses the symmetry plane
        for index in 0..self.positions.len() {
            if self.positions[index].0 < 0.0 {
                self.positions[index].0 = 0.0;
            }
        }
    }
}


/// Converts a floating point number in the open interval (0, 1) to a normal
/// distribution percentile, with an interval.
fn convert_to_gaussian(f: f64, mean: f64, sigma: f64,
                       start: f64, end: f64) -> f64 {
    use utils::cdf::{cdf, inverse_cdf};
    
    let start =
        if start == f64::NEG_INFINITY {
            0.0
        }
        else {
            cdf((start - mean) / sigma)
        };
    let end =
        if end == f64::INFINITY {
            1.0
        }
        else {
            cdf((end - mean) / sigma)
        };
    
    // Now, we map the given float (which is between 0 and 1, excluding
    // ends) into the interval (start, end)
    let f = interpolate(start, end, f);
    
    // We then use the inverse_cdf of the standard normal distribution
    // to convert the `f` (which comes from a uniform distribution)
    // into another value (from the standard normal distribution)
    let normal = inverse_cdf(f);

    // Finally, convert the standard distribution to the one with
    // the given parameters (mean and standard deviation)
    normal * sigma + mean
}


