pub mod consts;
mod delaunay;

use std::f64;
use vector::Vector3D;
use utils::math::convert_to_gaussian;
use genome::Genome;

static PI: f64 = f64::consts::PI;


/// Represents a point on the surface of the unit sphere
#[derive(Copy,Clone)]
pub struct SurfacePoint {
    long: f64,
    lat: f64
}

impl SurfacePoint {
    pub fn new(longitude: f64, latitude: f64) -> SurfacePoint {
        assert!(-PI / 2.0 <= latitude && latitude <= PI / 2.0);
        SurfacePoint { long: longitude, lat: latitude }
    }
    
    pub fn new_degree(longitude: f64, latitude: f64) -> SurfacePoint {
        let longitude = longitude * PI / 180.0;
        let latitude = latitude * PI / 180.0;
        SurfacePoint::new(longitude, latitude)
    }
}

impl From<SurfacePoint> for Vector3D {
    
    #[inline(always)]
    fn from(sp: SurfacePoint) -> Vector3D {
        let SurfacePoint { long, lat } = sp;
        
        let x = lat.cos() * long.cos();
        let y = lat.cos() * long.sin();
        let z = lat.sin();
        
        Vector3D(x, y, z)
    }
}

impl From<Vector3D> for SurfacePoint {
    
    #[inline(always)]
    fn from(v: Vector3D) -> SurfacePoint {
        let Vector3D(x, y, z) = v.normalize();
        
        let longitude = y.atan2(x);
        let latitude = z.asin();
        
        SurfacePoint { long: longitude, lat: latitude }
    }
}

/// Represents a particle on the surface of a sphere, with a certain position, force and spread
struct Particle {
    position: SurfacePoint,
    attraction: f64,
    spread: f64
}

impl Particle {
    fn new(position: SurfacePoint, attraction: f64, spread: f64) -> Particle {
        assert!(spread > 0.0);
        Particle { position: position, attraction: attraction, spread: spread }
    }
}


pub struct Surface {
    particles: Vec<Particle>,
    pub points: Vec<SurfacePoint>,
    velocities: Vec<Vector3D>,
    pub vertices_to_edges: Vec<Vec<usize>>
}

impl Surface {
    
    pub fn new() -> Surface {
        Surface {
            particles: vec![],
            points: vec![],
            velocities: vec![],
            vertices_to_edges: vec![],
        }
    }
    
    pub fn attach_genome(&mut self, genome: &mut Genome) {
        // FIXME: Add particles and points symmetrically around the plane x = 0
        
        // Add particles that affect the points on the surface
        for _ in 0..consts::NUM_PARTICLES {
            let long = genome.get_float() * 2.0 * PI;
            let lat = genome.get_float() * PI - PI / 2.0;
            let attraction = convert_to_gaussian(genome.get_float(), 0.0, 1.0, -1.0, 1.0);
            let spread = convert_to_gaussian(genome.get_float(), 0.0, 1.0, 1.0, f64::INFINITY) *
                         consts::MIN_SPREAD;
            self.add_particle(SurfacePoint::new(long, lat), attraction, spread);
        }
        
        // Add the points that will move around the surface
        for _ in 0..consts::NUM_POINTS {
            let long = genome.get_float() * 2.0 * PI;
            let lat = genome.get_float() * PI - PI / 2.0;
            self.add_point(SurfacePoint::new(long, lat));
        }
        
        // We know need to determine the neighbuorhood of each point.
        // We do this by executing a Delaunay triangulation on the points.
        let neighbours = delaunay::triangulate(&self.points);
        for (i, set) in neighbours.iter().enumerate() {
            for &j in set.iter() {
                self.vertices_to_edges[i].push(j);
            }
        }
    }
    
    fn add_particle(&mut self, position: SurfacePoint, attraction: f64, spread: f64) {
        self.particles.push(Particle::new(position, attraction, spread));
    }
    
    fn add_point(&mut self, point: SurfacePoint) -> usize {
        self.points.push(point);
        self.velocities.push(Vector3D(0.0, 0.0, 0.0));
        self.vertices_to_edges.push(vec![]);
        self.points.len() - 1
    }
    
    fn add_edge(&mut self, p1: usize, p2: usize) {
        let (p1, p2) =
            if p1 < p2 { (p1, p2) }
            else       { (p2, p1) };
        
        self.vertices_to_edges[p1].push(p2);
        self.vertices_to_edges[p2].push(p1);
    }
    
    pub fn step(&mut self) {
        for index in 0..self.points.len() {
            self.update_velocities(index);
        }
        
        self.move_points();
    }
    
    fn update_velocities(&mut self, index: usize) {
        let here = Vector3D::from(self.points[index]);
        
        let mut acc = Vector3D(0.0, 0.0, 0.0);
        
        // Effect of the particles
        for particle in self.particles.iter() {
            let there = Vector3D::from(particle.position);
            let k1 = here * there;
            let vector = there - k1 * here;
            let dist = k1.acos();
            let value = -particle.attraction *
                        (-dist * dist / (particle.spread * particle.spread * 2.0)).exp();
            let factor = consts::PARTICLE_STRENGTH * value / (particle.spread * particle.spread);
            acc = acc + factor * vector;
        }
        
        // Effect of neighbours
        // for &neighbour in self.vertices_to_edges[index].iter() {
        for other in 0..self.points.len() {
            if other == index { continue; }
            
            let there = Vector3D::from(self.points[other]);
            let k1 = here * there;
            let vector = there - k1 * here;
            let dist = k1.acos();
            let value = dist - PI;
            let factor = value * consts::NEIGHBORS_STRENGTH / (dist * dist * dist);
            acc = acc + factor * vector;
        }
        
        self.velocities[index] = consts::DRAG * (self.velocities[index] + acc * consts::TIME_STEP);
    }
    
    fn move_points(&mut self) {
        for (pos, &vel) in Iterator::zip(self.points.iter_mut(), self.velocities.iter()) {
            *pos = SurfacePoint::from(Vector3D::from(*pos) + vel * consts::TIME_STEP);
        }
    }
}
