use std::f64;
use vector::Vector3D;
use utils::math::interpolate;

// TODO: Attention!
// The field is hardcoded to be a repetition of the cube from (0, 0, 0) to (1, 1, 1)
// With a symmetry around the plane `x = 0` For this reason, all particles **must** be located
// within these bounds.

#[derive(Debug)]
struct Particle {
    position: Vector3D,
    attraction: f64,
    spread_square: f64
}

impl Particle {
    pub fn new(position: Vector3D, attraction: f64, spread: f64) -> Particle {
        let Vector3D(x, y, z) = position;
        assert!(0.0 <= x && x < 1.0);
        assert!(0.0 <= y && y < 1.0);
        assert!(0.0 <= z && z < 1.0);
        
        Particle {
            position: position,
            attraction: attraction,
            spread_square: spread.powi(2)
        }
    }
    
    /// For a given particule, compute its effect (dx, dy and dz) at the given position
    /// NOTE: private method. Only gets called with positive x value
    fn get_particle_nabla(&self, point: Vector3D) -> Vector3D {
        // We must first find the closest instance of this particle to the location `point` given
        // as argument.
        
        let Vector3D(x, y, z) = point;
        let (xd, yd, zd) = (x - x.floor(), y - y.floor(), z - z.floor());
        
        // If we are exactly in between two instances, the nabla is (0, 0, 0)
        if xd == 0.5 || yd == 0.5 || zd == 0.5 { return Vector3D(0.0, 0.0, 0.0); }
        
        // Determine which cells around this position we need to test. This minimizes the number of
        // square roots that need to ba calculated from 27 to at most 8 (sometimes 4)
        #[allow(non_upper_case_globals)]
        static slice_0: [f64; 1] = [ 0.0];
        #[allow(non_upper_case_globals)]
        static slice_n: [f64; 2] = [-1.0, 0.0];
        #[allow(non_upper_case_globals)]
        static slice_p: [f64; 2] = [ 0.0, 1.0];
        
        let x_iters: &[_] = {
            if x < 1.0 { if xd < 0.5 { &slice_0 } else { &slice_p } }
            else       { if xd < 0.5 { &slice_n } else { &slice_p } }
        };
        let y_iters = if yd < 0.5 { &slice_n } else { &slice_p };
        let z_iters = if zd < 0.5 { &slice_n } else { &slice_p };
        
        let mut instance = Vector3D(0.0, 0.0, 0.0);
        let mut dist_square = f64::INFINITY;
        
        for &x0 in x_iters.iter() {
            for &y0 in y_iters.iter() {
                for &z0 in z_iters.iter() {
                    let this_instance = self.position + Vector3D(x0, y0, z0);
                    let this_dist_square = (point - this_instance).norm_square();
                    if dist_square == f64::INFINITY || this_dist_square < dist_square {
                        instance = this_instance;
                        dist_square = this_dist_square;
                    }
                }
            }
        };
        
        let factor = {
            let Vector3D(xi, yi, zi) = instance;
            let mut factor = (x - xi).abs();
            factor = f64::max(factor, (y - yi).abs());
            factor = f64::max(factor, (z - zi).abs());
            0.5 - factor
        };
        
        let vector = instance - point;
        let value = self.attraction *
                    (-dist_square / self.spread_square / 2f64).exp() / self.spread_square;
        
        factor * value * vector
    }
}

/// Contains a 3 dimensional infinite and symmetric field (around the plane `x = 0`, where each
/// point in the field is associated with a field value and the three derivatives (with respect to
/// each of the three axes) of that value. The field is an infinitely repeating pattern of the cube
/// that goes from (0, 0, 0) to (1, 1, 1). Points outside this cube on the positive side of the
/// symmetry plane have the same value and derivative as the correspoding value on the cube (using
/// modulo airithmetic); points on the negative side of the symmetry plane have the same value and
/// derivatives as their mirror image (with a symmetric x component to account for symmetry).
///
/// The values in the field are dictated by a set of particles, which have a position, an attraction
/// and a spread value, and follow a gaussian decay from that position to infinity.
///
/// To account for the infinitude (and repeatability) of the field, a point only senses the effect
/// of the closest instance of each particle. For example, a point that is exaclty in between two
/// instances of the same particle senses no effect from it.
///
/// The field is divided in a 3D grid, and the derivatives at each grid intersection are calculated
/// and stored. Requests for the derivatives on any point are then interpolations of these values.
/// This makes the requests much faster (at most 7 interpolations) since they do not have to perform
/// any significant mathematical operations between the points and the particles in the field.
pub struct Field {
    /// The number of cells that are used to divde the cube field in a regular 3D grid
    ndivisions: usize,
    
    /// The size of each of the cells in the 3D grid
    dsize: f64,
    
    /// The derivatives at the exact points that make up the grid, on the cube that has a face
    /// on the symmetry plane
    ds_symmetry: Vec<Vector3D>,
    
    /// The derivatives at the exact points that make up the grid, on the cube that does not have
    /// a face on the symmetry plane
    ds_elsewhere: Vec<Vector3D>,
}

impl Field {
    pub fn new(ndivisions: u32) -> Field {
        let ndivisions = ndivisions as usize;
        let npoints = (ndivisions + 1).pow(3);
        let dsize = 1.0 / (ndivisions as f64);
        let ds_symmetry = vec![Vector3D(0.0, 0.0, 0.0); npoints];
        let ds_elsewhere = vec![Vector3D(0.0, 0.0, 0.0); npoints];
        
        Field {
            ndivisions: ndivisions,
            dsize: dsize,
            ds_symmetry: ds_symmetry,
            ds_elsewhere: ds_elsewhere,
        }
    }
    
    /// This function receives a particle and adds its effect to the field.
    pub fn add_particle(&mut self, position: Vector3D, attraction: f64, spread: f64) {
        let particle = Particle::new(position, attraction, spread);
        
        for x0 in 0usize..self.ndivisions {
            let x = x0 as f64 * self.dsize;
            for y0 in 0usize..self.ndivisions {
                let y = y0 as f64 * self.dsize;
                for z0 in 0usize..self.ndivisions {
                    let z = z0 as f64 * self.dsize;
                    let index = x0 * self.ndivisions * self.ndivisions + y0 * self.ndivisions + z0;
                    
                    // Find the nabla on the cube adjacent to the symmetry plane and on the
                    // other cube
                    self.ds_symmetry[index] =
                        self.ds_symmetry[index] +
                        particle.get_particle_nabla(Vector3D(x, y, z));
                    self.ds_elsewhere[index] =
                        self.ds_elsewhere[index] +
                        particle.get_particle_nabla(Vector3D(x + 1.0, y, z));
                }
            }
        }
    }
        
    /// Get the value of the field at the given position, interpolating as necessary.
    pub fn get_nabla(&self, Vector3D(x, y, z): Vector3D) -> Vector3D {
        // If the point is on the negative side of the symmetry plane, mirror it, find the
        // nabla, and return a mirror of the nabla
        if x < 0.0 {
            let result = self.get_nabla(Vector3D(-x, y, z));
            let Vector3D(x, y, z) = result;
            return Vector3D(-x, y, z);
        }
        
        // Determine whether we want to take the nabla from the cube near the symmetry plane or not
        let ds_store =
            if x < 1.0 { &self.ds_symmetry }
            else { &self.ds_elsewhere };
        
        // Modulo arithmetic: remove the integer part
        let x = x - x.floor();
        let y = y - y.floor();
        let z = z - z.floor();
        
        // These are the values normalized between 0 and ndivisions
        let x = x / self.dsize;
        let y = y / self.dsize;
        let z = z / self.dsize;
        
        // These are the integer parts of that. If the integer part is equal to
        // the actual normalized value, we are over an axis and we don't need
        // to interpolate; otherwise, we do need to interpolate.
        let (x0, y0, z0) = (x.floor() as usize, y.floor() as usize, z.floor() as usize);
        
        // The next grid index, for interpolation. We need to take the remainder since we may end
        // up spilliong into the next cube.
        let x1 = (x0 + 1) % self.ndivisions;
        let y1 = (y0 + 1) % self.ndivisions;
        let z1 = (z0 + 1) % self.ndivisions;
        
        // The amounts to interpolate
        let (xd, yd, zd) = (x - x0 as f64, y - y0 as f64, z - z0 as f64);
        
        // The indices in the grid
        let i000 = x0 * self.ndivisions * self.ndivisions + y0 * self.ndivisions + z0;
        let i001 = x0 * self.ndivisions * self.ndivisions + y0 * self.ndivisions + z1;
        let i010 = x0 * self.ndivisions * self.ndivisions + y1 * self.ndivisions + z0;
        let i011 = x0 * self.ndivisions * self.ndivisions + y1 * self.ndivisions + z1;
        let i100 = x1 * self.ndivisions * self.ndivisions + y0 * self.ndivisions + z0;
        let i101 = x1 * self.ndivisions * self.ndivisions + y0 * self.ndivisions + z1;
        let i110 = x1 * self.ndivisions * self.ndivisions + y1 * self.ndivisions + z0;
        let i111 = x1 * self.ndivisions * self.ndivisions + y1 * self.ndivisions + z1;
        
        // The values. When x1 == 0, we need to take the values from ds_elsewhere unconditionally,
        // since the values on the ds_symmetry for x1 == 0 are different.
        let v000 = ds_store[i000];
        let v001 = ds_store[i001];
        let v010 = ds_store[i010];
        let v011 = ds_store[i011];
        let v100 = if x1 != 0 { ds_store } else { &self.ds_elsewhere }[i100];
        let v101 = if x1 != 0 { ds_store } else { &self.ds_elsewhere }[i101];
        let v110 = if x1 != 0 { ds_store } else { &self.ds_elsewhere }[i110];
        let v111 = if x1 != 0 { ds_store } else { &self.ds_elsewhere }[i111];
        
        let start = {
            let start = interpolate(v000, v100, xd);
            let end = interpolate(v010, v110, xd);
            interpolate(start, end, yd)
        };
        
        let end = {
            let start = interpolate(v001, v101, xd);
            let end = interpolate(v011, v111, xd);
            interpolate(start, end, yd)
        };
        
        interpolate(start, end, zd)
    }
    
}
