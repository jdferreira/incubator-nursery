use vector::Vector3D;
use std::f64;

// TODO: Attention!
// The field is hardcoded to be a repetition of the cube from (0, 0, 0) to (1, 1, 1)
// With a symmetry around the plane `x = 0` For this reason, all particles **must** be located
// within these bounds.

#[derive(Debug)]
pub struct Particle {
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
    fn get_particle_nabla(&self, point: Vector3D) -> Vector3D {
        let Vector3D(x, y, z) = point;
        
        // If the point is on the negative side of the symmetry plane, mirror it, find the nabla,
        // and return a mirror of the nabla
        if x < 0.0 {
            let result = self.get_particle_nabla(Vector3D(-x, y, z));
            let Vector3D(x, y, z) = result;
            return Vector3D(-x, y, z);
        }
        
        // We must first find the closest instance of this particle to the location `point` given
        // as argument.
        
        let (x_d, y_d, z_d) = (x - x.floor(), y - y.floor(), z - z.floor());
        
        // If we are exactly in between two instances, the nabla is 0
        if x_d == 0.5 || y_d == 0.5 || z_d == 0.5 { return Vector3D(0.0, 0.0, 0.0); }
        
        // Determine which cells around this position we need to test. This minimizes the number of
        // square roots that need to ba calculated from 27 to at most 8 (sometimes 4)
        let slice_0 = &[ 0.0];
        let slice_n = &[-1.0, 0.0];
        let slice_p = &[ 0.0, 1.0];
        
        let x_iters: &[f64] = {
            // This must be a generic slice without size, since it can be of different sizes.
            // The other iters can have size associated with them. The compiler could possibly
            // optimize some stuff with that infromation (one can wish :)
            if x < 1.0 { if x_d < 0.5 { slice_0 } else { slice_p } }
            else       { if x_d < 0.5 { slice_n } else { slice_p } }
        };
        let y_iters = if y_d < 0.5 { slice_n } else { slice_p };
        let z_iters = if z_d < 0.5 { slice_n } else { slice_p };
        
        let mut instance = Vector3D(0.0, 0.0, 0.0);
        let mut min_distance = f64::INFINITY;
        
        for &x0 in x_iters.iter() {
            for &y0 in y_iters.iter() {
                for &z0 in z_iters.iter() {
                    let this_instance = self.position + Vector3D(x0, y0, z0);
                    let dist = (point - this_instance).norm();
                    if min_distance == f64::INFINITY || dist < min_distance {
                        instance = this_instance;
                        min_distance = dist;
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
        
        let vector = self.position - point;
        let dist_square = vector.norm_square();
        let value = self.attraction *
                    (-dist_square / self.spread_square / 2f64).exp();
        
        factor * value / self.spread_square * vector
    }
}

// /// Contains a 3 dimensional field, where each point in the field can be
// /// associated with a field value and the three derivatives (with respect to
// /// each of the three axes). The field always represents a cube from
// /// (-bound, -bound, -bound) to (bound, bound, bound). Points outside this cube
// /// do not have a value nor derivatives
// pub struct Field {
//     /// The number of cells that are used to divde the cube field in a regular
//     /// 3D grid
//     ndivisions: usize,
    
//     /// Half the size of the cube that this field is contained in. Since the
//     /// cube is centered at (0, 0, 0), the coordinates of the cube are thus
//     /// (-bound, -bound, -bound) to (bound, bound, bound).
//     bound: f64,
    
//     /// The size of each of the cells in the 3D grid
//     dsize: f64,
    
//     /// The value of the field at the exact points that make up the grid
//     /// plus the three derivatives
//     derivatives: Vec<Vector3D>
// }

// impl Field {
//     pub fn new(ndivisions: u32, bound: f64) -> Field {
//         assert!(ndivisions > 0);
//         assert!(bound > 0.0);
        
//         let ndivisions = ndivisions as usize;
//         let npoints = (ndivisions + 1).pow(3);
//         let dsize = 2.0 * bound / (ndivisions as f64);
//         let derivatives = vec![Vector3D(0.0, 0.0, 0.0); npoints];
        
//         Field {
//             ndivisions: ndivisions,
//             bound: bound,
//             dsize: dsize,
//             derivatives: derivatives
//         }
//     }
    
//     /// This function receives a particle and adds its effect to the field.
//     pub fn add_particle(&mut self, particle: &Particle) {
//         for x0 in 0usize..(self.ndivisions + 1) {
//             let x = -self.bound + (x0 as f64) * self.dsize;
//             for y0 in 0usize..(self.ndivisions + 1) {
//                 let y = -self.bound + (y0 as f64) * self.dsize;
//                 for z0 in 0usize..(self.ndivisions + 1) {
//                     let z = -self.bound + (z0 as f64) * self.dsize;
//                     let index = x0 * self.ndivisions * self.ndivisions +
//                                 y0 * self.ndivisions + z0;
//                     self.derivatives[index] =
//                         self.derivatives[index] +
//                         particle.get_particle_nabla(Vector3D(x, y, z));
//                 }
//             }
//         }
//     }
        
//     /// Get the value of the field at the given position, interpolating as
//     /// appropriate. This returns a tuple (value, dx, dy, dz).
//     pub fn get_nabla(&self, Vector3D(x, y, z): Vector3D) -> Vector3D {
        
//         if x < -self.bound || y < -self.bound || z < -self.bound ||
//            x >  self.bound || y >  self.bound || z >  self.bound {
//             return Vector3D(0.0, 0.0, 0.0);
//         }
        
//         // These are the values normalized between 0 and ndivisions
//         let x = (x + self.bound) / self.dsize;
//         let y = (y + self.bound) / self.dsize;
//         let z = (z + self.bound) / self.dsize;
        
//         // These are the integer parts of that. If the integer part is equal to
//         // the actual normalized value, we are over an axis and we don't need
//         // to interpolate; otherwise, we do need to interpolate.
//         let x0 = x.floor() as usize;
//         let y0 = y.floor() as usize;
//         let z0 = z.floor() as usize;
        
//         // Are we on a grid edge?
//         let (on_x, on_y, on_z) = (x == x0 as f64, y == y0 as f64, z == z0 as f64);
        
//         // The next grid index, for interpolation
//         let (x1, y1, z1) = (x0 + 1, y0 + 1, z0 + 1);
        
//         // The amounts to interpolate
//         let (x_d, y_d, z_d) = (x - x0 as f64, y - y0 as f64, z - z0 as f64);
        
//         if on_x && on_y && on_z {
//             self.derivate0(x0, y0, z0)
//         }
//         else if on_x && on_y {
//             let start = self.derivate0(x0, y0, z0);
//             let end   = self.derivate0(x0, y0, z1);
//             interpolate(start, end, z_d)
//         }
//         else if on_x && on_z {
//             let start = self.derivate0(x0, y0, z0);
//             let end   = self.derivate0(x0, y1, z0);
//             interpolate(start, end, y_d)
//         }
//         else if on_y && on_z {
//             let start = self.derivate0(x0, y0, z0);
//             let end   = self.derivate0(x1, y0, z0);
//             interpolate(start, end, x_d)
//         }
//         else if on_x {
//             let start = {
//                 let start = self.derivate0(x0, y0, z0);
//                 let end   = self.derivate0(x0, y1, z0);
//                 interpolate(start, end, y_d)
//             };
//             let end = {
//                 let start = self.derivate0(x0, y0, z1);
//                 let end   = self.derivate0(x0, y1, z1);
//                 interpolate(start, end, y_d)
//             };
//             interpolate(start, end, z_d)
//         }
//         else if on_y {
//             let start = {
//                 let start = self.derivate0(x0, y0, z0);
//                 let end   = self.derivate0(x0, y0, z1);
//                 interpolate(start, end, z_d)
//             };
//             let end = {
//                 let start = self.derivate0(x1, y0, z0);
//                 let end   = self.derivate0(x1, y0, z1);
//                 interpolate(start, end, z_d)
//             };
//             interpolate(start, end, x_d)
//         }
//         else if on_z {
//             let start = {
//                 let start = self.derivate0(x0, y0, z0);
//                 let end   = self.derivate0(x1, y0, z0);
//                 interpolate(start, end, x_d)
//             };
//             let end = {
//                 let start = self.derivate0(x0, y1, z0);
//                 let end   = self.derivate0(x1, y1, z0);
//                 interpolate(start, end, x_d)
//             };
//             interpolate(start, end, y_d)
//         }
//         else {
//             let start = {
//                 let start = {
//                     let start = self.derivate0(x0, y0, z0);
//                     let end   = self.derivate0(x1, y0, z0);
//                     interpolate(start, end, x_d)
//                 };
//                 let end = {
//                     let start = self.derivate0(x0, y1, z0);
//                     let end   = self.derivate0(x1, y1, z0);
//                     interpolate(start, end, x_d)
//                 };
//                 interpolate(start, end, y_d)
//             };
            
//             let end = {
//                 let start = {
//                     let start = self.derivate0(x0, y0, z1);
//                     let end   = self.derivate0(x1, y0, z1);
//                     interpolate(start, end, x_d)
//                 };
//                 let end = {
//                     let start = self.derivate0(x0, y1, z1);
//                     let end   = self.derivate0(x1, y1, z1);
//                     interpolate(start, end, x_d)
//                 };
//                 interpolate(start, end, y_d)
//             };
            
//             interpolate(start, end, z_d)
//         }
//     }
    
//     #[inline(always)]
//     fn derivate0(&self, x0: usize, y0: usize, z0: usize) -> Vector3D {
//         self.derivatives[x0 * self.ndivisions * self.ndivisions +
//                     y0 * self.ndivisions + z0]
//     }
// }


#[test]
fn test() {
    let p1 = Particle::new(Vector3D(0.3, 0.8, 0.0), 1.0, 1.0);
    let p2 = Particle::new(Vector3D(0.1, 0.5, 0.0), 1.0, 1.0);
    let p3 = Particle::new(Vector3D(0.5, 0.1, 0.0), 1.0, 1.0);
    let p4 = Particle::new(Vector3D(0.8, 0.4, 0.0), 1.0, 1.0);
    
    let nabla = p1.get_particle_nabla(Vector3D(0.7, 0.3, 0.0)) +
                p2.get_particle_nabla(Vector3D(0.7, 0.3, 0.0)) +
                p3.get_particle_nabla(Vector3D(0.7, 0.3, 0.0)) +
                p4.get_particle_nabla(Vector3D(0.7, 0.3, 0.0));
    println!("nabla = {:?}", nabla);
}

