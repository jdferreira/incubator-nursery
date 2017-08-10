use vector::Vector3D;
use utils::math::interpolate;

#[derive(Debug)]
pub struct Particle {
    position: Vector3D,
    attraction: f64,
    spread_square: f64
}

impl Particle {
    pub fn new(position: Vector3D, attraction: f64, spread: f64) -> Particle {
        Particle {
            position: position,
            attraction: attraction,
            spread_square: spread.powi(2)
        }
    }
    
    /// For a given particule, compute its effect (dx, dy and dz) at
    /// the given 3D point
    fn get_particle_nabla(&self, p: Vector3D) -> Vector3D {
        let v = p - self.position;
        let dist_square = v.norm_square();
        let value = self.attraction *
                    (-dist_square / self.spread_square / 2f64).exp();
        
        -value / self.spread_square * v
        
        // let dx = -value * v.0 / self.spread_square;
        // let dy = -value * v.1 / self.spread_square;
        // let dz = -value * v.2 / self.spread_square;
        // Vector3D(dx, dy, dz)
    }
}

/// Contains a 3 dimensional field, where each point in the field can be
/// associated with a field value and the three derivatives (with respect to
/// each of the three axes). The field always represents a cube from
/// (-bound, -bound, -bound) to (bound, bound, bound). Points outside this cube
/// do not have a value nor derivatives
pub struct Field {
    /// The number of cells that are used to divde the cube field in a regular
    /// 3D grid
    ndivisions: usize,
    
    /// Half the size of the cube that this field is contained in. Since the
    /// cube is centered at (0, 0, 0), the coordinates of the cube are thus
    /// (-bound, -bound, -bound) to (bound, bound, bound).
    bound: f64,
    
    /// The size of each of the cells in the 3D grid
    dsize: f64,
    
    /// The value of the field at the exact points that make up the grid
    /// plus the three derivatives
    derivatives: Vec<Vector3D>
}

impl Field {
    pub fn new(ndivisions: u32, bound: f64) -> Field {
        assert!(ndivisions > 0);
        assert!(bound > 0.0);
        
        let ndivisions = ndivisions as usize;
        let npoints = (ndivisions + 1).pow(3);
        let dsize = 2.0 * bound / (ndivisions as f64);
        let derivatives = vec![Vector3D(0.0, 0.0, 0.0); npoints];
        
        Field {
            ndivisions: ndivisions,
            bound: bound,
            dsize: dsize,
            derivatives: derivatives
        }
    }
    
    /// This function receives a particle and adds its effect to the field.
    pub fn add_particle(&mut self, particle: &Particle) {
        for x0 in 0usize..(self.ndivisions + 1) {
            let x = -self.bound + (x0 as f64) * self.dsize;
            for y0 in 0usize..(self.ndivisions + 1) {
                let y = -self.bound + (y0 as f64) * self.dsize;
                for z0 in 0usize..(self.ndivisions + 1) {
                    let z = -self.bound + (z0 as f64) * self.dsize;
                    let index = x0 * self.ndivisions * self.ndivisions +
                                y0 * self.ndivisions + z0;
                    self.derivatives[index] =
                        self.derivatives[index] +
                        particle.get_particle_nabla(Vector3D(x, y, z));
                }
            }
        }
    }
        
    /// Get the value of the field at the given position, interpolating as
    /// appropriate. This returns a tuple (value, dx, dy, dz).
    pub fn get_nabla(&self, Vector3D(x, y, z): Vector3D) -> Vector3D {
        
        if x < -self.bound || y < -self.bound || z < -self.bound ||
           x >  self.bound || y >  self.bound || z >  self.bound {
            return Vector3D(0.0, 0.0, 0.0);
        }
        
        // These are the values normalized between 0 and ndivisions
        let x = (x + self.bound) / self.dsize;
        let y = (y + self.bound) / self.dsize;
        let z = (z + self.bound) / self.dsize;
        
        // These are the integer parts of that. If the integer part is equal to
        // the actual normalized value, we are over an axis and we don't need
        // to interpolate; otherwise, we do need to interpolate.
        let x0 = x.floor() as usize;
        let y0 = y.floor() as usize;
        let z0 = z.floor() as usize;
        
        // Are we on a grid edge?
        let (on_x, on_y, on_z) = (x == x0 as f64, y == y0 as f64, z == z0 as f64);
        
        // The next grid index, for interpolation
        let (x1, y1, z1) = (x0 + 1, y0 + 1, z0 + 1);
        
        // The amounts to interpolate
        let (x_d, y_d, z_d) = (x - x0 as f64, y - y0 as f64, z - z0 as f64);
        
        if on_x && on_y && on_z {
            self.derivate0(x0, y0, z0)
        }
        else if on_x && on_y {
            let start = self.derivate0(x0, y0, z0);
            let end   = self.derivate0(x0, y0, z1);
            interpolate(start, end, z_d)
        }
        else if on_x && on_z {
            let start = self.derivate0(x0, y0, z0);
            let end   = self.derivate0(x0, y1, z0);
            interpolate(start, end, y_d)
        }
        else if on_y && on_z {
            let start = self.derivate0(x0, y0, z0);
            let end   = self.derivate0(x1, y0, z0);
            interpolate(start, end, x_d)
        }
        else if on_x {
            let start = {
                let start = self.derivate0(x0, y0, z0);
                let end   = self.derivate0(x0, y1, z0);
                interpolate(start, end, y_d)
            };
            let end = {
                let start = self.derivate0(x0, y0, z1);
                let end   = self.derivate0(x0, y1, z1);
                interpolate(start, end, y_d)
            };
            interpolate(start, end, z_d)
        }
        else if on_y {
            let start = {
                let start = self.derivate0(x0, y0, z0);
                let end   = self.derivate0(x0, y0, z1);
                interpolate(start, end, z_d)
            };
            let end = {
                let start = self.derivate0(x1, y0, z0);
                let end   = self.derivate0(x1, y0, z1);
                interpolate(start, end, z_d)
            };
            interpolate(start, end, x_d)
        }
        else if on_z {
            let start = {
                let start = self.derivate0(x0, y0, z0);
                let end   = self.derivate0(x1, y0, z0);
                interpolate(start, end, x_d)
            };
            let end = {
                let start = self.derivate0(x0, y1, z0);
                let end   = self.derivate0(x1, y1, z0);
                interpolate(start, end, x_d)
            };
            interpolate(start, end, y_d)
        }
        else {
            let start = {
                let start = {
                    let start = self.derivate0(x0, y0, z0);
                    let end   = self.derivate0(x1, y0, z0);
                    interpolate(start, end, x_d)
                };
                let end = {
                    let start = self.derivate0(x0, y1, z0);
                    let end   = self.derivate0(x1, y1, z0);
                    interpolate(start, end, x_d)
                };
                interpolate(start, end, y_d)
            };
            
            let end = {
                let start = {
                    let start = self.derivate0(x0, y0, z1);
                    let end   = self.derivate0(x1, y0, z1);
                    interpolate(start, end, x_d)
                };
                let end = {
                    let start = self.derivate0(x0, y1, z1);
                    let end   = self.derivate0(x1, y1, z1);
                    interpolate(start, end, x_d)
                };
                interpolate(start, end, y_d)
            };
            
            interpolate(start, end, z_d)
        }
    }
    
    #[inline(always)]
    fn derivate0(&self, x0: usize, y0: usize, z0: usize) -> Vector3D {
        self.derivatives[x0 * self.ndivisions * self.ndivisions +
                    y0 * self.ndivisions + z0]
    }
}
