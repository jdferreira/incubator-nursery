use std::f64;
use proj;

static PI: f64 = f64::consts::PI;

#[derive(Copy,Clone,Debug)]
pub struct Vector3d(pub f64, pub f64, pub f64);

pub struct Particle {
    pub long: f64,
    pub lat: f64,
    pub f: f64,
    pub s: f64
}

impl Particle {
    
    pub fn new(long: f64, lat: f64, f: f64, s: f64) -> Particle {
        assert!(-1.0 <= f && f <= 1.0);
        assert!(s > 0.0);
        
        Particle { long: long, lat: lat, f: f, s: s }
    }
    
}

pub struct Surface {
    particles: Vec<Particle>,
    points: Vec<(f64, f64)>
}

impl Surface {
    
    pub fn new(particles: Vec<Particle>, points: Vec<(f64, f64)>) -> Surface {
        Surface { particles: particles, points: points }
    }
    
    pub fn img_buffer(&self, width: u32, height: u32) -> Vec<f64> {
        let mut result = vec![0.0; (width * height) as usize];
        
        let y_top = height / 16;
        let y_bottom = 15 * y_top;
        let chunk_width = width as f64 / 32.0;
        
        for y in 0..height {
            for x in 0..width {
                
                let index = (y * width + x) as usize;
                let mut long = x as f64 / width as f64 * 2.0 * PI - PI;
                let lat = PI / 2.0 - y as f64 / height as f64 * PI;

                // If this is a top or bottom pixel, we need to change the picture a bit, so that
                // the final image has triangles on the poles.
                if y < y_top || y >= y_bottom {
                    let how_far = if y < y_top {
                        (y_top - y) as f64 / y_top as f64
                    }
                    else {
                        (y - y_bottom) as f64 / y_top as f64
                    };
                    
                    let x_chunk = 32 * x / width;
                    let x_chunk_start = x_chunk as f64 * chunk_width;
                    let x_triangle_start = (x_chunk_start as f64 + how_far * chunk_width / 2.0) as u32;
                    let x_triangle_end   = (x_chunk_start as f64 + chunk_width * (1.0 - how_far / 2.0)) as u32;
                    
                    if x < x_triangle_start || x > x_triangle_end {
                        result[index] = f64::NEG_INFINITY;
                        continue;
                    }
                    
                    let f = ((x as f64 - x_triangle_start as f64) /
                             (x_triangle_end as f64 - x_triangle_start as f64));
                    let new_x = x_chunk_start as f64 + f * width as f64 / 32.0;
                    long = new_x / width as f64 * 2.0 * PI - PI;
                }
                                
                // And then the corresponding x, y, z coordinates
                let Vector3d(xc, yc, zc) = proj::to_cartesian(long, lat);
                
                // We need to get the combined effect of all particles sensed at this location
                let mut combined = 0.0;
                for p in self.particles.iter() {
                    let Vector3d(xp, yp, zp) = proj::to_cartesian(p.long, p.lat);
                    let cos_angle = xp * xc + yp * yc + zp * zc;
                    let dist = cos_angle.acos();
                    
                    combined += p.f * (-dist * dist / (2.0 * p.s * p.s)).exp();
                }
                
                result[index] = combined;
            }
        }
        
        result
    }
}
