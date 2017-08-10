use std::f64;
use std::ops::{Add,Sub,Mul,Div};

use image;
use std::fs::File;
use std::path::Path;

static PI: f64 = f64::consts::PI;


/// Represents a point in 3D space
#[derive(Clone,Copy,Debug)]
pub struct Vector3d(f64, f64, f64);

impl Vector3d {
    #[inline(always)]
    pub fn norm(self) -> f64 {
        self.norm_square().sqrt()
    }
    
    #[inline(always)]
    pub fn norm_square(self) -> f64 {
        self.0 * self.0 + self.1 * self.1 + self.2 * self.2
    }
    
    #[inline(always)]
    pub fn normalize(self) -> Vector3d {
        self / self.norm()
    }
}

impl Add for Vector3d {
    type Output = Vector3d;
    
    #[inline(always)]
    fn add(self, other: Vector3d) -> Vector3d {
        Vector3d(self.0 + other.0, self.1 + other.1, self.2 + other.2)
    }
}

impl Sub for Vector3d {
    type Output = Vector3d;
    
    #[inline(always)]
    fn sub(self, other: Vector3d) -> Vector3d {
        Vector3d(self.0 - other.0, self.1 - other.1, self.2 - other.2)
    }
}

impl Mul<f64> for Vector3d {
    type Output = Vector3d;
    
    #[inline(always)]
    fn mul(self, other: f64) -> Vector3d {
        Vector3d(self.0 * other, self.1 * other, self.2 * other)
    }
}

impl Mul<Vector3d> for f64 {
    type Output = Vector3d;
    
    #[inline(always)]
    fn mul(self, other: Vector3d) -> Vector3d {
        Vector3d(self * other.0, self * other.1, self * other.2)
    }
}

impl Mul<Vector3d> for Vector3d {
    type Output = f64;
    
    #[inline(always)]
    fn mul(self, other: Vector3d) -> f64 {
        self.0 * other.0 + self.1 * other.1 + self.2 * other.2
    }
}

impl Div<f64> for Vector3d {
    type Output = Vector3d;
    
    #[inline(always)]
    fn div(self, other: f64) -> Vector3d {
        Vector3d(self.0 / other, self.1 / other, self.2 / other)
    }
}


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

impl From<SurfacePoint> for Vector3d {
    
    #[inline(always)]
    fn from(sp: SurfacePoint) -> Vector3d {
        let SurfacePoint { long, lat } = sp;
        
        let x = lat.cos() * long.cos();
        let y = lat.cos() * long.sin();
        let z = lat.sin();
        
        Vector3d(x, y, z)
    }
}

impl From<Vector3d> for SurfacePoint {
    
    #[inline(always)]
    fn from(v: Vector3d) -> SurfacePoint {
        let Vector3d(x, y, z) = v.normalize();
        
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
    points: Vec<SurfacePoint>,
    velocities: Vec<Vector3d>,
    vertices_to_edges: Vec<Vec<usize>>
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
    
    pub fn add_particle(&mut self, position: SurfacePoint, attraction: f64, spread: f64) {
        self.particles.push(Particle::new(position, attraction, spread));
    }
    
    pub fn add_point(&mut self, point: SurfacePoint) -> usize {
        self.points.push(point);
        self.velocities.push(Vector3d(0.0, 0.0, 0.0));
        self.vertices_to_edges.push(vec![]);
        self.points.len() - 1
    }
    
    pub fn add_edge(&mut self, p1: usize, p2: usize) {
        let (p1, p2) =
            if p1 < p2 { (p1, p2) }
            else       { (p2, p1) };
        
        self.vertices_to_edges[p1].push(p2);
        self.vertices_to_edges[p2].push(p1);
    }
    
    pub fn step(&mut self, time_step: f64) {
        for index in 0..self.points.len() {
            self.update_velocities(index, time_step);
        }
        
        self.move_points(time_step);
    }
    
    fn update_velocities(&mut self, index: usize, time_step: f64) {
        let here = Vector3d::from(self.points[index]);
        
        let mut force = Vector3d(0.0, 0.0, 0.0);
        
        // Effect of the particles
        for particle in self.particles.iter() {
            let there = Vector3d::from(particle.position);
            let k1 = here * there;
            let vector = there - k1 * here;
            let dist = k1.acos();
            let value = -particle.attraction *
                        (-dist * dist / (particle.spread * particle.spread * 2.0)).exp();
            force = force + value * vector / (particle.spread * particle.spread);
        }
        
        // Effect of neighbours
        for &neighbour in self.vertices_to_edges[index].iter() {
            let there = Vector3d::from(self.points[neighbour]);
            let k1 = here * there;
            let vector = there - k1 * here;
            let dist = k1.acos();
            let value = dist - PI / 2.0;
            force = force + value * vector;
        }
        
        self.velocities[index] = (self.velocities[index] + force * time_step) * 0.90;
    }
    
    fn move_points(&mut self, time_step: f64) {
        let it = Iterator::zip(self.points.iter_mut(), self.velocities.iter());
        for (pos, &vel) in it {
            *pos = SurfacePoint::from(Vector3d::from(*pos) + vel * time_step);
        }
    }
}

pub fn get_background(surface: &Surface, width: u32) -> Vec<[u8; 3]> {
    let height = width / 2;
    let mut buf = vec![0.0; (width * height) as usize];
    let mut min = f64::INFINITY;
    let mut max = f64::NEG_INFINITY;
    
    for y in 0..height {
        for x in 0..width {
            let long = x as f64 / width as f64 * 2.0 * PI - PI;
            let lat = PI / 2.0 - y as f64 / height as f64 * PI;
            let point = SurfacePoint::new(long, lat);
            
            // And then the corresponding x, y, z coordinates
            let here = Vector3d::from(point);
            
            // We need to get the combined effect of all particles sensed at this location
            let mut combined = 0.0;
            for particle in surface.particles.iter() {
                let there = Vector3d::from(particle.position);
                let dist = (here * there).acos();
                
                combined += particle.attraction *
                            (-dist * dist / (particle.spread * particle.spread * 2.0)).exp();
            }
            
            buf[(y * width + x) as usize] = combined;
            
            if combined > max { max = combined; }
            if combined < min { min = combined; }
        }
    }
    
    if -min > max {
        max = -min;
    }
    else {
        min = -max
    }
    
    buf.iter().map(|&val| {
        [((val - min) / (max - min) * 255.0) as u8; 3]
    }).collect()
}

fn wrap(val: u32, inc: i32, max: u32) -> u32 {
    if inc < 0 {
        let dec = (-inc) as u32;
        if dec > val {
            max - dec + val
        }
        else {
            val - dec
        }
    }
    else {
        let inc = inc as u32;
        if val + inc >= max {
            val + inc - max
        }
        else {
            val + inc
        }
    }
}

pub fn draw(surface: &Surface, background: &[[u8; 3]], filename: &str) {
    let width = (background.len() as f64 * 2.0).sqrt() as u32;
    let height = (width / 2) as u32;

    let mut imgbuf = image::ImageBuffer::new(width, height);
    
    // Draw the background
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        *pixel = image::Rgb { data: background[(y * width + x) as usize] };
    }
    
    // Draw the points
    for p in surface.points.iter() {
        let x = ((p.long + PI) / (2.0 * PI) * width as f64) as u32;
        let y = ((PI / 2.0 - p.lat) / PI * height as f64) as u32;
        
        for xi in -1..2 {
            let xi = wrap(x, xi, width);
            for yi in -1..2 {
                let yi = wrap(y, yi, height);
                *imgbuf.get_pixel_mut(xi, yi) = image::Rgb { data: [0u8, 0, 255] };
            }
        }
    }
    
    // Draw the lines
    // for (point_index, neighbours) in surface.vertices_to_edges.iter().enumerate() {
    //     let here = Vector3d::from(surface.points[point_index]);
    //     for &neighbour_index in neighbours.iter() {
    //         if neighbour_index < point_index { continue; }
            
    //         let there = Vector3d::from(surface.points[neighbour_index]);
    //         let r = here * there;
    //         let h = there - r * here;
            
    //         let b1 = here;
    //         let b2 = h.normalize();
    //         let final_angle = r.acos();
    //         let n = 200u32;
            
    //         let p = surface.points[point_index];
    //         let x = ((p.long + PI) / (2.0 * PI) * (width - 1) as f64) as u32;
    //         let y = ((PI / 2.0 - p.lat) / PI * (height - 1) as f64) as u32;
    //         let (mut prev_x, mut prev_y) = (x as i32, y as i32);
            
    //         for index in 1..(n + 1) {
    //             let a = (index as f64 / n as f64) * final_angle;
    //             let k = b1 * a.cos() + b2 * a.sin();
                
    //             let p = SurfacePoint::from(k);
    //             let x = ((p.long + PI) / (2.0 * PI) * (width - 1) as f64) as i32;
    //             let y = ((PI / 2.0 - p.lat) / PI * (height - 1) as f64) as i32;
                
    //             if (x - prev_x).abs() > (y - prev_y).abs() {
    //                 if prev_x < x {
    //                     for xi in prev_x..(x + 1) {
    //                         let slope = (y - prev_y) as f64 / (x - prev_x) as f64;
    //                         let yi = prev_y + (slope * (xi - prev_x) as f64) as i32;
    //                         *imgbuf.get_pixel_mut(xi as u32, yi as u32) =
    //                             image::Rgb { data: [255u8, 0, 0] };
    //                     }
    //                 }
    //             }
                
    //             prev_x = x;
    //             prev_y = y;
    //         }
    //     }
    // }
    
    let mut fout = &mut File::create(&Path::new(filename)).unwrap();
    image::ImageRgb8(imgbuf).save(fout, image::PNG).unwrap();
}
