extern crate image;

use std::fs::File;
use std::path::Path;
use std::f64;

/// Notice: Latitude and longitude given in radians
fn to_xyz(lat: f64, long: f64) -> (f64, f64, f64) {
    let x = lat.cos() * long.cos();
    let y = lat.cos() * long.sin();
    let z = lat.sin();
    (x, y, z)
}

fn to_lat_long(x: f64, y: f64, z: f64) -> (f64, f64) {
    let lat = z.asin();
    let long = y.atan2(x);
    (lat, long)
}

fn to_cartesian(lat: f64, long: f64) -> (f64, f64) {
    (long, 2.0 * lat.sin())
}

fn main() {
    let pi = f64::consts::PI;
    
    // let (lat, long) = (32.0, -9.0);
    // let (lat, long) = (lat * pi / 180.0, long * pi / 180.0);
    // println!("{} {}", lat, long);
    
    // let (x, y, z) = to_xyz(lat, long);
    // println!("{} {} {}", x, y, z);
    // let (lat, long) = to_lat_long(x, y, z);
    // println!("{} {}", lat, long);
    // let (x, y) = to_cartesian(lat, long);
    // println!("{} {}", x, y);
    
    let scale = 200.0;
    
    let width = (scale * 2.0 * pi) as u32;
    let height = (scale * 4.0) as u32;

    let mut imgbuf = image::ImageBuffer::new(width, height);
    
    let pivot = (-9.0 * pi / 180.0, 32.0 * pi / 180.0); // Approx. location of Portugal on Earth
    let (xp, yp, zp) = to_xyz(pivot.0, pivot.1);
    
    // Iterate over the coordiantes and pixels of the image
    for (x, y, pixel) in imgbuf.enumerate_pixels_mut() {
        let lat = x as f64 / scale - pi; // [-pi, pi[
        let long = (1.0 - y as f64 / scale / 2.0).asin(); // ]-pi/2, pi/2]
        
        let (xm, ym, zm) = to_xyz(lat, long);
        let dist = (xp * xm + yp * ym + zp * zm).acos();
        
        let color = 255 - (dist / pi * 255.0) as u8;
        *pixel = image::Rgb::<u8> { data: [color; 3] };
    }
    
    let mut fout = &mut File::create(&Path::new("map.png")).unwrap();
    image::ImageRgb8(imgbuf).save(fout, image::PNG).unwrap();
}
