#![allow(dead_code)]
#![feature(slice_patterns)]

mod body_symmetric;
mod field;
mod genome;
mod mesh;
mod rand_get;
mod rand_map;
mod surface;
mod utils;
mod vector;

use std::io::Write;
use vector::Vector3D;


macro_rules! stderr (
    ($($arg:tt)*) => {
        match writeln!(&mut ::std::io::stderr(), $($arg)*) {
            Ok(_) => (),
            Err(err) => panic!("Unable to write to stderr: {}", err),
        }
    }
);

fn main0() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    
    if args.len() != 3 {
        return Err("USAGE: body_creator GENOME ITERATIONS".to_string());
    }
    
    let genome_filename = &args[1];
    let iterations = match args[2].parse::<u32>() {
        Ok(iterations) => iterations,
        Err(err) => {
            return Err(format!("{}", err));
        }
    };
    
    let mut my_genome = match genome::Genome::new(genome_filename) {
        Ok(my_genome) => my_genome,
        Err(err) => {
            return Err(format!("{}", err));
        }
    };
    
    let mut body = body_symmetric::Body::new();
    body.attach_genome(&mut my_genome);
    
    for &face in body.mesh.faces.iter() {
        println!("{:?}", face);
    }
    println!("///");
    
    for &Vector3D(x, y, z) in body.mesh.positions.iter() {
        println!("{:.8} {:.8} {:.8}", x, y, z);
    }
    println!("///");
    
    for i in 0..iterations {
        stderr!("{}", i);
        body.step();
        for &Vector3D(x, y, z) in body.mesh.positions.iter() {
            println!("{:.8} {:.8} {:.8}", x, y, z);
        }
        println!("///");
    }
    
    Ok(())
}

fn main1() -> Result<(), String> {
    let args: Vec<_> = std::env::args().collect();
    
    if args.len() != 3 {
        return Err("USAGE: body_creator GENOME ITERATIONS".to_string());
    }
    
    let genome_filename = &args[1];
    let iterations = match args[2].parse::<u32>() {
        Ok(iterations) => iterations,
        Err(err) => {
            return Err(format!("{}", err));
        }
    };
    
    let mut my_genome = match genome::Genome::new(genome_filename) {
        Ok(my_genome) => my_genome,
        Err(err) => {
            return Err(format!("{}", err));
        }
    };
    
    let mut s = surface::Surface::new();
    s.attach_genome(&mut my_genome);
    
    for neighbours in s.vertices_to_edges.iter() {
        println!("{:?}", neighbours);
    }
    println!("///");

    for &point in s.points.iter() {
        let Vector3D(x, y, z) = Vector3D::from(point);
        println!("{:.8} {:.8} {:.8}", x, y, z);
    }
    println!("///");
    
    for _ in 0..iterations {
        s.step();
        for &point in s.points.iter() {
            let Vector3D(x, y, z) = Vector3D::from(point);
            println!("{:.8} {:.8} {:.8}", x, y, z);
        }
        println!("///");
    }
    
    Ok(())
}

fn main() {
    match main1() {
        Ok(()) => (),
        Err(err) => {
            stderr!("{:?}", err);
            std::process::exit(1);
        }
    }
}
