mod vector;
mod mesh;

use std::error::Error;
use std::fs::File;
use std::io::{self,BufRead,BufReader,Read};
use std::collections::{HashSet,HashMap};

use mesh::Mesh;
use vector::Vector3D;

const MAX_STRING_LENGTH : usize = 5;

fn get_mesh(path: &str) -> io::Result<Mesh> {
    let mut buf = String::new();
    
    let mut f = try!(File::open(path));
    try!(f.read_to_string(&mut buf));
    
    let parts: Vec<_> = buf.split("///\n").collect();
    
    let mut faces = HashSet::new();
    for line in parts[0].lines() {
        let fields: Vec<_> =
            line[1 .. line.len() - 1]
            .split(',')
            .map(|s| s.trim().parse::<usize>().unwrap())
            .collect();
        faces.insert((fields[0], fields[1], fields[2]));
    }
    
    let mut positions = vec![];
    for line in parts[parts.len() - 2].lines() {
        let numbers: Vec<_> =
            line
            .split_whitespace()
            .map(|s| s.parse::<f64>().unwrap())
            .collect();
        positions.push(Vector3D(numbers[0], numbers[1], numbers[2]));
    }
    
    Ok(Mesh::new(positions, faces))
}

fn step1(m: &Mesh, t: f64) {
    let mut sharpness = HashMap::new();
    
    for (idx, map) in m.neighbours.iter().enumerate() {
        for (&n, &(left, right)) in map {
            // The vertex with index `idx` has neighbour `n` with
            // `left` and `right` has the third vertex on each face adjacent
            // to the edge from `idx` to `n`
            
            // We want to process each edge only once, and we do it
            // when `idx < n`
            if idx > n { continue; }
            
            let here = m.positions[idx];
            let there = m.positions[n];
            let left = m.positions[left];
            let right = m.positions[right];
            
            let edge = there - here;
            let nr = (right - here).cross(edge).normalize();
            let nl = (left  - here).cross(edge).normalize();
            
            let cos_a = nr * nl;
            sharpness.insert((idx, n), (cos_a + 1.0) / 2.0);
        }
    }
    
    let t = (t * sharpness.len() as f64) as usize;
    let sharp_edges: Vec<(usize, usize)> = {
        let mut tmp: Vec<_> =
            sharpness.iter().map(|(&(a, b), &s)| (s, a, b)).collect();
        tmp.sort_by(|a, b| b.partial_cmp(a).unwrap());
        tmp = tmp.iter().take(t).cloned().collect();
        tmp.iter().map(|&(s, a, b)| (a, b)).collect()
    };
    
    // The vertices that are neighbors to some feature and as such cannot be
    // part of one themselves
    let feature_neighbors = HashSet::new();
    
    for (a, b) in sharp_edges {
        // expand the feature curve
        for &(a, b) in [(a, b), (b, a)].into_iter() {
            let (best, _) = make_curve(b, &vec![a], MAX_STRING_LENGTH,
                                       &m, &feature_neighbors, &sharpness);
            println!("{:?}", best);
        }
    }
}

fn make_curve(vertex: usize, before: &[usize], max: usize,
              m: &Mesh, fns: &HashSet<usize>,
              sharpness: &HashMap<(usize, usize), f64>)
        -> (Vec<usize>, f64) {
    if max == 0 {
        return (vec![vertex], 0.0);
    }
    
    // Grab this vertex's neighbors that can still be used to continue the
    // feature
    let ns = m.neighbours[vertex].keys()
        .filter(|n| !fns.contains(&n) && !before.contains(&n));
    let mut before = before.to_vec();
    before.push(vertex);
    
    let mut max_sharpness = 0.0;
    let mut best = None;
    for &n in ns {
        // println!("Vertex {} is neighbor to {}", vertex, n);
        let (this, mut this_sharpness) =
            make_curve(n, &before, max - 1, m, fns, sharpness);
        
        let key = if n > vertex { (vertex, n) } else { (n, vertex) };
        this_sharpness += *sharpness.get(&key).unwrap();
        if this_sharpness > max_sharpness {
            let mut the_best = vec![vertex];
            the_best.extend_from_slice(&this);
            best = Some(the_best);
            max_sharpness = this_sharpness;
        }
    }
    
    match best {
        None => (vec![vertex], 0.0),
        Some(best) => (best, max_sharpness),
    }
}

fn main() {
    let args: Vec<_> = std::env::args().collect();
    let m = match get_mesh(&args[1]) {
        Ok(m) => m,
        Err(why) => panic!("{}", why),
    };
    let t = args[2].parse().unwrap();
    
    let s = step1(&m, t);
}
