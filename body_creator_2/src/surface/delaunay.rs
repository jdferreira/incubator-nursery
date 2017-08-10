use super::SurfacePoint;
use std::collections::{HashSet,HashMap};
use vector::Vector3D;

struct Delaunay {
    points: Vec<Vector3D>,
    faces: HashSet<(usize, usize, usize)>,
}

impl Delaunay {
    fn new(points: &[SurfacePoint]) -> Delaunay {
        Delaunay {
            points: points.iter().cloned().map(Vector3D::from).collect(),
            faces: HashSet::new(),
        }
    }
    
    fn process(&mut self) {
        // Make the first triangle using the three first surface points
        self.faces.insert((0, 1, 2));
        self.faces.insert((0, 2, 1));
        
        // Keep inserting each new point sequentially
        for index in 3..self.points.len() {
            let new_position = self.points[index];
            let mut in_conflict = vec![];
            let mut ring_arena = RingArena::new();
            
            
            for &(i, j, k) in self.faces.iter() {
                let pi = self.points[i];
                let pj = self.points[j];
                let pk = self.points[k];
                
                let circumcenter = Vector3D::cross(pj - pi, pk - pi).normalize();
                let radius = (circumcenter * pi).acos();
                let d = (circumcenter * new_position).acos();
                
                if d < radius {
                    in_conflict.push((i, j, k));
                }
            }
            
            for face in in_conflict {
                self.faces.remove(&face);
                ring_arena.add_face(face);
            }
            
            for (&i, &j) in ring_arena.rings[0].iter() {
                self.faces.insert((i, j, index));
            }
        }
    }
}

struct RingArena {
    rings: Vec<HashMap<usize, usize>>,
}

impl RingArena {
    fn new() -> RingArena {
        RingArena { rings: vec![] }
    }
    
    fn add_face(&mut self, (i, j, k): (usize, usize, usize)) {
        let adjacents = self.get_adjacents(i, j, k);
        
        match adjacents {
            [None, None, None] => {
                let mut map = HashMap::new();
                map.insert(i, j);
                map.insert(j, k);
                map.insert(k, i);
                self.rings.push(map);
            },
            
            [Some(a), None, None] |
            [None, Some(a), None] |
            [None, None, Some(a)] => {
                let (ring_idx_a, start_a, end_a, other_a) = a;
                self.rings[ring_idx_a].insert(start_a, other_a);
                self.rings[ring_idx_a].insert(other_a, end_a);
            },
            
            [Some(a), Some(b), None] |
            [None, Some(a), Some(b)] |
            [Some(b), None, Some(a)] => {
                let (ring_idx_a, start_a, end_a, other_a) = a;
                let (ring_idx_b, _, _, _) = b;
                
                let mut curr = start_a;
                while curr != other_a {
                    let next = self.rings[ring_idx_b][&curr];
                    self.rings[ring_idx_a].insert(curr, next);
                    curr = next;
                }
                self.rings[ring_idx_a].insert(other_a, end_a);
                self.rings.remove(ring_idx_b);
            },
            
            [Some(a), Some(b), Some(c)] => {
                let (ring_idx_a, start_a, end_a, other_a) = a;
                let (ring_idx_b, _, _, _) = b;
                let (ring_idx_c, _, _, _) = c;
                
                let mut curr = start_a;
                while curr != other_a {
                    let next = self.rings[ring_idx_b][&curr];
                    self.rings[ring_idx_a].insert(curr, next);
                    curr = next;
                }
                while curr != end_a {
                    let next = self.rings[ring_idx_c][&curr];
                    self.rings[ring_idx_a].insert(curr, next);
                    curr = next;
                }
                
                self.rings.remove(ring_idx_b);
                self.rings.remove(ring_idx_c);
            }
        }
    }
    
    #[inline(always)]
    fn get_adjacents(&self, i: usize, j: usize, k: usize) -> [Option<(usize, usize, usize, usize)>; 3] {
        let mut result = [None, None, None];
        for (ring_idx, ring) in self.rings.iter().enumerate() {
            if ring.get(&j) == Some(&i) {
                result[0] = Some((ring_idx, j, i, k));
            }
            if ring.get(&k) == Some(&j) {
                result[1] = Some((ring_idx, k, j, i));
            }
            if ring.get(&i) == Some(&k) {
                result[2] = Some((ring_idx, i, k, j));
            }
        }
        result
    }
}

pub fn triangulate(points: &[SurfacePoint]) -> Vec<HashSet<usize>> {
    let mut d = Delaunay::new(points);
    d.process();
    
    let mut result = vec![HashSet::new(); points.len()];
    for &(i, j, k) in d.faces.iter() {
        result[i].insert(j);
        result[j].insert(i);
        result[i].insert(k);
        result[k].insert(i);
        result[j].insert(k);
        result[k].insert(j);
    }
    result
}






#[cfg(test)]
impl ::std::fmt::Display for RingArena {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        try!(write!(f, "RingArena {{"));
        for ring in self.rings.iter() {
            let start = ring.keys().cloned().min().unwrap();
            try!(write!(f, " {}", start));
            
            let mut curr = ring[&start];
            while curr != start {
                try!(write!(f, " -> {}", curr));
                curr = ring[&curr];
            }
            try!(write!(f, " -> {},", curr));
        }
        try!(write!(f, " }}"));
        Ok(())
    }
}

#[test]
fn test() {
    let mut r = RingArena::new();
    r.add_face((0, 1, 4));
    println!("{}", r);
    r.add_face((1, 2, 4));
    println!("{}\n", r);
    
    let mut r = RingArena::new();
    r.add_face((0, 1, 4));
    println!("{}", r);
    r.add_face((4, 2, 3));
    println!("{}\n", r);
    
    let mut r = RingArena::new();
    r.add_face((0, 1, 4));
    println!("{}", r);
    r.add_face((4, 2, 3));
    println!("{}", r);
    r.add_face((1, 2, 4));
    println!("{}\n", r);
    
    let mut r = RingArena::new();
    r.add_face((0, 1, 4));
    println!("{}", r);
    r.add_face((4, 2, 3));
    println!("{}", r);
    r.add_face((1, 5, 2));
    println!("{}", r);
    r.add_face((1, 2, 4));
    println!("{}\n", r);
    
}
