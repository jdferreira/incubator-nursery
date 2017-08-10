use super::*;
// use std::collections::hash_map::Entry;

// trait ExtrusionAlgorithm {
//     type Result: ExtrusionResult;
    
//     fn hook<'a>(&'a ExtrusionData);
//     fn rewire<'a>(&'a ExtrusionData);
// }

// trait ExtrusionResult {
    
// }

pub struct ExtrusionData<'a> {
    mesh: &'a mut Mesh,
    ring: &'a Ring,
    dir: Vector3D,
    inside_vertices: Vec<usize>,
    extrud_vertices: HashMap<usize, usize>,
    middle_vertices: HashMap<usize, usize>,
}

impl <'a> ExtrusionData<'a> {
    fn rewire(&mut self) {
        
        for v in self.ring.vertices() {
            let mut map_clone = self.mesh.neighbours.get_mut(v).unwrap().clone();
            
            for (&n, &mut EdgeData { left: ref mut l, right: ref mut r }) in &mut map_clone {
                if !self.inside_vertices.contains(&n) { continue; }
                
                if self.ring.contains_vertex(n) && self.ring.has_edge(v, n) {
                    let v_lifted = self.get_lifted(v);
                    let n_lifted = self.get_lifted(n);
                    
                    if self.ring.next(v) == n {
                        // clone edge (v → n)
                        let l_lifted = self.get_lifted(*l);
                        let r_lifted = self.middle_vertices[&*l];
                        
                        let ref mut map = self.mesh.neighbours.get_mut(v_lifted).unwrap();
                        map.insert(n_lifted, EdgeData { left: l_lifted, right: r_lifted });
                        
                        // change the left of this edge
                        *l = self.middle_vertices[&v];
                    }
                    else {
                        // clone edge (v → n)
                        let l_lifted = self.middle_vertices[&*r];
                        let r_lifted = self.get_lifted(*r);
                        
                        let ref mut map = self.mesh.neighbours.get_mut(v_lifted).unwrap();
                        map.insert(n_lifted, EdgeData { left: l_lifted, right: r_lifted });
                        
                        // change the right of this edge
                        *l = self.middle_vertices[&n];
                    }
                }
            }
        }
    }
    
    #[inline(always)]
    fn get_lifted(&self, vertex: usize) -> usize {
        match self.extrud_vertices.get(&vertex) {
            None => vertex,
            Some(&lifted) => lifted,
        }
    }
}

pub fn extrude<'a>(mesh: &'a mut Mesh, ring: &'a Ring, dir: Vector3D) -> ExtrusionData<'a> {
    
    let extrud_vertices = make_extrud_vertices(mesh, ring, dir);
    let middle_vertices = make_middle_vertices(mesh, ring, dir);
    let inside_vertices = mesh.get_inside_vertices(ring);
    
    translate_inside_vertices(mesh, &inside_vertices, dir);
    
    let mut data = ExtrusionData {
        mesh: mesh,
        ring: ring,
        dir: dir,
        inside_vertices: inside_vertices,
        extrud_vertices: extrud_vertices,
        middle_vertices: middle_vertices,
    };
    
    data.rewire();
    data
}

fn make_extrud_vertices(mesh: &mut Mesh, ring: &Ring, dir: Vector3D) -> HashMap<usize, usize> {
    let mut new_vertices = HashMap::new();
    
    for v in ring.vertices() {
        let new_pos = mesh.positions[v] + dir;
        let new_vertex_idx = mesh.positions.len();
        
        new_vertices.insert(v, new_vertex_idx);
        
        mesh.positions.push(new_pos);
        mesh.neighbours.push(HashMap::new());
    }
    
    new_vertices
}

fn make_middle_vertices(mesh: &mut Mesh, ring: &Ring, dir: Vector3D) -> HashMap<usize, usize> {
    let mut new_vertices = HashMap::new();
    
    for (v1, v2) in ring.edges() {
        let new_pos = (mesh.positions[v1] + mesh.positions[v2]) / 2.0 + dir / 2.0;
        let new_vertex_idx = mesh.positions.len();
        
        new_vertices.insert(v1, new_vertex_idx);
        
        mesh.positions.push(new_pos);
        mesh.neighbours.push(HashMap::new());
    }
    
    new_vertices
}

fn translate_inside_vertices(mesh: &mut Mesh, inside_vertices: &[usize], dir: Vector3D) {
    for &idx in inside_vertices {
        mesh.positions[idx] += dir;
    }
}
