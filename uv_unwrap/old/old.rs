fn get_mesh(path: &str) -> Mesh {
    
    let mut f = match File::open(path) {
        Ok(f) => BufReader::new(f),
        Err(why) => panic!("couldn't open {}: {}", path, Error::description(&why)),
    };
    
    
    let mut reading_triangles = true;
    let mut faces = HashSet::new();
    let mut max_vertex = 0;
    let mut positions = vec![];
    let mut pos_i = 0;
    
    for line in f.lines() {
        let mut line = match line {
            Ok(line) => line,
            Err(why) => panic!("couldn't read line {}", Error::description(&why)),
        };
        
        if reading_triangles {
            if line == "///" {
                reading_triangles = false;
                positions = vec![Vector3D(0.0, 0.0, 0.0); max_vertex + 1];
                pos_i = 0;
                continue;
            }
            
            let fields: Vec<_> =
                line[1 .. line.len() - 1]
                .split(',')
                .map(|s| s.trim().parse::<usize>().unwrap())
                .collect();
            faces.insert((fields[0], fields[1], fields[2]));
            
            let tmp = fields.iter().cloned().max().unwrap();
            if tmp > max_vertex { max_vertex = tmp; }
        }
        else {
            if line == "///" {
                positions = vec![Vector3D(0.0, 0.0, 0.0); max_vertex + 1];
                pos_i = 0;
                continue;
            }
            
            let numbers: Vec<_> = line.split(" ").map(|s| s.parse::<f64>().unwrap()).collect();
            positions[pos_i] = Vector3D(numbers[0], numbers[1], numbers[2]);
            pos_i += 1;
        }
    }
    
    Mesh::new(positions, faces)
}
