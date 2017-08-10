use std::collections::HashMap;

fn main() {
    let mut x: HashMap<usize, usize> = HashMap::new();
    x.insert(10, 20);
    x.insert(11, 21);
    x.insert(12, 22);
    x.insert(12, 0);
    
    println!("{:?}", x);
    println!("{:?}", x[&10]);
}
