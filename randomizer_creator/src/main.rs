/// Returns a vector of integers, where the value in position `i` is the number
/// of distinct bit sequences of size `nbits` whose value, as interpreted by
/// the `bit_interpretion` function, is equal to `i`.
fn count(nbits: usize, per_group: usize) -> Vec<usize> {
    // First find the maximum sum possible so that we can initialize the vector
    // The maximum in each group is (2.pow(per_group) - 1)
    // The maximum overall is that multipled by the number of groups
    let max_sum: usize = {
        let nbits = nbits as usize;
        let per_group = per_group as usize;
        (nbits / per_group) * (2usize.pow(per_group as u32) - 1)
    };
    
    // Initialize the result with exactly `max_sum` 0s. They will be increased
    // in the next loop.
    let mut result = vec![0; max_sum + 1];
    
    for seq in BitIterator::new(nbits) {
        let value = bit_interpretion(seq, per_group);
        result[value] += 1;
    }
    
    return result;
}


/// Interpret the given bit sequence and assign it a value, based on the number
/// of bits in each sub-group. To do this, divide the bit sequence into groups
/// of size `per_group`, convert each group into an integer and return the sum
/// of these values
fn bit_interpretion(seq: Vec<u8>, per_group: usize) -> usize {
    let mut result = 0;
    for chunk in seq.chunks(per_group) {
        result += to_int(chunk);
    }
    result
}


/// Given a slice of bits, interpret them in the usual binary way. A catch here
/// is that the first bit in the slice represents the least significant bit
/// of the number, so to_int(&[0, 1, 1]) returns 6 instead of 3.
fn to_int(slice: &[u8]) -> usize {
    slice.iter().rev().fold(0, |acc, &b| acc * 2 + (b as u32)) as usize
}



/// This respresents an iterator that produces all the possible bit sequences
/// of a specified size.
struct BitIterator {
    current: Vec<u8>,
    untouched: bool
}

impl BitIterator {
    /// Create a BitIterator with `nbits` bits
    fn new(nbits: usize) -> BitIterator {
        BitIterator { current: vec![0; nbits], untouched: true }
    }
}

impl Iterator for BitIterator {
    type Item = Vec<u8>;
    
    /// Iterate over this BitIterator in order to find the next bit pattern
    fn next(&mut self) -> Option<Vec<u8>> {
        if self.untouched {
            self.untouched = false;
            return Some(self.current.clone());
        }
        
        let mut i = 0;
        while i < self.current.len() && self.current[i] == 1 {
            self.current[i] = 0;
            i += 1
        }
        
        if i == self.current.len() {
            None
        }
        else {
            self.current[i] = 1;
            Some(self.current.clone())
        }
    }
}


fn join(vec: &Vec<f64>, sep: &str) -> String {
    let mut iter = vec.iter();
    match iter.next() {
        Some(elem) => {
            let mut output = format!("{:.8}", elem);
            for elem in iter {
                output.push_str(sep);
                output.push_str(&format!("{:.8}", elem));
            }
            output
        }
        None => String::new()
    }
}


fn main() {
    const NBITS: usize = 12;
    const PER_GROUP: usize = 4;
    const TOTAL: u32 = 1 << NBITS;
    
    let count_sequence = count(NBITS, PER_GROUP);
    let mut cumm = 0u32;
    let map: Vec<_> = count_sequence.iter()
        .map(|&x| { cumm += x as u32; cumm as f64 / TOTAL as f64 })
        .collect();
    
    println!("pub const NBITS: usize = {};", NBITS);
    println!("pub const PER_GROUP: usize = {};", PER_GROUP);
    println!("pub const MAP: [f64; {}] = [\n\t{}];",
             map.len(),
             join(&map, ",\n\t"));
}
