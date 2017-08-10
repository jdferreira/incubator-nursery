use rand_map::{NBITS,PER_GROUP,MAP};

pub const OVER: usize = 4;
pub const BIT_LENGTH: usize = NBITS + OVER;
pub const SECONDARY_MAX: f64 = (1 << OVER) as f64;


#[inline(always)]
fn classic(seq: &[u8]) -> usize {
    seq.iter().rev().fold(0, |acc, &b| acc * 2 + (b as u32)) as usize
}


/// Returns the value associated with a bit sequence
pub fn get_value(seq: &[u8]) -> f64 {
    assert!(seq.len() == BIT_LENGTH);
    let (primary, secondary) = seq.split_at(NBITS);
    
    let primary_index = {
        let mut index = 0;
        for chunk in primary.chunks(PER_GROUP) {
            index += classic(chunk);
        }
        index
    };
    let prev = if primary_index == 0 { 0.0 } else { MAP[primary_index - 1] };
    let next = MAP[primary_index];
    
    let weight = classic(secondary) as f64 / SECONDARY_MAX as f64;
    prev * weight + next * (1.0 - weight)
}
