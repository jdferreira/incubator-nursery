use std::fs::File;
use std::io::{Error, ErrorKind};
use std::io::Read;

use rand_get;

pub struct Genome {
    file: File,
    leftover: Vec<u8>
}

impl Genome {
    
    pub fn new(path: &str) -> Result<Genome, Error> {
        let file = try!(File::open(path));
        Ok(Genome { file: file, leftover: vec![] })
    }
    
    fn get_bits(&mut self, amount: usize) -> Result<Vec<u8>, Error> {
        let mut result: Vec<u8> = vec![];
        
        if self.leftover.len() > 0 {
            if amount < self.leftover.len() {
                let result: Vec<_> = self.leftover.iter().take(amount)
                                         .cloned().collect();
                self.leftover = self.leftover.iter().skip(amount)
                                             .cloned().collect();
                return Ok(result);
            }
            
            result = self.leftover.clone();
            self.leftover.clear();
        }
        
        let bytes_to_read = {
            let bits_to_read = amount - result.len();
            if bits_to_read % 8 == 0 {
                bits_to_read / 8
            }
            else {
                bits_to_read / 8 + 1
            }
        };
        
        let mut buf = vec![0; bytes_to_read];
        let bytes_read = try!(self.file.read(&mut buf));
        result.extend(decode(&buf));
                
        // If we read less than the required bytes, it means that the file
        // does not contain more bytes, and as such we return an error
        if bytes_read < bytes_to_read {
            self.leftover = result;
            Err(Error::new(ErrorKind::Other, "Not enough bits to read"))
        }
        // Otherwise, everything is peachy, and we can return the result. Notice
        // that this result can have more bits than requested, in which case the
        // excess must be transferred to the leftovers
        else {
            self.leftover = result.iter().skip(amount).cloned().collect();
            result = result.iter().take(amount).cloned().collect();
            Ok(result)
        }
    }
    
    pub fn get_float(&mut self) -> f64 {
        let bits = match self.get_bits(rand_get::BIT_LENGTH) {
            Ok(bits) => bits,
            Err(error) => {
                panic!(format!("ERROR: genome::get_float {}", error));
            }
        };
        
        // if bits.iter().all(|&x| x == 0) {
        //     bits[bits.len() - 1] = 1;
        // }
        let result = rand_get::get_value(&bits);
        if result < 0f64 || result > 1f64 {
            panic!("WARNING: genome::get_float returns {}", result);
        }
        else {
            result
        }
    }
}

fn decode(buf: &[u8]) -> Vec<u8> {
    let mut result = vec![0; buf.len() * 8];
    for buf_index in 0..buf.len() {
        let mut tmp = buf[buf_index];
        for inner_index in 0..8 {
            let inner_index = buf_index * 8 + 7 - inner_index;
            result[inner_index] = tmp % 2;
            tmp = tmp / 2;
        }
    }
    result
}
