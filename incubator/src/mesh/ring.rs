use std::collections::HashMap;

/// Represents a closed ring of vertices in a `Mesh`.
///
/// A ring is a sequence of vertices that form a closed loop. For example, a `Mesh` can have a ring
/// defined by the sequence `1 → 4 → 9 → 3`, as long as there is an edge from `1` to `4`, etc. The
/// ring divides a `Mesh` in two regions: an inside and an outside. The inside of the ring is
/// defined in such a way that, if `a` is a vertex in the ring and `b` is the next vertex in the
/// ring, then the inside of the ring contains the vertex that is to the left of the edge `a → b`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Ring {
    inner: HashMap<usize, RingData>,
}

#[derive(Clone, PartialEq, Eq, Debug)]
struct RingData {
    prev: usize,
    next: usize,
}

impl Ring {
    /// Creates a `Ring` from a slice. This requires that the slice contains at least two elements,
    /// and that all elements are distinct.
    pub fn new(slice: &[usize]) -> Ring {
        assert!(slice.len() >= 2);
        
        let n = slice.len();
        let window = {
            let mut window = vec![];
            window.push(slice[n - 1]);
            window.extend(slice);
            window.push(slice[0]);
            window
        };
        
        let prev = window[0..n].iter();
        let this = window[1..n + 1].iter();
        let next = window[2..n + 2].iter();
        
        // Zip through the given slice by forming pairs of adjacent vertices in the slice
        Ring { inner:
            izip!(prev, this, next)
            .map(|(&prev, &this, &next)| { (this, RingData { prev: prev, next: next, }) })
            .collect(),
        }
    }
    
    /// Determines the size of this ring
    pub fn len(&self) -> usize {
        self.inner.len()
    }
    
    /// Gets a random vertex from the ring
    pub fn get_random_vertex(&self) -> usize {
        *self.inner.keys().take(1).next().unwrap()
    }
    
    /// Given a vertex in the ring, returns the next vertex
    #[inline(always)]
    pub fn next(&self, vertex: usize) -> usize {
        self.inner[&vertex].next
    }
    
    /// Given a vertex in the ring, returns the next vertex
    #[inline(always)]
    pub fn prev(&self, vertex: usize) -> usize {
        self.inner[&vertex].prev
    }
    
    /// Determines whether the ring contains the specified vertex
    pub fn contains_vertex(&self, v: usize) -> bool {
        self.inner.contains_key(&v)
    }
    
    /// Determines whether the two given vertices are part of the ring and whether they are
    /// neighbors on the ring.
    pub fn has_edge(&self, a: usize, b: usize) -> bool {
        match self.inner.get(&a) {
            None => false,
            Some(&ref data) => data.next == b || data.prev == b,
        }
    }
    
    /// Reroutes the ring so that an existing edge is split into two. This method assumes that the
    /// ring contains the `src` vertex and that it does not contain the `dest` vertex.
    pub fn insert_edge(&mut self, src: usize, dest: usize) {
        // assert!(self.inner.contains_key(&src));
        // assert!(!self.inner.contains_key(&dest));
        
        let next = self.next(src);
        
        self.inner.get_mut(&src).unwrap().next = dest;
        self.inner.get_mut(&next).unwrap().prev = dest;
        
        self.inner.insert(dest, RingData { prev: src, next: next });
    }
    
    /// Removes a vertex from the ring. This method assumes that the size of the ring is at least 3
    /// and that the ring contains the given vertex
    pub fn remove_vertex(&mut self, vertex: usize) {
        // assert!(self.len() >= 3);
        // assert!(self.inner.contains_key(&vertex));
        
        let prev = self.prev(vertex);
        let next = self.next(vertex);
        
        self.inner.get_mut(&prev).unwrap().next = next;
        self.inner.get_mut(&next).unwrap().prev = prev;
        
        self.inner.remove(&vertex);
    }
    
    /// Returns an iterator that yields the vertices of the ring, starting in a random vertex and
    /// iterating over the ring in order
    pub fn vertices<'a>(&'a self) -> impl Iterator<Item = usize> + 'a {
        RingVertexIter::new(self)
    }
    
    /// Returns an iterator that yields the edges of the ring, in a random order.
    pub fn edges<'a>(&'a self) -> impl Iterator<Item = (usize, usize)> + 'a {
        self.inner.iter().map(|(&x, &RingData { prev: _, next })| (x, next))
    }
    
}

struct RingVertexIter<'a> {
    current: usize,
    start: usize,
    done: bool,
    ring: &'a Ring
}

impl<'a> RingVertexIter<'a> {
    fn new(ring: &'a Ring) -> RingVertexIter<'a> {
        let current = ring.get_random_vertex();
        
        RingVertexIter {
            ring: ring,
            current: current,
            start: current,
            done: false,
        }
    }
}

impl<'a> Iterator for RingVertexIter<'a> {
    type Item = usize;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.done { return None; }
        
        let result = self.current;
        self.current = self.ring.next(result);
        
        if self.current == self.start {
            self.done = true;
        }
        
        Some(result)
    }
}

// /// Creates a vector containing the vertices of the ring, starting with a random vertex. The vertex
// /// next to the first is placed in the second position in the vector, the vertex next to the second
// /// in the third position, and so on.
// impl<'a> From<&'a Ring> for Vec<usize> {
//     fn from(ring: &Ring) -> Vec<usize> {
//         let sentinel = ring.get_random_vertex();
        
//         let mut result = Vec::with_capacity(ring.len());
//         let mut current = sentinel;
//         while ring.next(current) != sentinel {
//             result.push(current);
//             current = ring.next(current);
//         }
        
//         result.push(current);
//         result
//     }
// }

// /// Creates a Ring from a slice of `usize` integers, where each integer in the slice is mapped into
// /// the next one, and the last integer is mapped to the first one.
// ///
// /// This conversion assumes the given slice contains at least two elements, and that all elements
// /// are distinct.
// impl<'a> From<&'a [usize]> for Ring {
//     fn from(slice: &[usize]) -> Ring {
//         Ring::new(slice)
//     }
// }

#[cfg(test)]
mod test {
    use super::*;
    
    #[test]
    #[should_panic]
    fn test_new_ring_empty() {
        Ring::new(&[]);
    }
    
    #[test]
    #[should_panic]
    fn test_new_ring_size_1() {
        Ring::new(&[1]);
    }
    
    #[test]
    fn test_slice_to_ring_small() {
        let ring = Ring::new(&[1, 2]);
        
        assert_eq!(ring.next(1), 2);
        assert_eq!(ring.next(2), 1);
        
        assert_eq!(ring.prev(1), 2);
        assert_eq!(ring.prev(2), 1);
    }
    
    #[test]
    fn test_slice_to_ring() {
        let ring = Ring::new(&[1, 2, 3, 4, 5]);
        
        assert_eq!(ring.next(1), 2);
        assert_eq!(ring.next(2), 3);
        assert_eq!(ring.next(3), 4);
        assert_eq!(ring.next(4), 5);
        assert_eq!(ring.next(5), 1);
        
        assert_eq!(ring.prev(1), 5);
        assert_eq!(ring.prev(2), 1);
        assert_eq!(ring.prev(3), 2);
        assert_eq!(ring.prev(4), 3);
        assert_eq!(ring.prev(5), 4);
    }
    
    #[test]
    #[should_panic]
    fn test_ring_next_on_non_existing_vertex() {
        let ring = Ring::new(&[1, 2, 3, 4, 5]);
        ring.next(10);
    }
    
    fn shift_slice(slice: &mut [usize], start: usize) {
        
        let shift = match slice.iter().position(|&el| el == start) {
            Some(idx) => idx,
            None => return,
        };
        
        let new_vec = {
            let iter = slice.iter().skip(shift).chain(slice.iter().take(shift));
            iter.cloned().collect::<Vec<_>>()
        };
        
        for (idx, &i) in new_vec.iter().enumerate() {
            slice[idx] = i;
        }
    }
    
    fn test_ring_to_slice_unit(slice: &[usize]) {
        let ring: Ring = Ring::new(slice);
        let mut vec = ring.vertices().collect::<Vec<_>>();
        shift_slice(&mut vec, slice[0]);
        assert_eq!(vec, slice);
    }
    
    #[test]
    fn test_ring_to_slice() {
        test_ring_to_slice_unit(&[1, 2]);
        test_ring_to_slice_unit(&[1, 2, 3]);
        test_ring_to_slice_unit(&[1, 2, 3, 4]);
        test_ring_to_slice_unit(&[1, 2, 3, 4, 5]);
    }
}
