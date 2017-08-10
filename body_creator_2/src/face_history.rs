
/// Contains the history of a face, which is needed to properly do UV unwrapping of the whole
/// mesh. `history` is a vector that contains events in the history f this face. In particular,
/// it describes, for each iteration, the size of the triangle in the mesh (width and height),
/// the relative position of the vertices so that the triangle lays flat on its largest side,
/// the actual position of the vertices in the local UV space of this triangle, and whether the
/// face exists or is absent in each iteration.
///
/// See also `FaceStep`.
pub struct FaceHistory {
    history: Vec<FaceStep>,
    current_corner: usize,
    current_removed: bool,
}

/// The events that can happen in the history of a face.
/// A face history starts implicitly with an `Added` event on iteration 0, unless the first event
/// is actually an `Added` event. When iterating through a `FaceHistory`, the iteration number is
/// only changed with the `Sizes` and `Added` events.
pub enum FaceStep {
    /// The orientation of the face changes. The leftmost vertex is now the one specified by the
    /// `corner` field; since the faces are always stored in counter-clockwise fashion, this is
    /// enough to identify which positions the two other vertices of the face take (rightmost and
    /// top)
    Orientation {
        corner: usize,
    },
    /// Each iteration needs to know the UV coordinates of the vertices in the triangle. Since the
    /// leftmost vertex always lives in the position (0, 0), and the v-coordinate of the rightmost
    /// vertex is also 0, we need only to store the `width` (u-coordinate of the rightmost vertex),
    /// `height` (v-coordinate of the top vertex) and the `mid` (u-coordinate of the top vertex).
    /// Each `Sizes` event represents a different iteration, and each iteration that the face
    /// exists must have a `Sizes` event.
    Sizes {
        width: f64,
        height: f64,
        mid: f64
    },
    /// The face was removed from the mesh
    Removed,
    /// The face was added (or re-added) to the mesh on iteration `iteration`.
    Added {
        iteration: usize,
    },
}

impl FaceHistory {
    #[inline(always)]
    pub fn new() -> FaceHistory {
        FaceHistory::new_on_iteration(0)
    }
    
    pub fn new_on_iteration(iteration: usize) -> FaceHistory {
        FaceHistory {
            history: vec![FaceStep::Added { iteration: iteration }],
            current_corner: 0,
            current_removed: false
        }
    }
    
    pub fn add_step(&mut self, iteration: usize, corner: usize, width: f64, height: f64, mid: f64) {
        if self.current_removed {
            self.history.push(FaceStep::Added { iteration: iteration });
        }
        self.current_removed = false;
        
        if self.history.len() == 0 || self.current_corner != corner {
            self.history.push(FaceStep::Orientation{ corner: corner });
        }
        self.current_corner = corner;
        
        self.history.push(FaceStep::Sizes { width: width, height: height, mid: mid });
    }
    
    pub fn remove(&mut self) {
        self.history.push(FaceStep::Removed);
        self.current_removed = true;
    }
    
    pub fn latest_event(&self) -> &FaceStep {
        &self.history[self.history.len() - 1]
    }
}
