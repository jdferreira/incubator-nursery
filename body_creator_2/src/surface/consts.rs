/// The points on the surface respond to the particles. This constant
/// defines how strong this responsiveness is
pub const PARTICLE_STRENGTH: f64 = 0.1;

/// Number of points in the sphere
pub const NUM_PARTICLES: u32 = 20;

/// This constant defines the strength governing the movement of a point
/// towards its neighbours
pub const NEIGHBORS_STRENGTH: f64 = 0.0002;

/// Points experience drag, which reduces their velocity
pub const DRAG: f64 = 0.8;

/// Time step
pub const TIME_STEP: f64 = 0.2;

/// Number of points in the sphere
pub const NUM_POINTS: u32 = 20;

/// We also want each particle to have at least some spread factor
pub const MIN_SPREAD: f64 = 0.3;
