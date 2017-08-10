/// We want the field to be calculated with a certain number of divisions.
/// If more accuracy is necessary, this can be increased.
pub const FIELD_DIVISION: u32 = 20;

/// These number of attraction/repulsion particles to create in the field
pub const NUM_PARTICLES: u32 = 150;

/// We also want each particle to have at least some spread factor
pub const MIN_SPREAD: f64 = 0.1;

/// This decides how many divisions on the starting embryo should be
/// calculated. 0 means use the simple icosahedron; higher values lead to
/// further subdivision steps
pub const EMBRYO_DIVISIONS: u32 = 3;

/// The distance from the center of the field, which is where the vertices
/// are initialized
pub const RADIUS: f64 = 0.05;

/// The vertices respond to the field. This constant defines how strong
/// this responsiveness is
pub const FIELD_STRENGTH: f64 = 0.125;

/// The vertices also respond to the presence of other vertices nearby.
/// This constant defines how strong that responsiveness is
pub const PROXIMITY_STRENGTH: f64 = 0.00026875;

/// This constant defines the strength governing the movement of a vertex
/// towards the point in between its neighbours
pub const NEIGHBOURS_STRENGTH: f64 = 30.0;

/// The maximum distance between a vertex and the field's edge to activate
/// the edge repulsion
pub const SYMMETRY_EDGE_THICKNESS: f64 = 0.025;

/// Vertices experience drag, which reduces their velocity
pub const DRAG: f64 = 0.8;

/// The maximum allowed velocity in each time frame
pub const MAX_VELOCITY: f64 = 0.05;

/// Time step
pub const TIME_STEP: f64 = 0.1;

/// The minimum spread for the central particle in the body
pub const CENTRAL_SPREAD: f64 = MIN_SPREAD * 0.2;
