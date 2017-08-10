// use std::f64;
// use surface::Vector3d;

// #[inline(always)]
// pub fn to_cartesian(long: f64, lat: f64) -> Vector3d {
//     let x = lat.cos() * long.cos();
//     let y = lat.cos() * long.sin();
//     let z = lat.sin();
//     Vector3d(x, y, z)
// }

// #[inline(always)]
// pub fn to_rad(val: f64) -> f64 {
//     val * f64::consts::PI / 180.0
// }

// #[inline(always)]
// pub fn to_degree(val: f64) -> f64 {
//     val * 180.0 / f64::consts::PI
// }
