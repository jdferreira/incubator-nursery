use std::ops::{Add,Sub,Mul,Div};

#[derive(Debug,Copy,Clone)]
pub struct Vector3D(pub f64, pub f64, pub f64);

impl Vector3D {
    #[inline(always)]
    pub fn norm(self) -> f64 {
        self.norm_square().sqrt()
    }
    
    #[inline(always)]
    pub fn norm_square(self) -> f64 {
        self.0 * self.0 + self.1 * self.1 + self.2 * self.2
    }
    
    #[inline(always)]
    pub fn normalize(self) -> Vector3D {
        self / self.norm()
    }
}

impl Add for Vector3D {
    type Output = Vector3D;
    
    #[inline(always)]
    fn add(self, other: Vector3D) -> Vector3D {
        Vector3D(self.0 + other.0, self.1 + other.1, self.2 + other.2)
    }
}

impl Sub for Vector3D {
    type Output = Vector3D;
    
    #[inline(always)]
    fn sub(self, other: Vector3D) -> Vector3D {
        Vector3D(self.0 - other.0, self.1 - other.1, self.2 - other.2)
    }
}

impl Mul<f64> for Vector3D {
    type Output = Vector3D;
    
    #[inline(always)]
    fn mul(self, other: f64) -> Vector3D {
        Vector3D(self.0 * other, self.1 * other, self.2 * other)
    }
}

impl Mul<Vector3D> for f64 {
    type Output = Vector3D;
    
    #[inline(always)]
    fn mul(self, other: Vector3D) -> Vector3D {
        Vector3D(self * other.0, self * other.1, self * other.2)
    }
}

impl Mul<Vector3D> for Vector3D {
    type Output = f64;
    
    #[inline(always)]
    fn mul(self, other: Vector3D) -> f64 {
        self.0 * other.0 + self.1 * other.1 + self.2 * other.2
    }
}

impl Div<f64> for Vector3D {
    type Output = Vector3D;
    
    #[inline(always)]
    fn div(self, other: f64) -> Vector3D {
        Vector3D(self.0 / other, self.1 / other, self.2 / other)
    }
}

