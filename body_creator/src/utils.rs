pub mod cdf {
    use std::f64;
    
    pub fn cdf(p: f64) -> f64 {
        return 0.5 * (1.0 + erf(p / f64::consts::SQRT_2))
    }
    
    /// Adapted from
    /// http://home.online.no/~pjacklam/notes/invnorm/impl/field/ltqnorm.txt
    pub fn inverse_cdf(p: f64) -> f64 {
        const A0: f64 = -3.969683028665376e+01;
        const A1: f64 =  2.209460984245205e+02;
        const A2: f64 = -2.759285104469687e+02;
        const A3: f64 =  1.383577518672690e+02;
        const A4: f64 = -3.066479806614716e+01;
        const A5: f64 =  2.506628277459239e+00;
        
        const B0: f64 = -5.447609879822406e+01;
        const B1: f64 =  1.615858368580409e+02;
        const B2: f64 = -1.556989798598866e+02;
        const B3: f64 =  6.680131188771972e+01;
        const B4: f64 = -1.328068155288572e+01;
        
        const C0: f64 = -7.784894002430293e-03;
        const C1: f64 = -3.223964580411365e-01;
        const C2: f64 = -2.400758277161838e+00;
        const C3: f64 = -2.549732539343734e+00;
        const C4: f64 =  4.374664141464968e+00;
        const C5: f64 =  2.938163982698783e+00;
        
        const D0: f64 =  7.784695709041462e-03;
        const D1: f64 =  3.224671290700398e-01;
        const D2: f64 =  2.445134137142996e+00;
        const D3: f64 =  3.754408661907416e+00;
        
        const LOW_THRESHOLD: f64 = 0.02425;
        const HIGH_THRESHOLD: f64 = 1.0 - LOW_THRESHOLD;
        
        if p < LOW_THRESHOLD || p > HIGH_THRESHOLD {
            let q = if p < LOW_THRESHOLD {
                // Rational approximation for lower region
                (-2.0 * p.ln()).sqrt()
            }
            else {
                (-2.0 * (1.0 - p).ln()).sqrt()
            };
            
            let num = ((((C0 * q + C1) * q + C2) * q + C3) * q + C4) * q + C5;
            let den = (((D0 * q + D1) * q + D2) * q + D3) * q + 1.0;
            
            if p < LOW_THRESHOLD {
                num / den
            }
            else {
                -num / den
            }
        }
        else {
            // Rational approximation for central region
            let q = p - 0.5;
            let r = q * q;
            
            let num = ((((A0 * r + A1) * r + A2) * r + A3) * r + A4) * r + A5;
            let den = ((((B0 * r + B1) * r + B2) * r + B3) * r + B4) * r + 1.0;
            
            q * num / den
        }
    }
    
    
    /// Adapted from http://www.johndcook.com/blog/cpp_erf/
    fn erf(x: f64) -> f64 {
        // constants
        const A1: f64 =  0.254829592;
        const A2: f64 = -0.284496736;
        const A3: f64 =  1.421413741;
        const A4: f64 = -1.453152027;
        const A5: f64 =  1.061405429;
        const P:  f64 =  0.3275911;
     
        // Save the sign of x
        let (sign, x) = if x < 0.0 { (-1.0, -x) } else { (1.0, x) };
        
        // A & S formula 7.1.26
        let t = 1.0 / (1.0 + P * x);
        let y = 1.0 - (((((A5 * t + A4) * t) + A3) * t + A2) * t + A1) *
                      t * (-x * x).exp();
        
        return sign * y;
    }
}

pub mod math {
    use std::ops::{Add,Mul};
    use std::f64;
    
    #[inline(always)]
    pub fn interpolate<T>(start: T, end: T, fraction: f64) -> T
    where T: Add<T, Output=T> + Mul<f64, Output=T> {
        start * (1.0 - fraction) + end * fraction
    }
    
    /// Converts a floating point number in the open interval (0, 1) to a normal
    /// distribution percentile, with an interval.
    pub fn convert_to_gaussian(f: f64, mean: f64, sigma: f64, start: f64, end: f64) -> f64 {
        use utils::cdf::{cdf, inverse_cdf};
        
        let start =
            if start == f64::NEG_INFINITY { 0.0 }
            else { cdf((start - mean) / sigma) };
        let end =
            if end == f64::INFINITY { 1.0 }
            else { cdf((end - mean) / sigma) };
        
        // Now, we map the given float (which is between 0 and 1, excluding
        // ends) into the interval (start, end)
        let f = interpolate(start, end, f);
        
        // We then use the inverse_cdf of the standard normal distribution
        // to convert the `f` (which comes from a uniform distribution)
        // into another value (from the standard normal distribution)
        let normal = inverse_cdf(f);

        // Finally, convert the standard distribution to the one with
        // the given parameters (mean and standard deviation)
        normal * sigma + mean
    }

}
