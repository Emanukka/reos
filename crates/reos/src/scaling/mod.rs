pub mod dehlouz;

const KB:f64 = 1.380649e-23;
const NAV:f64 = 6.02214076e23;


pub fn rosenfeld_viscosity(t:f64, d:f64, mw:f64) -> f64 {

    let m = mw / NAV / 1000.;
    
    (d * NAV).powf(2./3.) * (m * KB * t).sqrt()

}