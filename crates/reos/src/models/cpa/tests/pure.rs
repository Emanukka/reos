use approx::assert_relative_eq;
use ndarray::array;

use crate::{models::cpa::SCPA, residual::Residual};

pub fn water()->SCPA{
    let water = super::recipes::water4c();
    super::recipes::scpa(vec![water], vec![]).unwrap()
}

const TOL:f64 = 1e-12;
#[test]
fn test_scpa_helmholtz() {
    
    let model = water();
    let t = 298.15;
    let d= 58_000.;
    let x = array![1.0];
    let val = model.r_helmholtz(t, d, &x);
    assert_relative_eq!(val, -9.705428401736341 , epsilon = TOL)
}
#[test]
fn test_scpa_entropy() {
    
    let model = water();
    let t = 298.15;
    let d= 58_000.;
    let x = array![1.0];
    let val = model.r_entropy(t, d, &x);
    assert_relative_eq!(val, -6.954623177034006, epsilon = TOL)
}
#[test]
fn test_scpa_chem_pot() {
    
    let model = water();
    let t = 298.15;
    let d= 58_000.;
    let x = array![1.0];
    let val = model.r_chemical_potential(t, d, &x);
    assert_relative_eq!(val[0], -9.76898881609111 , epsilon = TOL)
}

#[test]
fn test_scpa_pressure() {
    
    let model = water();
    let t = 298.15;
    let d= 58_000.;
    let x = array![1.0];
    let val = model.r_pressure(t, d, &x) / d;
    
    
    assert_relative_eq!(val, (0.9364395856452267 - 1.) , epsilon = TOL)
}