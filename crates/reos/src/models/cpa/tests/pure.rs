use approx::assert_relative_eq;
use ndarray::array;

use crate::{models::cpa::SCPA, residual::Residual};

pub fn water()->SCPA{
    let water = super::recipes::water4c();
    super::recipes::scpa(vec![water], vec![]).unwrap()
}

#[test]
fn test_scpa_helmholtz() {
    
    let model = water();
    let t = 298.15;
    let d= 1_000.;
    let x = array![1.0];
    let val = model.r_helmholtz(t, d, &x);
    assert_relative_eq!(val, -0.058144295861 + -1.597921023379 , epsilon = 1e-9)
}
#[test]
fn test_scpa_entropy() {
    
    let model = water();
    let t = 298.15;
    let d= 1_000.;
    let x = array![1.0];
    let val = model.r_entropy(t, d, &x);
    assert_relative_eq!(val, -0.041951593945 + -4.713659269705, epsilon = 1e-9)
}
#[test]
fn test_scpa_chem_pot() {
    
    let model = water();
    let t = 298.15;
    let d= 1_000.;
    let x = array![1.0];
    let val = model.r_chemical_potential(t, d, &x);
    assert_relative_eq!(val[0], -0.115660251059 + -2.54386196979185 , epsilon = 1e-10)
}

#[test]
fn test_scpa_pressure() {
    
    let model = water();
    let t = 298.15;
    let d= 1_000.;
    let x = array![1.0];
    let val = model.r_pressure(t, d, &x);
    assert_relative_eq!(val, -57.5159551979349 + -945.9409464127781, epsilon = 1e-8)
}