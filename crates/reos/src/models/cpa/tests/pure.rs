use approx::assert_relative_eq;
use ndarray::array;

use crate::{residual::Residual};

use super::recipes;

pub fn water()->super::super::CPA{

    let water = super::recipes::water4c();
    recipes::scpa(vec![water], vec![]).unwrap()
}

pub fn water_cs()->super::super::CPA{

    let water = super::recipes::water4c();
    recipes::cpa_cs(vec![water], vec![]).unwrap()
}

const TOL:f64 = 1e-12;

const T:f64 = 298.15;
const RHO:f64 = 58_000.;

#[test]
fn test_scpa_helmholtz() {
    
    let model = water();
    let x = array![1.0];
    let val = model.helmholtz(T, RHO, &x);
    assert_relative_eq!(val, -9.705428401736341 , epsilon = TOL)
}
#[test]
fn test_scpa_df_dt() {
    
    let model = water();
    
    let x = array![1.0];
    let val = model.df_dt(T, RHO, &x);
    assert_relative_eq!(val, 0.05592847092141594, epsilon = TOL)
}
#[test]
fn test_scpa_df_dn() {
    
    let model = water();
    let x = array![1.0];
    let val = model.df_dn(T, RHO, &x);
    assert_relative_eq!(val[0], -9.76898881609112 , epsilon = TOL)
}

#[test]
fn test_scpa_compressibility() {
    
    let model = water();
    let x = array![1.0];
    let val = - model.df_dv(T, RHO, &x) / RHO;
    
    
    assert_relative_eq!(val, -0.06356041435477426 , epsilon = TOL)
}

#[test]
fn test_cpa_carnahan_starling_helmholtz() {
    
    let model = water_cs();
    let x = array![1.0];
    let val = model.helmholtz(T, RHO, &x);
    assert_relative_eq!(val, -9.866702450522459 , epsilon = TOL)
}
#[test]
fn test_cpa_carnahan_starling_df_dt() {
    
    let model = water_cs();
    
    let x = array![1.0];
    let val = model.df_dt(T, RHO, &x);
    assert_relative_eq!(val,  0.056068888300127705, epsilon = TOL)
}
#[test]
fn test_cpa_carnahan_starling_df_dn() {
    
    let model = water_cs();
    let x = array![1.0];
    let val = model.df_dn(T, RHO, &x);
    assert_relative_eq!(val[0], -9.96984220269814 , epsilon = TOL)
}

#[test]
fn test_scpa_carnahan_starling_compressibility() {
    
    let model = water_cs();
    let x = array![1.0];
    let val = - model.df_dv(T, RHO, &x) / RHO;
    
    
    assert_relative_eq!(val, -0.10313975217568408 , epsilon = TOL)
}