use crate::{models::cubic::{alpha::AlphaRecord, parameters::CubicPureRecord}, residual::Residual};

use super::recipes;

use approx::assert_relative_eq;
use ndarray::array;

#[test]
fn test_der_record() {
    
    // let cub = recipes::hexane_pr78_twu91(); 
    let s = r#"{
    "molar_weight": 0.0,
    "name": "1,1-dimethoxyethane",
    "tc": 507.8,
    "c": -1.5928999999999998e-6,
    "pc": 3772999.9999999995,
    "alpha":
        {
            "type":"twu91",
            "l": 0.1447,
            "m": 0.8299,
            "n": 2.9049
        }
    }"#;
    // let s = r#"{
    // "molar_weight": 0.0,
    // "name": "1,1-dimethoxyethane",
    // "tc": 507.8,
    // "c": -1.5928999999999998e-6,
    // "pc": 3772999.9999999995,
    // "l": 0.1447,
    // "m": 0.8299,
    // "n": 2.9049
    // }"#;
    let r:CubicPureRecord = serde_json::from_str(s).unwrap();

    println!("{r}")
    // let s = serde_json::to_string_pretty(&cub.parameters).unwrap();

    // println!("{}", s);


    // println!("{}", cub.parameters);
    
    // assert_relative_eq!(cub.helmholtz(T, RHO, &array![1.]), -0.37074849150969297, epsilon = TOL)
}
#[test]
fn test_der_record2() {
    
    // let cub = recipes::hexane_pr78_twu91(); 
    let s = r#"{
        "type": "soave_regressed",
        "c1": 1.0
        }"#;
    // let s = r#"{
    //     "type": "soave",
    //     "l": 0.1447,
    //     "m": 0.8299,
    //     "n": 2.9049
    //     }"#;
    // let s = r#"{
    // "molar_weight": 0.0,
    // "name": "1,1-dimethoxyethane",
    // "tc": 507.8,
    // "c": -1.5928999999999998e-6,
    // "pc": 3772999.9999999995,
    // "l": 0.1447,
    // "m": 0.8299,
    // "n": 2.9049
    // }"#;
    let r:AlphaRecord = serde_json::from_str(s).unwrap();

    println!("{r}")
    // let s = serde_json::to_string_pretty(&cub.parameters).unwrap();

    // println!("{}", s);


    // println!("{}", cub.parameters);
    
    // assert_relative_eq!(cub.helmholtz(T, RHO, &array![1.]), -0.37074849150969297, epsilon = TOL)
}