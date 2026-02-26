
use approx::assert_relative_eq;
use ndarray::{Array1, array};

use crate::models::associative::parameters::AssociativeParameters;
use crate::parameters::Parameters;

use super::super::Associative;

use super::recipes;

const TOL:f64 = 1e-12;

fn water_carbon_dioxide() -> Associative {

    let w = recipes::water();
    let co2 = recipes::carbon_dioxide();
    let bin = recipes::water_carbon_dioxide();

    let parameters = AssociativeParameters::new(vec![w,co2], vec![bin], ()).unwrap();
    Associative{parameters}

}
fn water_inert_carbon_dioxide() -> Associative {

    let w = recipes::water();
    let co2 = recipes::carbon_dioxide();

    let parameters = AssociativeParameters::new(vec![w,co2], vec![], ()).unwrap();
    Associative{parameters}

}
fn water_acetic_acid() -> Associative {

    let w = recipes::water();
    let acoh = recipes::acetic_acid();
    let bin = recipes::water_acetic_acid();

    let parameters = AssociativeParameters::new(vec![w,acoh], vec![bin], ()).unwrap();
    Associative{parameters}

}

// mod unbonded_sites {

//     use approx::assert_relative_eq;
//     use ndarray::array;

//     use crate::models::associative::sites::CombiningRule;

//     use super::*;
    
//     #[test]
//     fn test_4c_inert() {

//         // t=298.15, d=55_000, x=[0.5,0.5]

//         let asc = water_inert_carbon_dioxide();
//         let p = &asc.parameters;
//         let interactions = &p.interactions;
        
//         assert_eq!(interactions.len(), 1);

//         let self_water = &interactions[0];
//         // let solv_water_co2 = &interactions[1];

        
//         assert_eq!(self_water.epsilon, 166.55e2);
//         assert_eq!(self_water.kappa, 0.0692);
//         assert_eq!(self_water.combining_rule, CombiningRule::CR1);
        
//         let k = array![[ 0.              , 25.04896563685188,  0.              ],
//        [25.04896563685188,  0.              ,  0.              ],
//        [ 0.              ,  0.              ,  0.              ]];
        
//         let x = array![0.5, 0.5];
//         let m = &asc.sites_mole_frac(&x);
//         let xtan = asc.x_tan(m, &k).unwrap();

//         assert_relative_eq!(xtan, array![0.09503654234105, 0.09503654234105, 1.              ], epsilon = TOL);

//     }

//     #[test]
//     fn test_4c_solvation() {

//         // t=298.15, d=55_000, x=[0.5,0.5]

//         let asc = water_carbon_dioxide();
//         let p = &asc.parameters;
//         let interactions = &p.interactions;
        
//         assert_eq!(interactions.len(), 2);

//         let self_water = &interactions[0];
//         let solv_water_co2 = &interactions[1];

        
//         assert_eq!(self_water.epsilon, 166.55e2);
//         assert_eq!(self_water.kappa, 0.0692);
//         assert_eq!(solv_water_co2.epsilon, 0.5 * 166.55e2);
//         assert_eq!(solv_water_co2.kappa, 0.1836);
//         assert_eq!(solv_water_co2.combining_rule, CombiningRule::CR1);

//         let k = array![[0., 25.048965636851882, 3.210256575812608],
//                                 [25.048965636851882, 0., 0.],
//                                 [3.210256575812608, 0., 0.]];

//         let x = array![0.5, 0.5];
//         let m = &asc.sites_mole_frac(&x);
//         let xtan = asc.x_tan(m, &k).unwrap();

//         assert_relative_eq!(xtan, array![0.03874121673562, 0.20484626829382, 0.66778989678062], epsilon = TOL);

//     }   

//     #[test]
//     fn test_4c_1a_elliot_cr() {
        
//         // t=298.15, d=1_000, x=[0.5,0.5]

//         let asc = water_acetic_acid();
//         let p = &asc.parameters;
//         let interactions = &p.interactions;
        
//         assert_eq!(interactions.len(), 4);

//         let self_water = &interactions[0];
//         let self_acoh = &interactions[3];
//         let cross_water_acoh = &interactions[1];

//         assert_eq!(self_water.epsilon, 166.55e2);
//         assert_eq!(self_water.kappa, 0.0692);
//         assert_eq!(self_acoh.epsilon, 40323.);
//         assert_eq!(self_acoh.kappa,0.0045);

//         assert_eq!(cross_water_acoh.epsilon, 28489.);
//         assert_eq!(cross_water_acoh.kappa, 0.017646529403823292);
//         assert_eq!(cross_water_acoh.combining_rule, CombiningRule::ECR);

//         let k = array![
//             [  0.              ,   0.47210091350445,   8.53435728830102],
//             [  0.47210091350445,   0.              ,   8.53435728830102],
//             [  8.53435728830102,   8.53435728830102, 154.27899468296786]];

//         let x = array![0.75, 0.25];
//         let m = &asc.sites_mole_frac(&x);
//         let xtan = asc.x_tan(m, &k).unwrap();

//         assert_relative_eq!(xtan, array![0.54759463290882, 0.54759463290882, 0.0120204045421 ], epsilon = TOL);

//     }    

    
// }

