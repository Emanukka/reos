use approx::assert_relative_eq;
use ndarray::array;


use crate::models::cubic::Cubic;
use crate::parameters::{PureRecord, BinaryRecord,Parameters};
use crate::residual::Residual;
use super::super::parameters::{CubicPureRecord, CubicBinaryRecord, CubicParameters};
use super::super::options::*;


use super::recipes::water_co2_bip;


#[test]
fn water_co2_pr78_from_json() {

    let s = r#"{
    "tc": 647.1,
    "pc": 220.55e5,
    "w" : 0.345
}"#;
    // println!("{}",s);
    
    let water: CubicPureRecord = serde_json::from_str(s).unwrap();
    
    println!("{}", &water);

    let s = r#"{
    "tc": 304.2,
    "pc": 73.83e5,
    "w" : 0.224
}"#;
    // println!("{}",s);
    
    let co2: CubicPureRecord = serde_json::from_str(s).unwrap();
    
    println!("{}", &co2);

    let s = r#"{
    "kij": 0.1,
    "lij": 0.05
}"#;
    let mb: CubicBinaryRecord = serde_json::from_str(s).unwrap();
    println!("{}",&mb);

    let pr1 = PureRecord::new(0., "water", water);
    let pr2 = PureRecord::new(0., "carbon dioxide", co2);

    let br = BinaryRecord::new(mb, "water", "carbon dioxide");
    
    // let options = CubicOptions::new(PR78.into(), Alpha::soave(), MixingRule::default());
    let options = CubicOptions::classic_soave(CubicModelOption::PR78);

    let p = CubicParameters::new(vec![pr1, pr2], vec![br], options).unwrap();

    let cub = Cubic::from_parameters(p);
    let reff = water_co2_bip(); 
    let t = 298.15;
    let d= 38.082099077791675;
    let x = &array![0.6, 0.4];


    assert_relative_eq!(cub.r_entropy(t, d, x), reff.r_entropy(t, d, x))
    // }
}

#[test]
fn options(){

    let options = CubicOptions::classic_soave(CubicModelOption::SRK);

    let s = r#"
    {
        "cubic_model": "srk",
        "alpha_model": "soave"
    }
    "#;

    assert_eq!(options, serde_json::from_str(s).unwrap());

    let options = CubicOptions::classic_soave(CubicModelOption::PR78);

    let s = r#"
    {
        "cubic_model": "pr78",
        "alpha_model": "soave"
    }
    "#;

    assert_eq!(options, serde_json::from_str(s).unwrap());

    let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);

    let s = r#"
    {
        "cubic_model": "pr78",
        "alpha_model": "twu91"
    }
    "#;

    assert_eq!(options, serde_json::from_str(s).unwrap());

    // assert_eq!(options, classic_soave_srk); 
}
// #[test]
// fn options_from_json() {

// //     let s = r#"{
// //     "tc": 647.1,
// //     "pc": 220.55e5,
// //     "w" : 0.345
// // }"#;
// //     // println!("{}",s);
    
// //     let water: CubicPureRecord = serde_json::from_str(s).unwrap();
    
// //     println!("{}", &water);

// //     let s = r#"{
// //     "tc": 304.2,
// //     "pc": 73.83e5,
// //     "w" : 0.224
// // }"#;
// //     // println!("{}",s);
    
// //     let co2: CubicPureRecord = serde_json::from_str(s).unwrap();
    
// //     println!("{}", &co2);

// //     let s = r#"{
// //     "kij": 0.1,
// //     "lij": 0.05
// // }"#;
// //     let mb: CubicBinaryRecord = serde_json::from_str(s).unwrap();
// //     println!("{}",&mb);

// //     let pr1 = PureRecord::new(0., "water", water);
// //     let pr2 = PureRecord::new(0., "carbon dioxide", co2);

// //     let br = BinaryRecord::new(mb, "water", "carbon dioxide");
    
//     // let options = CubicOptions::new(PR78.into(), Alpha::soave(), MixingRule::default());
//     let s = r#"
// "mixingrule":{
//     "quadratic":"quadratic"
// }"#;
//     let mix: MixingRule = serde_json::from_str(s).unwrap();
//     // let mix = MixingRule::default();

//     let s = serde_json::to_string(&mix).unwrap();
//     // println!("{}",s);

//     // let mix: MixingRule = serde_json::from_str(&s).unwrap();

//     // let options = CubicOptions::classic(PR78.into());

//     // let p = CubicParameters::new(vec![pr1, pr2], vec![br], options);

//     // let cub = Cubic::from_parameters(p);
//     // let reff = water_co2_bip(); 
//     // let t = 298.15;
//     // let d= 38.082099077791675;
//     // let x = &array![0.6, 0.4];


//     // assert_relative_eq!(cub.r_entropy(t, d, x), reff.r_entropy(t, d, x))
//     // }
// }