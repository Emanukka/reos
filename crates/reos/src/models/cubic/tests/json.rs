use approx::assert_relative_eq;
use ndarray::array;

use crate::{models::cubic::{Cubic, mixing_rule::MixingRule, models::PR78, options::CubicOptions, parameters::{CubicBinaryRecord, CubicParameters, CubicPureRecord}}, parameters::{BinaryRecord, Parameters, PureRecord}, residual::Residual};

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
    let options = CubicOptions::classic_soave(PR78.into());

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
fn bla(){
        use ::serde::{ Deserialize, Serialize };

        #[derive(Deserialize, Serialize)]
        struct X;

        #[derive(Deserialize, Serialize)]
        struct Y;


    
        #[derive(Serialize,Deserialize)]
        // #[serde(untagged)]
        enum Z{
            X(X),
            Y(Y)
        }

        #[derive(Deserialize, Serialize)]
        struct W{
            z:Z,
            i:usize
        }
        // let ty = List { elements: vec![42] };
        let z = Z::X(X);
        let serialized = serde_json::to_string(&z).unwrap();

        dbg!(serialized);

        let w = W{z,i:1};
        let serialized = serde_json::to_string(&w).unwrap();

        println!("{}",serialized);
        // assert_eq!(serialized, r#"[42]"#);

        // #[derive(Deserialize, Serialize)]
        // struct ListObject<T> {
        //     elements: Vec<T>,
        // }

        // let ty = ListObject { elements: vec![42] };
        // let serialized = serde_json::to_string(&ty).unwrap();

        // assert_eq!(serialized, r#"{"elements":[42]}"#);
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