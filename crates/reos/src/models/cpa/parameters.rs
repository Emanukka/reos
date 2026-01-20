use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::models::associative::parameters::{AssociationBinaryRecord, AssociationPureRecord, AssociativeParameters};

use crate::models::cubic::models::CubicModels;
use crate::models::cubic::parameters::{CubicBinaryRecord, CubicParameters, CubicPureRecord};
use crate::parameters::{Parameters};
use crate::parameters::records::{BinaryParameter, BinaryRecord, PureRecord};


#[derive(Serialize,Deserialize, Debug, Clone)]
pub struct CPAPureRecord{
    #[serde(flatten)]
    pub c:CubicPureRecord,
    #[serde(flatten)]
    pub a:AssociationPureRecord

}

impl CPAPureRecord{

    pub fn new(c:CubicPureRecord,a:AssociationPureRecord)->Self{
        Self { c, a }
    }

}

impl std::fmt::Display for CPAPureRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "CPAPureRecord(c={}, a={})", self.c, self.a)
        
    }
}

#[derive(Clone)]
pub struct CPAParameters{
    pub cubic: CubicParameters, 
    pub assoc: AssociativeParameters,
}

#[derive(Serialize,Deserialize,Clone,Debug)]

pub struct CPABinaryRecord{
    #[serde(flatten)]
    pub a: Option<AssociationBinaryRecord>,
    #[serde(flatten)]
    pub c: Option<CubicBinaryRecord>

}

impl std::fmt::Display for CPABinaryRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let mut data = String::from_str("CPABinaryRecord(\n").unwrap();


        if let Some(c) = &self.c {
            data.push_str(concat!("cubic="));
            data.push_str(c.to_string().as_str());
            // let _ = write!(f, "CPABinaryRecord(cubic={:#?}, assoc={:#?})", c);
        } else {
            data.push_str("cubic=None");
        }
        data.push_str("\n");

        if let Some(a) = &self.a {
            data.push_str(concat!(" assoc="));
            data.push_str(a.to_string().as_str());            
        }
        
        else {
            data.push_str(" assoc=None");
        }

        data.push_str("\n)");

        write!(f, "{}", data)

    }
}

impl CPABinaryRecord{

    pub fn new(c:Option<CubicBinaryRecord>,a:Option<AssociationBinaryRecord>)->Self{
        Self { c, a }
    }

    pub fn full(c:CubicBinaryRecord, a:AssociationBinaryRecord)->Self{
        
        Self::new(Some(c), Some(a))
    }

    pub fn only_c(c:CubicBinaryRecord)->Self{
        
        Self::new(Some(c), None)
    }

    pub fn only_a(a:AssociationBinaryRecord)->Self{
        
        Self::new(None, Some(a))
    }

}


type Pure = CPAPureRecord;
type Binary = CPABinaryRecord;



impl Parameters<Pure, Binary, CubicModels> for CPAParameters {

    fn from_raw(pure:Vec<Pure>, binary: Vec<BinaryParameter<Binary>>, properties: Option<crate::parameters::Properties>, opt: CubicModels) -> Self {
    
        let n = pure.len();

        // CPAbinMap -> CPAmodelRecord -> C,A -> CParamets,AParamets
        let mut c_pure= Vec::with_capacity(n);
        let mut a_pure = Vec::with_capacity(n);
        let mut c_binary = vec![];
        let mut a_binary = vec![];
        
        for b in binary{

            let id1 = b.id1;
            let id2 = b.id2;

            if let Some(c) = b.model_record.c { 

                c_binary.push(
                    BinaryParameter::new(c, id1, id2)
                );
            }

            if let Some(a) = b.model_record.a {
                a_binary.push(
                    BinaryParameter::new(a, id1, id2)
                );
            }
                
        }
        for record in pure{
            c_pure.push(record.c);
            a_pure.push(record.a);
        }

        let cubic = CubicParameters::from_raw(c_pure, c_binary, properties, opt);
        let assoc = AssociativeParameters::from_raw(a_pure, a_binary, None, ());

        CPAParameters{
            cubic,
            assoc,
        }

    }
}

impl std::fmt::Display for CPAParameters {
    
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "CPAParameters(\n  {},\n  {}\n)", self.cubic, self.assoc)
    }

}
// impl<C: CubicModel>  CPAParameters<C> {

//     pub fn from_records(records:Vec<CPAPureRecord>,binary:Vec<CPABinaryRecord>)->Self {
        
//         let n = records.len();
//         let mut c_records= Vec::with_capacity(n);
//         let mut a_records = Vec::with_capacity(n);
//         let mut c_binary = vec![];
//         let mut a_binary = vec![];

//         for record in records{
//             c_records.push(record.c);
//             a_records.push(record.a);
//         }

//         for b in binary{

//             if let Some(c) = b.c { c_binary.push(c); }

//             if let Some(a) = b.a { a_binary.push(a); }
                
//         }

//         let cubic_params=CubicParameters::<C>::new_(C::model(), c_records,c_binary);
//         let asc_params=AssociativeParameters::from_records(a_records,a_binary);

//         CPAParameters{
//             cubic:cubic_params,
//             assoc:asc_params,
//             }
//     }
    

// }
pub mod readyto{

    pub use super::*;
    use super::super::super::associative::sites::CombiningRule;
    
    pub type Pure = PureRecord<CPAPureRecord>; 
    pub type Binary = BinaryRecord<CPABinaryRecord>; 

    #[macro_export]
    macro_rules! arr_eq {
        ($a:expr, $b:expr, $tol:expr) => {{
            if $a.len() != $b.len() {
                println!("dim wrong!");
                false
            } else {
                $a.iter().zip($b.iter()).all(|(x, y)| {

                  if (x - y).abs() < $tol + $tol * y.abs(){
                    true
                  }
                  else {
                    println!("left = {}, right = {}",x,y);
                    false
                  }
                })
              }}
          }
    }
    pub fn water4c()->Pure{

        let c = CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);

        let a = AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0]);

        let m = CPAPureRecord::new(c, a);
            
        PureRecord::new(0.0, "water".to_string(), m)

    }

    pub fn acetic1a()->Pure{
        let c=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);

        let a = AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1]);

        let m = CPAPureRecord::new(c, a);
        PureRecord::new(0.0, "acetic_acid".to_string(), m)

    }

    pub fn water4c_acetic1a()->Binary{

        let c = CubicBinaryRecord::TemperatureIndependent { kij: -0.222 };
        
        let a = AssociationBinaryRecord {  epsilon: None, kappa: None, combining_rule: CombiningRule::ECR };
        
        let b = CPABinaryRecord::full(c, a);

        BinaryRecord::new(b, "water".into(), "acetic_acid".into())

    }

    pub fn co2()->Pure{

        let c=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

        let a=AssociationPureRecord::solvate(
            [0,1,0]);

        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "co2".to_string(), m)

    }

    pub fn water4c_co2()->Binary{

        let c = CubicBinaryRecord::TemperatureDependent { aij: -0.15508 , bij: 0.000877 };
        
        let a = AssociationBinaryRecord {epsilon: None, kappa: Some(0.1836), combining_rule: CombiningRule::default() };
        
        let b = CPABinaryRecord::full(c, a);
        
        BinaryRecord::new(b, "water".into(), "co2".into())

    }

    pub fn octane()->Pure{
        let c=CubicPureRecord::new_set1(34.8750e-1, 0.1424e-3, 0.99415, 568.7);
        let a=AssociationPureRecord::inert();


        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "octane".to_string(), m)    
    }
    pub fn acoh_octane()->Binary{


        let c = CubicBinaryRecord::TemperatureIndependent { kij: 0.064 };
        let a = AssociationBinaryRecord {epsilon: None, kappa: None, combining_rule: CombiningRule::default() };
        let b = CPABinaryRecord::full(c, a);
        
        BinaryRecord::new(b, "acetic_acid".into(), "octane".into())    

    } 
    pub fn methanol3b()->Pure{

            let c=
            CubicPureRecord::new_set1(
                4.5897e-1, 
                0.0334e-3, 
                1.0068,
                513.);

            let a=AssociationPureRecord::associative(
                160.70e2, 
                34.4e-3, 
                [2,1,0],
            );

        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "methanol".to_string(), m)    

    } 

}

// #[cfg(test)]
// mod tests{
//     use std::{collections::HashMap, error::Error, path::Path};

//     use serde_json::{from_str, to_string, to_string_pretty};

//     use crate::{models::{associative::parameters::AssociativeParameters, cpa::{SCPA, parameters::{CPABinaryRecord, CPAParameters, CPAPureRecord}}, cubic::SRK}, parameters::{Parameters, records::{hashmap_from_file, pure_record_from_file, pure_records_from_file}}, residual::Residual};

//     use super::readyto::*;


    

    
//     #[test]
//     fn test_from_json_file(){

//         // let record: Pure = pure_record_from_file("src/models/cpa/water.json").unwrap();
//         let records:Vec<Pure>= pure_records_from_file(vec!["a".into()],"src/models/cpa/pure.json").unwrap();
//         // let records:HashMap<String,Pure> = hashmap_from_file(vec![],"src/models/cpa/pure.json").unwrap();

//         println!("{:#?}",&records);
//         // let records: Vec<Pure> = serde_json::from_reader("pure.json").unwrap();
        

        
//         // let binary: Vec<Binary> = serde_json::from_str(s).unwrap();
//         // // let bin: Vec<Binary> = s
//         // let p = CPAParameters::new(records, binary);

//         // let cpa = SCPA::from_parameters(p);
//         // let c = to_string_pretty(&cpa.cubic.parameters).unwrap();
//         // let a = to_string_pretty(&cpa.assoc.assoc.parameters).unwrap();

//         // // println!("{c}");
//         // // println!("{a}");

//         // let asc = cpa.assoc.assoc.parameters;

//         // assert_eq!(asc.interactions[1].epsilon, 166.55e2 /2.);
//         // assert_eq!(asc.interactions[1].kappa,   0.1836);

//         // cpa.assoc.assoc.parameters.interactions[]
//     }
//     #[test]
//     fn cpa_records_json(){

//         let data1 = r#"
//         {   
//             "name": "water",
//             "a0":   0.12277,
//             "b":    0.0145e-3, 
//             "kappa":0.6736, 
//             "tc":   647.14,
//             "na":   2,
//             "nb":   2,
//             "epsilon": 166.55e2,
//             "beta": 0.0692,
//             "molar_weight": 18.01528
//         }
//         "#;
        
//         let data2 = r#"
//         {   
//             "name": "co2",
//             "a0":   0.35079, 
//             "b":    0.0272e-3, 
//             "kappa":0.7602, 
//             "tc":   304.12,
//             "nb":   1
//         }
//         "#;

//         let c1:Pure = from_str(data1).unwrap();
//         let c2:Pure = from_str(data2).unwrap();

//         let c1_string = serde_json::to_string_pretty(&c1).unwrap();
//         let c2_string = serde_json::to_string_pretty(&c2).unwrap();

//         println!("{}",c1_string);

//         println!("{}",c2_string);

//         let p = CPAParameters::new(vec![c1,c2], vec![]);

//         let cpa = SCPA::from_parameters(p);

//         dbg!(cpa.molar_weight());
//         // let p = AssociativeParameters::from_records(vec![c1,c2]);

//         // let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
//         // let c2=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

//         // let a1=AssociationPureRecord::associative(
//         //     166.55e2, 
//         //     0.0692, 
//         //     [2,2,0],
//         // );
//         // let a2=AssociationPureRecord::solvate(
//         //     [0,1,0]);

//         // let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
//         // let mut parameters = CPAParameters::from_records(records);
//         // // parameters.cubic.set_kij(0, 1, -0.222);

//         // parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
//         // parameters.assoc.set_binary_from_owners(0, 1, None, Some(0.1836));

//         // let cpa = SCPAsrkCR1::from_parameters(parameters);
//         // //Create new State
//         // E::from_residual(cpa)
//     }

//     #[test]
//     fn cpa_binary_record(){

//         let data = r#"
//         {
//             "kappa": 0.1836,
//             "aij": -0.15508,
//             "bij": 0.000877
//         }
//         "#;
        

//         let b:CPABinaryRecord = from_str(data).unwrap();
//         println!("{}",serde_json::to_string_pretty(&b).unwrap())

//                 // parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
//         // parameters.assoc.set_binary_from_owners(0, 1, None, Some(0.1836));
//     }
// }