
use super::records::*;
use serde::de::{DeserializeOwned};

use std::{env, fs::File, io:: Read, path::Path};
use thiserror::Error;
use std::error::Error;

#[derive(Error,Debug,PartialEq)]
pub enum RecordError {
    
    #[error("component {0} not included in json")]
    NotIncluded(String),
    #[error("provide at least 1 component name")]
    NoComponents,


}

type PureRecords<M> = Vec<PureRecord<M>>;
type BinaryRecords<B> = Vec<BinaryRecord<B>>;


fn get_file<P: AsRef<Path>>(path:P)-> Result<File, Box< dyn Error>> {

    let mut current = env::current_dir()?;
    current.push(path);
    
    let val = File::open(current)?;

    Ok(val)

}

fn data_from_file(mut file: File)-> String {

    let mut buf = String::new();
    let _ = file.read_to_string(&mut buf).unwrap();
    buf

}

fn get_pure_records<M>(
    names: &[&str], 
    map: PureRecords<M>)->Result<PureRecords<M>, RecordError>{
    
    let n = names.len();

    if n == 0 { return  Err(RecordError::NoComponents) }

    let mut records = Vec::with_capacity(n);

    
    let mut map = map.into_iter();
    for &name in names {
        
        if let Some(r) = map.find(|r| r.name == name) {
            
            records.push(r);

        } else {
            
            return Err(RecordError::NotIncluded(name.to_string()))

        } 
    }

    Ok(records)

}

fn get_binary_records<M>(names: &[&str], map: BinaryRecords<M>) -> BinaryRecords<M>{

    let n = names.len();
    let mut records = Vec::with_capacity(n);
    let mut iter = map.into_iter();

    for i in 0..n-1 {
        for j in i+1..n {
            
            let name_i = names[i];
            let name_j = names[j];

            if let Some(r) = iter.find(|x| { 
                let id = x.get_id();
                id == (name_i, name_j) || id == (name_j, name_i)
            }) {
                records.push(r);
            }

        }

    }

    records
}

pub fn p_from_file<M: DeserializeOwned>(names: &[&str], path: &str) -> Result<Vec<PureRecord<M>>, Box<dyn Error>> {
    
    let file = get_file(path)?;
    let s = data_from_file(file);
    let map:PureRecords<M> = serde_json::from_str(&s)?;

    let records = get_pure_records(names, map)?;

    Ok(records)

}



pub fn b_from_file<B: DeserializeOwned>(names: &[&str], path: &str) -> Result<BinaryRecords<B>, Box<dyn Error>> {
    

    let file = get_file(path)?;
    let s = data_from_file(file);
    let map:BinaryRecords<B> = serde_json::from_str(&s)?;

    let records = get_binary_records(names, map);
    Ok(records)

}


#[cfg(test)]
mod tests{

    use super::*;   
    use super::super::Parameters;
    use crate::models::cpa::parameters::readyto::*;
    use crate::models::cubic::SRK;

    #[test]
    fn file_not_found(){

        let result = p_from_file::<Pure>(&vec!["none".into()],"A").unwrap_err();
        let err: &std::io::Error = result.downcast_ref().unwrap();
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);

    }

    #[test]
    fn comp_not_found(){

        let err = p_from_file::<Pure>(&vec!["none"],"src/parameters/pure.json").unwrap_err();
        let err: &RecordError = err.downcast_ref().unwrap();
        assert_eq!(err, &RecordError::NotIncluded("none".into()));

    }


    #[test]
    fn cpa_from_json(){

        let names = &["water","co2"];
        let ppath = "src/parameters/pure.json";
        let bpath = "src/parameters/bin.json";

        let p = CPAParameters::<SRK>::from_json(names, ppath, bpath.into()).unwrap();

        let aij = p.cubic.aij.as_slice().unwrap();
        let bij = p.cubic.bij.as_slice().unwrap();
        let interactions= p.assoc.interactions;

        assert_eq!(interactions.len(), 2);
        assert_eq!(interactions[1].epsilon, 166.55e2 /2.);
        assert_eq!(interactions[1].kappa,   0.1836);

        assert_eq!(aij, &[0.0, -0.15508, -0.15508, 0.0]);
        assert_eq!(bij, &[0.0, 0.000877, 0.000877, 0.0]);


    }
}

 // #[test]
    // fn test_from_json(){
    //     let s = r#"
    //         [

    //             {   
    //                 "name": "water",
    //                 "a0":   0.12277,
    //                 "b":    0.0145e-3, 
    //                 "c1":   0.6736, 
    //                 "tc":   647.14,
    //                 "na":   2,
    //                 "nb":   2,
    //                 "epsilon": 166.55e2,
    //                 "kappa": 0.0692,
    //                 "molar_weight": 18.01528
    //             },

    //             {   
    //                 "name": "co2",
    //                 "a0":   0.35079, 
    //                 "b":    0.0272e-3, 
    //                 "c1":   0.7602, 
    //                 "tc":   304.12,
    //                 "nb":   1
    //             }
    //         ]

    //     "#;
        
    //     let records: Vec<Pure> = serde_json::from_str(s).unwrap();
        
    //     let s = r#"[
    //         {
    //             "id1": "water",
    //             "id2": "co2",
    //             "aij": -0.15508,
    //             "bij": 0.000877,
    //             "kappa": 0.1836
    //         }

    //     ]"#;

    //     let binary: Vec<Binary> = serde_json::from_str(s).unwrap();
    //     // let bin: Vec<Binary> = s
    //     let p = CPAParameters::new(records, binary);

    //     let cpa = SCPA::from_parameters(p);
    //     let c = serde_json::to_string_pretty(&cpa.cubic.parameters).unwrap();
    //     let a = serde_json::to_string_pretty(&cpa.assoc.assoc.parameters).unwrap();

    //     // println!("{c}");
    //     // println!("{a}");

    //     let asc = cpa.assoc.assoc.parameters;

    //     assert_eq!(asc.interactions[1].epsilon, 166.55e2 /2.);
    //     assert_eq!(asc.interactions[1].kappa,   0.1836);

    //     // cpa.assoc.assoc.parameters.interactions[]
    // }
