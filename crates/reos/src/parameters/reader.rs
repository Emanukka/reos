

use super::records::*;
use serde::de::{DeserializeOwned};

use std::{collections::HashMap, env, fs::File, io::{ Read, Write}, iter::zip, path::Path};
use thiserror::Error;
use std::error::Error;

#[derive(Error,Debug,PartialEq)]
pub enum RecordError {
    
    #[error("component {0} not included in json")]
    NotIncluded(String),
    #[error("provide at least 1 component name")]
    NoComponents,
    #[error("components = ['{0}'] not found at {1}")]
    NotFound(String, String),
    #[error("number of binary paths ({0}) must be equal to: 0, 1, or n_sets = {1}")]
    BinaryJsonsUnmatchedSets(usize, usize),
    #[error("number of pure paths {0} must be equal to number of sets ({1})")]
    PureJsonsUnmatchedSets(usize, usize)
}

type PureRecords<M> = Vec<PureRecord<M>>;
type BinaryRecords<B> = Vec<BinaryRecord<B>>;


fn get_file<P: AsRef<Path>>(path:P)-> Result<File, Box< dyn Error>> {

    let mut current = env::current_dir()?;
    current.push(path);
    

    match File::open(&current) {
        Ok(f) => Ok(f),
        Err(e) => Err({
        println!("dir: {:?}", current.into_os_string());
            Box::new(e)
        }),
        
    }


}

fn data_from_file(mut file: File)-> String {

    let mut buf = String::new();
    let _ = file.read_to_string(&mut buf).unwrap();
    buf

}

fn get_pure_records<M: Clone, A: AsRef<str>>(
    names: &[A], 
    map: PureRecords<M>,
    path: &A)->Result<PureRecords<M>, RecordError>{
    
    let n = names.len();

    if n == 0 { return  Err(RecordError::NoComponents) }

    let hash: HashMap<&str, usize> = map.iter().enumerate().map(|(i, x)| (x.name.as_str(), i)).collect();
    let mut records = Vec::with_capacity(n);
    let mut not_found = Vec::new();

    for name in names {
        
        let name = name.as_ref();
        if let Some(&i) = hash.get(name){
            
            records.push(map[i].clone());

        } else {

            not_found.push(name);

        } 
    }

    if not_found.len() > 0 {

        // let mut s = String::with_capacity(cap);
        let s = not_found.join(", ");
        return Err(RecordError::NotFound(s, path.as_ref().to_string()))
    }

    Ok(records)

}

fn get_binary_records<B:Clone, A: AsRef<str>>(names: &[A], map: BinaryRecords<B>) -> BinaryRecords<B>{

    let n = names.len();
    let hash: HashMap<(&str, &str), usize> = map.iter().enumerate().map(|(i, x)| (x.get_id(), i)).collect();

    let mut records = Vec::with_capacity(n);

    for i in 0..n-1 {
        for j in i+1..n {
            
            let name_i = names[i].as_ref();
            let name_j = names[j].as_ref();
            let key1 = (name_i, name_j);
            let key2 = (name_j, name_i);

            if let Some(&i) = hash.get(&(key1)) {
                records.push(map[i].clone());

            } else if let Some(&i) = hash.get(&(key2)) {
                records.push(map[i].clone());
            }


        }

    }

    records
}

pub fn p_from_file<M: DeserializeOwned + Clone, A: AsRef<str>>(names: &[A], path: &A) -> Result<Vec<PureRecord<M>>, Box<dyn Error>> {
    
    let file = get_file(path.as_ref())?;
    let s = data_from_file(file);

    let map:PureRecords<M> = serde_json::from_str(&s)?;

    let records = get_pure_records(names, map, path)?;

    Ok(records)

}

pub fn p_from_files<M: DeserializeOwned + Clone, A: AsRef<str>>(sets: &[Vec<A>], ppaths: &[A], ) -> Result<Vec<PureRecord<M>>, Vec<Box<dyn Error>>> {

    let mut errors = vec![];
    let n_comp = sets.iter().fold(0, |acc, x| acc + x.len());
    let n_set = sets.len();
    let n_path = ppaths.len();
    
    if n_set != n_path {return Err(vec![RecordError::PureJsonsUnmatchedSets(n_path, n_set).into()])}

    let mut pure_records: Vec<PureRecord<M>> = Vec::with_capacity(n_comp);
    
    sets.iter().zip(ppaths.iter())
        .map(|(names, ppath)| {
            p_from_file::<M, A>(names, ppath)
        })
        .filter_map(|res| res.map_err(|e| errors.push(e)).ok())
        .for_each(|mut v| pure_records.append(&mut v));

    if errors.len() > 0 {

        // println!("{errors:?}");
        // let s = errors.to
        Err(errors)

    } else {
        
        Ok(pure_records)

    }
}



pub fn b_from_file<B: DeserializeOwned + Clone, A: AsRef<str>>(names: &[A], path: &A) -> Result<BinaryRecords<B>, Box<dyn Error>> {
    

    let file = get_file(path.as_ref())?;
    let s = data_from_file(file);
    let map:BinaryRecords<B> = serde_json::from_str(&s)?;

    let records = get_binary_records(names, map);
    Ok(records)

}
pub fn b_from_files<B: DeserializeOwned + Clone, A: AsRef<str> + Clone>(sets: &[Vec<A>], bpaths: &[A], ) -> Result<Vec<BinaryRecord<B>>, Vec<Box<dyn Error>>> {

    let n_set = sets.len();
    let n_path = bpaths.len();

    if n_path == 0 { Ok(vec![]) } 
    
    else if n_path == 1  {

        let res = b_from_file::<B, A>(&sets.concat(), &bpaths[0]);

        match res {
            Ok(records) => Ok(records),
            Err(e) => Err(vec![e])
        }
    
    } else if n_path == n_set {

        let mut records = Vec::with_capacity(n_path);
        let mut errors = vec![];

        zip(sets, bpaths)
            .map(|(names,path)| b_from_file::<B, A>(names, path) )
            .filter_map(|res|res.map_err(|e| errors.push(e)).ok())
            .for_each(|mut v| records.append(&mut v));
            
        if errors.len() > 0 { Err(errors) }
        else { Ok(records) }

    } else {
        Err(vec![RecordError::BinaryJsonsUnmatchedSets(n_path, n_set).into()])
    }

}

#[cfg(test)]
mod tests{

    use std::collections::HashMap;

    use serde::{Deserialize, Serialize};
    

    use super::*;   
    // use super::super::Parameters;
    use crate::{models::{associative::sites::CombiningRule, cpa::parameters::readyto::*, cubic::models::CubicModels}, parameters::Parameters};

    #[derive(Serialize,Deserialize, Clone)]
    pub struct PureModelTest {
        pub x: f64, 
        pub y: f64,
        pub z: f64 
    }

    #[derive(Serialize,Deserialize, Clone)]
    pub struct BinaryModelTest {
        pub kij: f64,
        pub lij: f64,
    }

    pub struct ParametersTest {
        pub pure: Vec<PureModelTest>,
        pub bin: HashMap<(usize,usize),BinaryModelTest>,
        pub prop: Properties,
    }
    impl std::fmt::Display for ParametersTest{
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            unimplemented!()
        }
    }

    impl std::fmt::Display for PureModelTest {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "Test(x={}, y={}, z={})", self.x, self.y, self.z)
        }
    }
    impl std::fmt::Display for BinaryModelTest {
        fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            write!(f, "Test(kij={}, lij={})", self.kij, self.lij)
        }
    }
    impl Parameters<PureModelTest, BinaryModelTest, ()> for ParametersTest {

        fn from_raw(pure:Vec<PureModelTest>, binary: Vec<BinaryParameter<BinaryModelTest>>, properties: Option<Properties>, _: ()) -> Self {

            // let pure = pure.iter().map(|x| x)
            Self {
                pure,
                bin: binary.into_iter().map(|x| ((x.id1, x.id2), x.model_record)).collect(),
                prop: properties.unwrap_or_default()
            }
        }
    }
    #[test]
    fn file_not_found(){

        let result = p_from_file::<Pure, &str>(&vec!["none"],&"A").unwrap_err();
        let err: &std::io::Error = result.downcast_ref().unwrap();
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);

    }

    #[test]
    fn comp_not_found(){

        let err = p_from_file::<Pure, &str>(&vec!["methane","propanol","co2"],&"src/parameters/tests/cpa/water_co2.json").unwrap_err();
        let err: &RecordError = err.downcast_ref().unwrap();
        assert_eq!(err, &RecordError::NotFound("methane, propanol".into(),"src/parameters/tests/cpa/water_co2.json".into()));
        // dbg!(err);
        // println!("{err}");

    }



    #[test]
    fn pure_jsons_unmatch_sets(){

        let names1 = vec!["water","co2"];
        let names2 = vec!["methane"];
        let sets = [names1, names2];
        let ppath1 = "src/parameters/tests/cpa/water_co2.json";

        let err = p_from_files::<Pure, &str>(&sets,&[ppath1]).unwrap_err();

        assert!(err.len() == 1);

        let err: &RecordError = err[0].downcast_ref().unwrap();
        
        assert_eq!(err, &RecordError::PureJsonsUnmatchedSets(1, 2));

        // dbg!(err);
        // println!("{err}");

    }

    #[test]
    fn binary_jsons_unmatch_sets(){

        let names1 = vec!["water","co2"];
        let names2 = vec!["methane"];
        let sets = &[names1, names2];

        let ppath1 = "src/parameters/tests/cpa/water_co2.json";

        let err = b_from_files::<Binary, &str>(sets,&[ppath1,ppath1,ppath1]).unwrap_err();

        assert!(err.len() == 1);
        let err: &RecordError = err[0].downcast_ref().unwrap();
        
        assert_eq!(err, &RecordError::BinaryJsonsUnmatchedSets(3, 2));

        // dbg!(err);
        // println!("{err}");

    }


    #[test]
    fn from_multiple_jsons_comps_not_found(){

        let names1 = vec!["water","co2"];
        let names2 = vec!["methane","propane","ethane","acetic acid"];
        // let names = vec![vec![names1], vec![names2]];
        
        let names = &[names1, names2];

        let ppath1 = "src/parameters/tests/cpa/water_co2.json";
        let ppath2 = "src/parameters/tests/cpa/acetic_acid.json";
        

        let res = CPAParameters
        ::from_multiple_jsons(names, &[ppath1,ppath2], None, CubicModels::default());

        // let res = CPAParameters::from_multiple_jsons(&names, &[ppath1,ppath2], &[], CubicModels::default());

        if let Err(e) = res {
            
            let s = format!("{:?}",e );
            // dbg!(e);
            println!("{s}");

        } else { panic!("should err")}

    }

    #[test]
    fn from_multiple_jsons() {

        let names1 = vec!["water".into()];
        let names2 = vec!["methane".into()];
        let names3 = vec!["propanol".into()];

        let sets = &[names1, names2, names3];

        let ppath1 = "src/parameters/tests/foo/water.json";
        let ppath2 = "src/parameters/tests/foo/methane.json";
        let ppath3 = "src/parameters/tests/foo/propanol.json";
        let ppaths = [ppath1, ppath2, ppath3];

        let bpath1 = "src/parameters/tests/foo/bin1.json";
        let bpath2 = "src/parameters/tests/foo/bin2.json";
        let bpath3 = "src/parameters/tests/foo/bin3.json";

        let bpaths = [bpath1, bpath2, bpath3];

        let bpaths = Some(bpaths.as_slice());

        let p = ParametersTest::from_multiple_jsons(sets, 
            &ppaths, 
            bpaths, 
            ()).unwrap();
        
        let water = "Test(x=1, y=2, z=3)".to_string();
        let methane = "Test(x=4, y=5, z=6)".to_string();
        let propanol = "Test(x=7, y=8, z=9)".to_string();
        
        let target = vec![water, methane, propanol];

        let pure = &p.pure;
        let bin = &p.bin;
        for i in 0..pure.len(){
            let s = format!("{}", pure[i]);
            assert_eq!(s, target[i]);
        }

        let water_methane = "Test(kij=1, lij=2)".to_string();
        let water_propanol = "Test(kij=3, lij=4)".to_string();
        let methane_propanol = "Test(kij=5, lij=6)".to_string();
        let target = vec![water_methane, water_propanol, methane_propanol];

        bin.iter().enumerate().for_each(|(i,((_id1,_id2), model))| {
            let s = format!("{}", model);
            assert_eq!(s, target[i]);
        });

    }
    #[test]
    fn cpa_from_json(){

        let ppath = "src/parameters/tests/cpa/water_co2.json";
        let bpath = "src/parameters/tests/cpa/bin.json";

        let p = CPAParameters::from_json(&["water","co2"], &ppath, Some(&bpath), CubicModels::default()).unwrap();

        let aij = p.cubic.aij.as_slice().unwrap();
        let bij = p.cubic.bij.as_slice().unwrap();
        let interactions= &p.assoc.interactions;

        assert_eq!(interactions.len(), 2);
        assert_eq!(interactions[1].epsilon, 166.55e2 /2.);
        assert_eq!(interactions[1].kappa,   0.1836);

        assert_eq!(aij, &[0.0, -0.15508, -0.15508, 0.0]);
        assert_eq!(bij, &[0.0, 0.000877, 0.000877, 0.0]);

        let s = p.to_string();
        println!("{}", s);
    }
    

    #[test]
    fn cpa_from_multiple_jsons(){

        let names1 = vec!["water"];
        let names2 = vec!["co2"];

        let sets = [names1, names2];

        let ppath1 = "src/parameters/tests/cpa/water_co2.json";
        let ppath2 = "src/parameters/tests/cpa/co2.json";
        let ppaths = [ppath1, ppath2];
        let bpaths = Some(["src/parameters/tests/cpa/bin.json"].as_slice());    

        let p = CPAParameters::from_multiple_jsons(&sets, 
            &ppaths, 
            bpaths, 
            CubicModels::default()).unwrap();
        
        let aij = p.cubic.aij.as_slice().unwrap();
        let bij = p.cubic.bij.as_slice().unwrap();
        let interactions= &p.assoc.interactions;

        assert_eq!(interactions.len(), 2);
        assert_eq!(interactions[1].epsilon, 166.55e2 /2.);
        assert_eq!(interactions[1].kappa,   0.1836);

        assert_eq!(aij, &[0.0, -0.15508, -0.15508, 0.0]);
        assert_eq!(bij, &[0.0, 0.000877, 0.000877, 0.0]);

        let s = p.to_string();
        println!("{}", s);
    }

    #[test]
    fn cpa_from_multiple_jsons2(){

        let names1 = vec!["water"];
        let names2 = vec!["acetic acid"];

        let sets = [names1, names2];

        let ppath1 = "src/parameters/tests/cpa/water_co2.json";
        let ppath2 = "src/parameters/tests/cpa/acetic_acid.json";
        let ppaths = [ppath1, ppath2];
        let bpaths = Some(["src/parameters/tests/cpa/bin.json"].as_slice());    

        let p = CPAParameters::from_multiple_jsons(&sets, 
            &ppaths, 
            bpaths, 
            CubicModels::default()).unwrap();
        
        // let aij = p.cubic.aij.as_slice().unwrap();
        // let bij = p.cubic.bij.as_slice().unwrap();

        let interactions= &p.assoc.interactions;

        // assert_eq!(interactions[1].epsilon, 166.55e2 /2.);
        // assert_eq!(interactions[1].kappa,   0.1836);

        // assert_eq!(aij, &[0.0, -0.15508, -0.15508, 0.0]);
        // assert_eq!(bij, &[0.0, 0.000877, 0.000877, 0.0]);
        let s = p.to_string();
        println!("{}", s);

        let cr = &[CombiningRule::CR1, CombiningRule::ECR, CombiningRule::ECR,CombiningRule::CR1];
        for (i, interaction) in interactions.iter().enumerate(){
            assert_eq!(interaction.combining_rule, cr[i]);
        }

        // assert_eq!(interactions[1].combining_rule, CombiningRule::ECR);


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
