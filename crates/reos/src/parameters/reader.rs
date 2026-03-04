

use super::records::*;
use serde::de::{DeserializeOwned};

use std::{collections::HashMap, env, fs::File, io::{ Read, Write}, iter::zip, path::Path};
use thiserror::Error;
use std::error::Error;

#[derive(Clone, Error, Debug,PartialEq)]
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

fn get_binary_records<B:Clone, A: AsRef<str>>(names: &[A], map: Vec<BinaryRecord<B>>) -> Vec<BinaryRecord<B>>{

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
    use crate::parameters::Parameters;

    #[derive(Serialize,Deserialize, Clone, Debug)]
    pub struct PureModelTest {
        pub x: f64, 
        pub y: f64,
        pub z: f64 
    }

    #[derive(Serialize,Deserialize, Clone, Debug)]
    pub struct BinaryModelTest {
        pub kij: f64,
        pub lij: Option<f64>,
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
            write!(f, "Test(kij={}, lij={})", self.kij, self.lij.unwrap_or_default())
        }
    }
    impl Parameters for ParametersTest {

        type Pure = PureModelTest;
        type Binary =  BinaryModelTest;
        type Options = ();
        fn from_raw(pure:Vec<PureModelTest>, binary: crate::parameters::BinaryMap<BinaryModelTest>, properties: Option<Properties>, _: ()) -> Result<Self, Box<dyn Error>> {
            // let pure = pure.iter().map(|x| x)
            Ok(Self {
                pure,
                bin: binary,
                prop: properties.unwrap_or_default()
            })
        }
    }
    type Pure = PureRecord<PureModelTest>;
    type Binary =  BinaryRecord<BinaryModelTest>;
    // type Options = ();

    #[test]
    fn file_not_found(){

        let result = p_from_file::<Pure, &str>(&vec!["none"],&"A").unwrap_err();
        let err: &std::io::Error = result.downcast_ref().unwrap();
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);

    }

    // #[cfg(feature = "cpa")]
    #[test]
    fn comp_not_found(){
        // use crate::models::cpa::parameters::CPAPureRecord;
        //  PureRecord<CPAPureRecord>
        let path = "src/parameters/dummy_data/parameters1.json";

        let err:Box<dyn Error> = p_from_file::<Pure, &str>(
            &vec!["methanol", "ethanol", "carbon dioxide"],
            &path
        ).unwrap_err();

        let err:Option<&RecordError> = err.downcast_ref();
        assert_eq!(err.unwrap(), &RecordError::NotFound("methanol, ethanol".into(), path.into()));

    }



    // #[cfg(feature = "cpaa")]
    #[test]
    fn pure_jsons_unmatch_sets(){

        let names1 = vec!["water", "carbon dioxide"];
        let names2 = vec!["methane"];
        let sets = [names1, names2];
        let ppath1 = "src/parameters/dummy_data/parameters1.json";

        let err = p_from_files::<Pure, &str>(&sets,&[ppath1]).unwrap_err();

        assert!(err.len() == 1);

        let err: &RecordError = err[0].downcast_ref().unwrap();
        
        assert_eq!(err, &RecordError::PureJsonsUnmatchedSets(1, 2));

        // dbg!(err);
        // println!("{err}");

    }

    #[test]
    fn binary_jsons_unmatch_sets(){

        let names1 = vec!["water","carbon dioxide"];
        let names2 = vec!["methane"];
        let sets = &[names1, names2];

        let ppath1 = "src/parameters/dummy_data/parameters3.json";

        let err = b_from_files::<Binary, &str>(sets,&[ppath1, ppath1, ppath1]).unwrap_err();

        assert!(err.len() == 1);
        let err: &RecordError = err[0].downcast_ref().unwrap();
        
        assert_eq!(err, &RecordError::BinaryJsonsUnmatchedSets(3, 2));

        // dbg!(err);
        // println!("{err}");

    }


    // #[cfg(feature = "cpa")]
    #[test]
    fn from_multiple_jsons_comps_not_found(){

        let names1 = vec!["water", "carbon dioxide"];
        let names2 = vec!["methane", "propane", "ethane","acetic acid"];
        
        let names = &[names1, names2];

        let ppath1 = "src/parameters/dummy_data/parameters1.json";
        let ppath2 = "src/parameters/dummy_data/parameters2.json";
        

        let res = ParametersTest
        ::from_multiple_jsons(names, &[ppath1,ppath2], None, ());

        // let res = CPAParameters::from_multiple_jsons(&names, &[ppath1,ppath2], &[], CubicModels::default());

        if let Err(e) = res {
            
            let s = format!("{:?}",e );
            // dbg!(e);
            println!("{s}");

        } else { panic!("should err")}

    }

    #[test]
    fn from_n_sets() {

        let names1 = vec!["water"];
        let names2 = vec!["carbon dioxide"];
        let names3 = vec!["methane"];

        let ppath = "src/parameters/dummy_data/parameters1.json";
        let bpath = "src/parameters/dummy_data/bin1.json";

        let p = ParametersTest::from_multiple_jsons(
            &[names1, names2, names3], 
            &[ppath,ppath,ppath], 
            Some(&[bpath]), 
            ()
        ).unwrap();
        
        let water = "Test(x=1, y=2, z=3)".to_string();
        let carbon_dioxide = "Test(x=4, y=5, z=6)".to_string();
        let methane = "Test(x=7, y=8, z=9)".to_string();
        
        let target = vec![water, carbon_dioxide, methane];

        let pure = &p.pure;
        let bin = &p.bin;

        for i in 0..pure.len(){
            let s = format!("{}", pure[i]);
            assert_eq!(s, target[i]);
        }

        
        let water_co2 = "Test(kij=1, lij=2)".to_string();
        let water_methane = "Test(kij=3, lij=4)".to_string();
        let co2_methane = "Test(kij=5, lij=6)".to_string();
        let target = vec![water_co2, water_methane, co2_methane];

        assert!(bin.len() == 3);

        for t in target {
            
            let (_, b) = bin.iter().find(|&(_, b)| format!("{}", b) == t).unwrap();
            assert!(format!("{}", b) == t)
        }

    }
    
    #[test]
    fn from_1_set() {

        let names1 = vec!["water"];
        let names2 = vec!["carbon dioxide"];
        let names3 = vec!["methane"];

        let sets = &[names1, names2, names3];
        let set = sets.concat();
        
        let ppath = "src/parameters/dummy_data/parameters1.json";
        let bpath = "src/parameters/dummy_data/bin1.json";
        
        let p = ParametersTest::from_multiple_jsons(
            &[set], 
            &[ppath], 
            Some(&[bpath]), 
            ()
        ).unwrap();
        
        let water = "Test(x=1, y=2, z=3)".to_string();
        let carbon_dioxide = "Test(x=4, y=5, z=6)".to_string();
        let methane = "Test(x=7, y=8, z=9)".to_string();
        
        let target = vec![water, carbon_dioxide, methane];

        let pure = &p.pure;
        let bin = &p.bin;
        for i in 0..pure.len(){
            let s = format!("{}", pure[i]);
            assert_eq!(s, target[i]);
        }

        let water_co2 = "Test(kij=1, lij=2)".to_string();
        let water_methane = "Test(kij=3, lij=4)".to_string();
        let co2_methane = "Test(kij=5, lij=6)".to_string();
        let target = vec![water_co2, water_methane, co2_methane];

        assert!(bin.len() == 3);

        for t in target {
            
            let (_, b) = bin.iter().find(|&(_, b)| format!("{}", b) == t).unwrap();
            assert!(format!("{}", b) == t)
        }

    }
    

}
