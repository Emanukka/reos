pub mod records;
pub mod chemical;
pub mod reader;
pub mod writer;

use std::{collections::{HashMap, HashSet}, error::Error, fmt::Display, vec};


use serde::de::DeserializeOwned;

pub use crate::parameters::records::{BinaryParameter, BinaryRecord, Properties, PureRecord};


pub trait Parameters: Display {

// pub trait Parameters<Self::Pure:DeserializeOwned + Clone, B:DeserializeOwned + Clone, T>: Display {
    type Pure: DeserializeOwned + Clone;
    type Binary: DeserializeOwned + Clone ;
    type Options;
    
    fn from_raw(pure:Vec<Self::Pure>, binary: BinaryMap<Self::Binary>, properties: Option<Properties>, opt: Self::Options) ->  Result<Self, Box<dyn Error>> 
    where Self: Sized ;


    fn new(pure_records:Vec<PureRecord<Self::Pure>>, binary_records: Vec<BinaryRecord<Self::Binary>>, opt: Self::Options) ->  Result<Self, Box<dyn Error>> 
        where Self: Sized  {

        let component_map = self::component_map(&pure_records);
        let properties = self::properties(&pure_records);

        let binary = self::binary_map(binary_records, component_map);
        let pure = self::only_model_record(pure_records);

        Self::from_raw(pure, binary, Some(properties), opt)
    }

    fn from_json<A: AsRef<str>>(names: &[A], ppath:&A, bpath:Option<&A>, opt: Self::Options) -> Result<Self, Box<dyn Error>> 
        where  Self: Sized {
        
        let pure_records = reader::p_from_file::<Self::Pure, A>(names, ppath)?;
        
        let binary_records:Vec<BinaryRecord<Self::Binary>>;

        if let Some(bpath) = bpath {
            
            binary_records = reader::b_from_file::<Self::Binary, A>(names, bpath)?
        
        } else {

            binary_records = Vec::with_capacity(0)    
        
        }

        // let p = Self::new(pure_records, binary_records, opt);

        Self::new(pure_records, binary_records, opt)
        // Ok(p)
    }

    fn from_multiple_jsons<A: AsRef<str> + Clone>(sets: &[Vec<A>], ppaths: &[A], bpaths: Option<&[A]>, opt: Self::Options)-> Result<Self, Vec<Box<dyn Error>>> where Self : Sized{
     
        let pure_records = reader::p_from_files::<Self::Pure, A>(sets, ppaths)?;
        let binary_records:Vec<BinaryRecord<Self::Binary>>;

        if let Some(bpaths) = bpaths {
            
            binary_records = reader::b_from_files::<Self::Binary, A>(sets, bpaths)?
        
        } else {

            binary_records = Vec::with_capacity(0)    
        
        }

        let res = Self::new(pure_records, binary_records, opt);

        match res {
            Ok(ok) => Ok(ok),
            Err(e) => {
                Err(vec![e])
            }
        }

        // Ok(p)
    }


    fn build_pure(pure_record:PureRecord<Self::Pure>, options:Self::Options) -> Result<Self,Box<dyn Error>> where Self: Sized{

        Self::new(vec![pure_record], vec![], options )

    }

}



fn component_map<M>(records: &Vec<PureRecord<M>>)->HashMap<String,usize>{
    records.iter().enumerate().map(|(i,r)|{
        (r.name.clone(), i)
    }).collect()
}


pub type BinaryMap<B> = HashMap<(usize,usize), B>;



fn binary_map<B>(records: Vec<BinaryRecord<B>>, component_map: HashMap<String,usize>)-> BinaryMap<B>{

    records.into_iter().map(|r|{
        
        let i = *component_map.get(&r.id1).expect("binary should have names belonging to pure records!");
        let j = *component_map.get(&r.id2).expect("binary should have names belonging to pure records!");
        let key: (usize,usize);
        if i < j{
            key = (i, j)

        } else if i > j {
            key = (j, i)    

        } else {panic!("comp_i = comp_j ?!")}

        (key, r.model_record)

    }).collect()
    // v
}


fn only_model_record<M>(pure_records:Vec<PureRecord<M>>) -> Vec<M> {
    pure_records.into_iter().map(|r|{
        r.model_record
    }).collect()
}

fn properties<M>(pure_records: &Vec<PureRecord<M>>)->Properties{

        let n = pure_records.len();
        let mut names = Vec::with_capacity(n);
        let mut molar_weight = Vec::with_capacity(n);
        pure_records.iter().for_each(|r| {
            names.push(r.name.clone());
            molar_weight.push(r.molar_weight);

        });

        Properties { names, molar_weight: ndarray::Array1::from_vec(molar_weight) }
}



