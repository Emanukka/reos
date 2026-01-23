pub mod records;
pub mod chemical;
pub mod reader;
pub mod writer;

use std::{collections::HashMap, error::Error, fmt::Display, vec};


use serde::de::DeserializeOwned;

// use crate::parameters::records::ModelRecords;
pub use crate::parameters::records::{BinaryParameter, BinaryRecord, Properties, PureRecord};


pub trait Parameters<M:DeserializeOwned + Clone, B:DeserializeOwned + Clone, T>: Display {



    /// Raw initializer function for a parameters object.
    /// This function receives processed parameters by the `new()` function.
    fn from_raw(pure:Vec<M>, binary: Vec<BinaryParameter<B>>, properties: Option<Properties>, opt: T) -> Self;


    fn new(pure_records:Vec<PureRecord<M>>, binary_records: Vec<BinaryRecord<B>>, opt: T) -> Self
        where Self: Sized  {

        let component_map = self::component_map(&pure_records);
        let properties = self::properties(&pure_records);

        let binary = self::binary_parameters(binary_records, component_map);
        let pure = self::only_model_record(pure_records);

        Self::from_raw(pure, binary, Some(properties), opt)
    }

    fn from_json<A: AsRef<str>>(names: &[A], ppath:&A, bpath:Option<&A>, opt: T) -> Result<Self, Box<dyn Error>> 
        where  Self: Sized {
        
        let pure_records = reader::p_from_file::<M, A>(names, ppath)?;
        
        let binary_records:Vec<BinaryRecord<B>>;

        if let Some(bpath) = bpath {
            
            binary_records = reader::b_from_file::<B, A>(names, bpath)?
        
        } else {

            binary_records = Vec::with_capacity(0)    
        
        }

        let p = Self::new(pure_records, binary_records, opt);

        Ok(p)
    }

    fn from_multiple_jsons<A: AsRef<str> + Clone>(sets: &[Vec<A>], ppaths: &[A], bpaths: Option<&[A]>, opt: T)-> Result<Self, Vec<Box<dyn Error>>> where Self : Sized{
     
        let pure_records = reader::p_from_files::<M, A>(sets, ppaths)?;
        let binary_records:Vec<BinaryRecord<B>>;

        if let Some(bpaths) = bpaths {
            
            binary_records = reader::b_from_files::<B, A>(sets, bpaths)?
        
        } else {

            binary_records = Vec::with_capacity(0)    
        
        }

        let p = Self::new(pure_records, binary_records, opt);

        Ok(p)
    }
    
}

fn component_map<M>(records: &Vec<PureRecord<M>>)->HashMap<String,usize>{
    records.iter().enumerate().map(|(i,r)|{
        (r.name.clone(), i)
    }).collect()
}


fn binary_parameters<B>(records: Vec<BinaryRecord<B>>,component_map: HashMap<String,usize>)-> Vec<BinaryParameter<B>>{
    let mut v = vec![];
    records.into_iter().for_each(|r|{
        
        let i = *component_map.get(&r.id1).expect("binary should have names belonging to pure records!");
        let j = *component_map.get(&r.id2).expect("binary should have names belonging to pure records!");
        let key: (usize,usize);
        if i < j{
            key = (i, j)
        } else {
            key = (j, i)    
        }
        v.push(
            BinaryParameter { model_record: r.model_record, id1: key.0, id2: key.1 }
        );
    });
    v
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


// struct A;
// struct 
// impl<T:> Parameters<f64,f64,T> for A {
    
//     fn from_records(pure_records: Vec<PureRecord<f64>>, binary_records: Vec<BinaryRecord<f64>>, opt: Option<T>) -> Self {
//         A
//     }

// }
// pub trait Parameters<M,B> 
//     where M: DeserializeOwned, B: Clone + DeserializeOwned
// {





//     fn from_records(pure_records:Vec<PureRecord<M>>,binary_records: Vec<BinaryRecord<B>>) -> ModelParameters<M, B>{

//         let n = pure_records.len();
//         let component_map = Self::component_map(&pure_records);
//         let properties = Self::properties(&pure_records);
//         let binary_map = Self::binary_map(binary_records, component_map);
        
//         let binary = Self::binary_parameters(n, binary_map);
//         let pure = Self::only_model_record(pure_records);

//         // Self::raw(records, binary, Some(properties))
//         ModelParameters { pure, binary, properties: Some(properties) }

//     }


//     // /// Initializer function for a parameters object. 
//     // /// This function receives processed parameters by the `new()` function. Therefore,
//     // /// this function must be implemented for the correspondent parameters object
//     // /// 
//     // /// We separate identifications - names and molar weights - in a Properties object.
//     // /// This enables the parameters object have fields for their specific parameters and 
//     // /// other only for identification, what enables a straightforward and generic implementation
//     // /// of a parameters object.
//     // fn from_raw(
//     //     pure: Vec<M>, 
//     //     binary: Vec<BinaryParameter<B>>, 
//     //     properties: Option<Properties>, 
//     //     model: Option<>)->Self;
    
//     // // fn get_properties(&self)->&Properties;

//     // /// Instantiate a parameters result from json files.
//     // fn mp_from_json(names:&[&str], ppath:&str, bpath:Option<&str>) -> Result<ModelParameters<M, B>, Box<dyn Error>> 
//     //     where  Self: Sized {
        
//     //     let pure_records = reader::p_from_file::<M>(names, ppath)?;
        
//     //     let binary_records:Vec<BinaryRecord<B>>;

//     //     if let Some(bpath) = bpath {
            
//     //         binary_records = reader::b_from_file::<B>(names, bpath)?
        
//     //     } else {

//     //         binary_records = Vec::with_capacity(0)    
        
//     //     }

//     //     let p = Self::mp_from_records(pure_records, binary_records);
//     //     Ok(p)

//     // }


// }


// pub trait ParametersOption{}
