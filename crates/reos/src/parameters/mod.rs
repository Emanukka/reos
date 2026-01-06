pub mod records;
pub mod reader;

use std::{collections::HashMap, error::Error, vec};


use serde::de::DeserializeOwned;

pub use crate::parameters::records::{BinaryParameter, BinaryRecord, Properties, PureRecord};

/// A Parameters trait for residual parameters object
/// 
/// All residual implementations must have a field for parameters, which
/// must implement this trait
/// 
/// This trait have the functionality of provide functions that
/// initialize a generic parameters object given as input components
/// and binary interaction records
///  
pub trait Parameters<M,B> 
    where M: DeserializeOwned, B: Clone + DeserializeOwned
{



    fn component_map(records: &Vec<PureRecord<M>>)->HashMap<String,usize>{

        records.iter().enumerate().map(|(i,r)|{
            (r.name.clone(), i)
        }).collect()
    }

    fn binary_map(records: Vec<BinaryRecord<B>>,component_map: HashMap<String,usize>)->HashMap<(usize,usize),B>{

        records.into_iter().map(|r|{
            
            let i = *component_map.get(&r.id1).expect("binary should have names belonging to pure records!");
            let j = *component_map.get(&r.id2).expect("binary should have names belonging to pure records!");
            
            let key: (usize,usize);

            if i < j{
                key = (i, j)
            } else {
                key = (j, i)    
            }

            (key, r.model_record)

        }).collect()

    }
    fn binary_map2(records: Vec<BinaryRecord<B>>,component_map: HashMap<String,usize>)->HashMap<(usize,usize),B>{

        records.into_iter().map(|r|{
            
            let i = *component_map.get(&r.id1).expect("binary should have names belonging to pure records!");
            let j = *component_map.get(&r.id2).expect("binary should have names belonging to pure records!");
            
            let key: (usize,usize);

            if i < j{
                key = (i, j)
            } else {
                key = (j, i)    
            }

            (key, r.model_record)

        }).collect()

    }

    fn binary_parameters(n:usize,binary_map: HashMap<(usize,usize),B>)->Vec<BinaryParameter<B>>{
        
        let mut v = vec![];
        for i in 0..n {
            for j in i+1..n {  
                
                if let Some(model) = binary_map.get(&(i,j)) {

                    v.push(
                        BinaryParameter { model_record: model.clone(), id1: i, id2: j }
                    )
                }

            }
        }
        v  
    }

    fn only_model_record(pure_records:Vec<PureRecord<M>>) -> Vec<M> {

        pure_records.into_iter().map(|r|{
            r.model_record
        }).collect()
    }

    fn properties(pure_records: &Vec<PureRecord<M>>)->Properties{

        let n = pure_records.len();
        let mut names = Vec::with_capacity(n);
        let mut molar_weight = Vec::with_capacity(n);
        pure_records.iter().for_each(|r| {
            names.push(r.name.clone());
            molar_weight.push(r.molar_weight);

        });

        Properties { names, molar_weight: ndarray::Array1::from_vec(molar_weight) }
    }

    /// Initialize a new parameters object from records.
    fn new(pure_records:Vec<PureRecord<M>>,binary_records: Vec<BinaryRecord<B>>)->Self where Self: Sized{

        let n = pure_records.len();
        let component_map = Self::component_map(&pure_records);
        let properties = Self::properties(&pure_records);

        let binary_map = Self::binary_map(binary_records, component_map);
        let binary = Self::binary_parameters(n, binary_map);
        let records = Self::only_model_record(pure_records);

        Self::raw(records, binary, Some(properties))
    }

    /// Initializer function for a parameters object. 
    /// This function receives processed parameters by the `new()` function. Therefore,
    /// this function must be implemented for the correspondent parameters object
    /// 
    /// We separate identifications - names and molar weights - in a Properties object.
    /// This enables the parameters object have fields for their specific parameters and 
    /// other only for identification, what enables a straightforward and generic implementation
    /// of a parameters object.
    fn raw(records: Vec<M>, binary: Vec<BinaryParameter<B>>, properties: Option<Properties>)->Self;
    
    fn get_properties(&self)->&Properties;

    /// Instantiate a parameters result from json files.
    fn from_json(names:&[&str], ppath:&str, bpath:Option<&str>) -> Result<Self, Box<dyn Error>> where  Self: Sized {
        
        let pure_records = reader::p_from_file::<M>(names, ppath)?;
        
        let binary_records:Vec<BinaryRecord<B>>;

        if let Some(bpath) = bpath {
            
            binary_records = reader::b_from_file::<B>(names, bpath)?
        
        } else {

            binary_records = Vec::with_capacity(0)    
        
        }

        let p = Self::new(pure_records, binary_records);
        Ok(p)

    }


}


