
use serde::{Deserialize, Serialize};

use ndarray::{Array1, array};

#[derive(Serialize,Deserialize, Debug, Clone)]
pub struct PureRecord<M>{
    #[serde(default)]
    pub molar_weight: f64, //
    #[serde(default)]
    pub name: String,
    #[serde(flatten)]
    pub model_record: M,

}


impl<M> PureRecord<M>{
    
    pub fn new(molar_weight:f64,name: String, model_record: M)->Self{
        
        Self{
            molar_weight,
            model_record,
            name
        }
    }

}

#[derive(Serialize,Deserialize,Debug)]
pub struct BinaryRecord<M>{
    #[serde(flatten)]
    pub model_record: M,
    pub id1: String,
    pub id2: String
}

impl<M> BinaryRecord<M>{
    
    pub fn new(model_record:M, id1:String, id2:String)->Self{
        Self { model_record, id1, id2 }
    }

    pub fn get_id(&self) -> (&str, &str) {
        
        (&self.id1, &self.id2)
    }
}
#[derive(Debug)]
pub struct BinaryParameter<M>{
    pub model_record: M,
    pub id1: usize,
    pub id2: usize
}

impl<M:Clone> BinaryParameter<M>{

    pub fn new(model_record:M,id1:usize,id2:usize,)->Self{
        BinaryParameter{
            model_record,
            id1,
            id2
        }
    }
}

#[derive(Debug)]
pub struct Properties{
    pub names: Vec<String>,
    pub molar_weight: Array1<f64>
    //other things...
}

impl Default for Properties{

    fn default() -> Self {
        Properties { names: vec![], molar_weight: array![] }
    }
}

