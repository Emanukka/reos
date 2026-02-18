
use std::{error::Error, fmt::Display, vec};

use serde::{Deserialize, Serialize, de::DeserializeOwned};

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



impl<M:Display> Display for PureRecord<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(f, "PureRecord(name={}, molar_weight={}, model_record={})", self.name, self.molar_weight, self.model_record)
        
    }
}

impl<M> PureRecord<M>{
    
    pub fn new<T:ToString>(molar_weight:f64,name: T, model_record: M)->Self{
        
        Self{
            molar_weight,
            model_record,
            name:name.to_string()
        }
    }


}


#[derive(Serialize,Deserialize,Debug,Clone)]
pub struct BinaryRecord<M>{
    #[serde(flatten)]
    pub model_record: M,
    pub id1: String,
    pub id2: String
}

impl<M> BinaryRecord<M>{
    
    pub fn new<T: ToString>(model_record:M, id1:T, id2:T)->Self{
        Self { model_record, id1:id1.to_string(), id2: id2.to_string() }
    }

    pub fn get_id(&self) -> (&str, &str) {
        
        (&self.id1, &self.id2)
    }
}
impl<M:Display> Display for BinaryRecord<M> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(f, "BinaryRecord(id1={}, id2={}, model_record={})", self.id1, self.id2, self.model_record)
        
    }
}

impl<B:Serialize> BinaryRecord<B> {

    pub fn to_json(vec: Vec<Self>) -> Result<String, Box<dyn Error>> {
        let json = serde_json::to_string_pretty(&vec)?;
        Ok(json)
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

#[derive(Debug,Clone, Serialize)]
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


#[cfg(test)] 
mod tests {
    use std::collections::HashMap;

    struct BinaryId(pub usize, pub usize);

    impl PartialEq for BinaryId {
        fn eq(&self, other: &Self) -> bool {

            let (i,j) = (self.0,self.1);
            let (k,l) = (other.0,other.1);

            let left = (i,j) == (k,l);
            let right = (i,j) == (l,k);
            dbg!(left,right,i,j,k,l);
            left || right

        }
    }

    
    impl Eq for BinaryId {}
    #[derive(Clone, PartialEq)]
    struct MyBinary(f64);
    
    impl BinaryId {
        fn from_ij(i:usize,j:usize)->Self{
            Self(i, j)
        }
    }
    // struct BinaryParameter<M>{
    //     pub model_record: M,
    //     pub id: BinaryId
    // }

    // impl<M:Clone> BinaryParameter<M>{

    //     pub fn new(model_record:M, id1:usize,id2:usize,)->Self{
    //         BinaryParameter{
    //             model_record,
    //             id: BinaryId(id1, id2)
    //         }
    //     }
    // }
    // #[test]
    // fn hash() {

    //     let mut hash: HashMap<BinaryId, MyBinary> = HashMap::new();
    //     hash.insert(BinaryId::from_ij(0, 1), MyBinary(1.));

    //     let left = hash.get(&BinaryId(1, 0)).unwrap();
    //     // assert!() == hash.get(&BinaryId(1, 0)  ).unwrap() );
    //     // hash.insert(BinaryId::from_ij(0, 1), MyBinary(1.));
    //     // dbg!(hash.)
    //     // let bin1 = BinaryParameter::new(MyBinary(1.), 0, 1);
    //     // let bin2 = BinaryParameter::new(MyBinary(1.), 1, 0);
        
    //     // let mut set = HashSet::new();
        
    //     // set.p
    // }
}



