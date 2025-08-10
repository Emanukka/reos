// imports
use ndarray::Array1;
use std::{collections::HashMap, io::Read};
use std::fs::File;

use serde::{Deserialize, Serialize};

use crate::models::cpa::parameters::{AssociationPureRecord, AssociationRule};
use crate::models::cubic::parameters::CubicPureRecord;

// use crate::eos::associative::{AssociationRule, SchemeType};
// use crate::eos::associative::{AssociationRule, SchemeType};

#[derive(Serialize, Deserialize,Debug,Clone)]
pub struct ParametersStruct{

    pub a0:f64,
    pub b:f64,
    pub tc:f64,
    pub c1:f64,
    pub epsilon: Option<f64>,
    pub beta: Option<f64>,
    pub scheme: Option<String>
}

#[derive(Serialize, Deserialize,Debug,Clone)]

pub struct CubicParameters{
    pub a0:f64,
    pub b: f64,
    pub tc:f64,
    pub c1:f64,
}
// #[derive(Serialize, Deserialize,Debug,Clone)]

// pub struct AssociativeParameters{

//     pub epsilon: Option<f64>,
//     pub beta:    Option<f64>,
//     pub scheme:  Option<String>
// }
#[derive(Serialize, Deserialize,Debug,Clone,PartialEq)]
#[serde(rename_all = "lowercase")]

pub enum AssociativeTerm{
    Inert,
    Solvate{
        #[serde(default="site_default")]
        na:f64,
        #[serde(default="site_default")]
        nb:f64,
        #[serde(default="site_default")]
        nc:f64},
    Associative{epsilon:f64,beta:f64,
        #[serde(default="site_default")]
        na:f64,
        #[serde(default="site_default")]
        nb:f64,
        #[serde(default="site_default")]
        nc:f64}
}

fn site_default()->f64{0.0}

#[derive(Serialize, Deserialize,Debug,Clone)]

pub struct Component{
    pub name:String,
    pub cubic: CubicParameters,
    pub assoc: Option<AssociativeTerm>,
}

#[derive(Serialize, Deserialize,Debug,Clone)]

pub struct CompRecord{

    pub name:String,
    pub cubic:CubicPureRecord,
    pub assoc:AssociationPureRecord,

}

impl CompRecord {
    
    pub fn map(json_path:&str)->HashMap<String,CompRecord>{
    // fn components(json_file_name:&str)->HashMap<String,CompRecord>{

        let mut path = std::env::current_dir().unwrap();
        path.push(json_path);
        // path.push("p.json");

        dbg!(&path);
        let mut file: File = File::open(path).unwrap();

        let mut content = String::new();

        file.read_to_string(&mut content).unwrap(); 
        
        let s = content.as_str();

        // let data:HashMap<String,CompRecord> = serde_json::from_str(s).unwrap();
        let data:HashMap<String,CompRecord> = serde_json::from_str(s).unwrap();
        // dbg!(data);
        data        
    }

}
pub struct BinaryRecord{
    pub assoc:AssociativeBinary,
}



#[derive(Serialize, Deserialize,Debug,Clone,PartialEq)]


#[serde(rename_all = "lowercase")]
pub enum AssociativeBinary{
    ///Normal:Interaction between two self-associative componentes.
    /// $a$ and $b$ are binary parameters: $\varepsilon^{\text{cross}} = \sqrt{\varepsilon_i \varepsilon_j} (1 - l_{ij})$, where $l_{ij} = aT + b$
    /// if $a=0$, then $l_{ij}=b=cte$
    Normal{rule:AssociationRule,a:Option<f64>,b:Option<f64>},

    Solvation{rule:AssociationRule,epsilon_cross:Option<f64>,beta_cross:f64}
}
#[derive(Serialize, Deserialize,Debug,Clone)]
pub struct Binary{
    
    /// kij a and b
    pub kij_a:  Option<f64>,
    pub kij_b:  Option<f64>,
    /// Associative Binary Parameters
    pub assoc: Option<AssociativeBinary>,

}
#[derive(Serialize, Deserialize,Debug,Clone)]
pub struct JsonStruct{
    pub components: HashMap<String,Component>,
    pub binary: HashMap<String,Binary>,
}

// pub struct Mixture{

//     components:Vec<Component>
// }


pub trait Parameters{

    fn new(data:JsonStruct,comps:Vec<&str>)->Self;

    fn from_json(json_str: &str,comps:Vec<&str>)->Self where Self: Sized{
        // let data = Self::from_json(json_str, comps.clone());
        let data: JsonStruct = serde_json::from_str(json_str).unwrap();
        // dbg!(data);
        Self::new(data, comps)
    }
    fn from_file(path:&str,comps:Vec<&str>)->Self where Self: Sized{
        let data = Self::get_data(path);
        Self::new(data, comps)

    }
    fn change_binary_keys(comps:&Array1<Component>,binary:&HashMap<String,Binary>,cross_map:&mut HashMap<(usize,usize),Binary>){

        for i in 0..comps.len(){
            for j in (i+1)..comps.len(){

                let name_i=comps[i].name.as_str();                    
                let name_j=comps[j].name.as_str();                    
                
                let i_and_j = name_i.to_owned() + "_" + name_j;
                let j_and_i = name_j.to_owned() + "_" + name_i;

        
                if binary.contains_key(&i_and_j){
        
                    let bin= binary.get(&i_and_j).unwrap().clone();
                    cross_map.insert((i,j), bin);
        
                }else if binary.contains_key(&j_and_i) {
        
                    let bin= binary.get(&j_and_i).unwrap().clone();
                    cross_map.insert((i,j), bin);
                }
            }
        }
    }

    fn get_data(dir:&str)->JsonStruct{

        let mut path = std::env::current_dir().unwrap();


        path.push(dir);
        // path.push("p.json");


        dbg!(&path);
        let mut file: File = File::open(path).unwrap();

        let mut content = String::new();

        file.read_to_string(&mut content).unwrap(); 
        
        let s = content.as_str();

        let data: JsonStruct = serde_json::from_str(s).unwrap();
        // dbg!(data);
        data
    }

    fn get_components_and_binary_parameters(data:JsonStruct,comps:&Vec<&str>)->(Array1<Component>,HashMap<(usize,usize),Binary>){

        // Get the data inside the JSON's file
        let ncomp = comps.len();
        let mut all_comps_from_json: Vec<Component>= Vec::with_capacity(ncomp);
        let components= data.components;
        let binary= data.binary;
        
        for i in 0..ncomp{

            let uppercase_i = comps[i].to_uppercase();
            let name_i = uppercase_i.as_str();

            if components.contains_key(name_i){

                let mut comp_i= components.get(name_i).unwrap().clone();
                comp_i.name= uppercase_i;

                all_comps_from_json.push(comp_i);

            }else {
                panic!("Component {name_i} not found.")
            }
        
        }

        let comps = Array1::from_vec(all_comps_from_json);
        let mut cross_map:HashMap<(usize,usize), Binary> = HashMap::new();

        // <String,Binary> --> <(i,j),binary>
        if ncomp>1{ Self::change_binary_keys(&comps,&binary,&mut cross_map) }

        (comps,cross_map)
    }



}
    
