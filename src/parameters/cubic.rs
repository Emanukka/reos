use core::fmt;
use std::collections::HashMap;

use ndarray::Array1;
use serde::{Deserialize, Serialize};

use crate::parameters::Parameters;



#[derive(Clone,Debug)]
pub struct CubicParameters{
    pub ncomp:  usize,
    pub va0:    Array1<f64>,
    pub vb:     Array1<f64>,
    pub vkappa: Array1<f64>,
    pub vtc:    Array1<f64>,
    pub pbinary:HashMap<(usize,usize),(f64,f64)>

}
#[derive(Clone,Debug,Serialize, Deserialize)]

pub struct CubicPureRecord{
    pub a0:f64,
    pub b:f64,
    pub c1:f64,
    pub tc:f64
}

impl CubicPureRecord {
    
    pub fn new(a0:f64,b:f64,c1:f64,tc:f64)->Self{
        Self{
            a0,b,c1,tc  
        }
    }
}
impl CubicParameters {
    
    pub fn from_records(records:Vec<CubicPureRecord>)->Self{

        let ncomp =records.len();

        let mut va0 = Array1::<f64>::zeros(ncomp);
        let mut vb = Array1::<f64>::zeros(ncomp);
        let mut vtc = Array1::<f64>::zeros(ncomp);
        let mut vkappa = Array1::<f64>::zeros(ncomp);

        let pbinary:HashMap<(usize,usize),(f64,f64)>= HashMap::new();

        for (i,record) in records.iter().enumerate(){
            va0[i]=record.a0;
            vb[i] =record.b;
            vkappa[i]=record.c1;
            vtc[i]=record.tc;
        }    
        Self{
            ncomp,
            va0,
            vb:vb.clone(),
            vkappa,
            vtc,
            pbinary
        }

    }

    pub fn set_binary(&mut self,i:usize,j:usize,kij_a:Option<f64>,kij_b:f64){

        self.pbinary.insert((i,j), (kij_a.unwrap_or(0.0),kij_b));
    } 
}


// impl Parameters for CubicParameters {
    
    
//     fn new(data:JsonStruct,comps:Vec<&str>)->Self {

//         let (components,binary)=Self::get_components_and_binary_parameters(data, &comps);
//         let ncomp = comps.len();

//         // Cubic Parameters
//         let mut va0 = Array1::<f64>::zeros(ncomp);
//         let mut vb = Array1::<f64>::zeros(ncomp);
//         let mut vtc = Array1::<f64>::zeros(ncomp);
//         let mut vkappa = Array1::<f64>::zeros(ncomp);

//         let mut pbinary:HashMap<(usize,usize),(f64,f64)>= HashMap::new();

//         for (i,comp) in components.iter().enumerate(){
//             va0[i]=comp.cubic.a0;
//             vb[i]=comp.cubic.b;
//             vkappa[i]=comp.cubic.c1;
//             vtc[i]=comp.cubic.tc;
//         }

//         for (i,j) in binary.keys(){

//             // if kij_t==true
//             let ij=&(*i,*j);
//             let ji=&(*j,*i);
            
//             let opt= binary.get(ij).or_else(||binary.get(ji)).unwrap();

//             let a=opt.kij_a;
//             let b=opt.kij_b;

//             if a.is_some()||b.is_some(){
//                 assert!(a.is_some(),"Interaction k{i}{j}: you forget a parameter");
//                 assert!(b.is_some(),"Interaction k{i}{j}: you forget b parameter"); 
//                 pbinary.insert((*i,*j), (a.unwrap(),b.unwrap())); 
//             }

            
//         }

        
//         Self{
//             ncomp,
//             va0,
//             vb:vb.clone(),
//             vkappa,
//             vtc,
//             pbinary
//         }
//     }
// }


impl fmt::Display for CubicParameters {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {

        writeln!(f, " === Cubic Parameters ===")?;

        writeln!(f, "  Number of Components (ncomp): {}", self.ncomp)?;

        writeln!(f, "  a0 (va0): {:?}", self.va0.to_vec())?;
        writeln!(f, "  b (vb):  {:?}", self.vb.to_vec())?;
        writeln!(f, "  kappa (c1):   {:?}", self.vkappa.to_vec())?;
        writeln!(f, "  T_crit (vtc): {:?}", self.vtc.to_vec())?;

        writeln!(f, "  Binary Interaction Map (kᵢⱼ (a,b) ): {:?}",self.pbinary)?;


        Ok(())
    }
}
