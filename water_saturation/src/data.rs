use std::collections::HashMap;

use reeos::parameters::CompRecord;
use serde::{Deserialize, Serialize};


#[derive(Deserialize,Serialize,Debug)]

pub struct Input{

    pub composição:HashMap<String,f64>,
    pub temperatura:f64,
    pub pressão:Vec<f64>,

}
#[derive(Deserialize,Serialize,Debug)]

pub enum Component {
    CO2(f64),
    CH4(f64),
    N2 (f64),
    H2S(f64),
    C2 (f64),
    C3 (f64),
    C4 (f64),
    IC4(f64),
    C5 (f64),
    IC5(f64),
    C6 (f64),
    C7 (f64),
    C8 (f64),
    C9 (f64),
    C10(f64),

}


//
// impl Component {
    
//     fn get_record(&self)->CompRecord{

//         match self {
            
//             Self::CO2(_)=>get_co2()
//         }
        
        
//     }
// }