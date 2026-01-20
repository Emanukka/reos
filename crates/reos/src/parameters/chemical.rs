// use serde::{Deserialize, Serialize};
// #[derive(Serialize, Deserialize)]
// pub struct ChemicalRecord{
//     pub cas_number: String,
//     pub name: String,
//     pub molar_weight: f64,
    
// }
// impl ChemicalRecord{
//     pub fn new(cas_number:String, name:String, molar_weight:f64)->Self{
//         Self { cas_number, name, molar_weight }
//     }
// }

// impl std::fmt::Display for ChemicalRecord{
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         write!(f, "ChemicalRecord(cas_number={}, name={}, molar_weight={})", self.cas_number, self.name, self.molar_weight)
//     }
// }
