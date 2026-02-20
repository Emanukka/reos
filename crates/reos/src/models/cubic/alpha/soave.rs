use serde::{Deserialize, Serialize};

use crate::models::cubic::{alpha::{AlphaError, AlphaModel, AlphaRecord, Alpha}, models::{CubicModel, CubicModels}, parameters};


#[derive(Clone,Debug, PartialEq)]
pub struct Soave(pub Vec<f64>);

impl AlphaModel for Soave {

    fn tag() -> Self {
        Soave(vec![])
    }
    
    fn build<A: AsRef<str>>(names:&[A], records:Vec<AlphaRecord>, model: &CubicModels) -> Result<Alpha, AlphaError> {
        // let mut errors = vec![];
        let mut errors = vec![];
        let mut kappa = Vec::with_capacity(records.len());

        records.into_iter().enumerate().for_each(|(i, r)| {
            
             if let AlphaRecord::Soave { w } = r {

                kappa.push(model.kappa_from_w(w))
            }

            else if let AlphaRecord::SoaveRegressed { c1 } = r {

                kappa.push(c1)

            } else {
                errors.push(names[i].as_ref());

                // Err(AlphaError::AssociatedRecord(names[i].as_ref().to_string()))
            }
        });

        // .filter_map(|res| res.map_err(| e| errors.push(e)).ok() )
        // .collect();

        if errors.len() != 0 {

            let e = errors.join(",");
            Err(AlphaError::AssociatedRecord(e))
        
        } else {

            Ok(Alpha::Soave(Soave(kappa)))
        }
    }

    fn alpha(&self,t:f64, tc:&[f64]) -> Vec<f64>{

        let kappa = &self.0;
        tc.iter().enumerate().map(|(i, tc)| {
            
            let tr = t / tc;
            (1.0 + kappa[i] * (1.0 - tr.sqrt())).powi(2)

        })
        .collect()
    }

    fn dalpha_dt(&self,t:f64,tc:&[f64]) -> Vec<f64> {
        
        let kappa = &self.0;
        let alpha = self.alpha(t, tc);

        tc.iter().enumerate().map(|(i, tc)| {
            
            let tr = t / tc;
            - kappa[i] * (alpha[i] * tr).sqrt() / t

        })
        .collect()
    }
    
    fn to_string(&self) -> String {

        format!("kappa={:?}", self.0)
    
    }

}

// #[cfg(test)]
// mod tests {
//     use serde::Deserialize;

//     use crate::models::cubic::{_alpha::{AlphaModel, AlphaOption, AlphaRecord, Alpha, CubicOptions, soave::SoaveFlavor}, models::{CubicModel, PR78}};



//     #[test]
//     fn options(){

//         let names = ["foo"];
//         let alpha = AlphaOption::Soave;
//         let model = PR78.into();
//         let options = CubicOptions::new(alpha, model);
        
//         let w = 0.275000;
        
//         let r = AlphaRecord::Soave(SoaveFlavor::Omega(w));

//         let alpha = Alpha::build(&names, vec![r], &options).unwrap();

//         let s = alpha.to_string();

//         assert_eq!(s, format!("kappa=[{}]",options.model.kappa_from_w(w)));

//     }

//     #[test]
//     fn from_json(){

//         let s = r#"
//             {
//                 "w": 0.275000
//             }
//         "#;

//         let r:AlphaRecord = serde_json::from_str(s).unwrap();

//         let names = ["foo"];
//         let alpha = AlphaOption::Soave;
//         let model = PR78.into();
//         let options = CubicOptions::new(alpha, model);
        
//         let w = 0.275000;
        
//         // let r = AlphaRecord::Soave(SoaveFlavor::Omega(w));

//         let alpha = Alpha::build(&names, vec![r], &options).unwrap();

//         let s = alpha.to_string();

//         assert_eq!(s, format!("kappa=[{}]",options.model.kappa_from_w(w)));

//     }
// }
// // 