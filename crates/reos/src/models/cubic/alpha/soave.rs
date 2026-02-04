use serde::{Deserialize, Serialize};

use crate::models::cubic::{alpha::{AlphaError, AlphaModel, AlphaRecord, Alphaa, CubicOptions}, models::CubicModel};


#[derive(Serialize, Deserialize, Clone, Debug)]
pub enum SoaveFlavor {
    #[serde(rename = "c1")]
    Kappa(f64),
    #[serde(rename = "w")]
    Omega(f64),
}
#[derive(Clone,Debug)]
pub struct Soave(pub Vec<f64>);

impl Soave {

    pub fn from_records<A: AsRef<str>>(names:&[A], records:Vec<AlphaRecord>, options: &CubicOptions)->Result<Alphaa, Vec<AlphaError>> {

        let cubic = &options.model;
        let mut errors = vec![];
        let kappa = records.into_iter().enumerate().map(|(i, r)| {
            
             if let AlphaRecord::Soave(flavor) = r {

                match flavor {
                    SoaveFlavor::Kappa(kp) => Ok(kp),
                    SoaveFlavor::Omega(w) => Ok(cubic.kappa_from_w(w))
                }
                
            } else {

                Err(AlphaError::AssociatedRecord(names[i].as_ref().to_string()))
            }
        })
        .filter_map(|res| res.map_err(| e| errors.push(e)).ok() )
        .collect();

        if errors.len() != 0 {

            Err(errors)
        
        } else {

            Ok(Soave(kappa).into())
        }
     
    }

}

impl AlphaModel for Soave {

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

#[cfg(test)]
mod tests {
    use serde::Deserialize;

    use crate::models::cubic::{alpha::{AlphaModel, AlphaOption, AlphaRecord, Alphaa, CubicOptions, soave::SoaveFlavor}, models::{CubicModel, PR78}};



    #[test]
    fn options(){

        let names = ["foo"];
        let alpha = AlphaOption::Soave;
        let model = PR78.into();
        let options = CubicOptions::new(alpha, model);
        
        let w = 0.275000;
        
        let r = AlphaRecord::Soave(SoaveFlavor::Omega(w));

        let alphaa = Alphaa::build(&names, vec![r], &options).unwrap();

        let s = alphaa.to_string();

        assert_eq!(s, format!("kappa=[{}]",options.model.kappa_from_w(w)));

    }

    #[test]
    fn from_json(){

        let s = r#"
            {
                "w": 0.275000
            }
        "#;

        let r:AlphaRecord = serde_json::from_str(s).unwrap();

        let names = ["foo"];
        let alpha = AlphaOption::Soave;
        let model = PR78.into();
        let options = CubicOptions::new(alpha, model);
        
        let w = 0.275000;
        
        // let r = AlphaRecord::Soave(SoaveFlavor::Omega(w));

        let alphaa = Alphaa::build(&names, vec![r], &options).unwrap();

        let s = alphaa.to_string();

        assert_eq!(s, format!("kappa=[{}]",options.model.kappa_from_w(w)));

    }
}
// 