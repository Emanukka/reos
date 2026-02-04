use serde::{Deserialize, Serialize};

use crate::models::cubic::{alpha::{AlphaError, AlphaModel, AlphaRecord, Alphaa, CubicOptions}, models::CubicModel};

#[derive(Clone,Debug)]
pub struct Twu91(pub Vec<[f64;3]>);


impl Twu91 {

    pub fn from_records<A: AsRef<str>>(names:&[A], records:Vec<AlphaRecord>, _: &CubicOptions)->Result<Alphaa, Vec<AlphaError>> {

        let mut errors = vec![];
        let lmn = records.into_iter().enumerate().map(|(i, r)| {
            
             if let AlphaRecord::Twu91 { l, m, n } = r {
                
                Ok([l, m, n])
                
            } else {

                Err(AlphaError::AssociatedRecord(names[i].as_ref().to_string()))
            }
        })
        .filter_map(|res| res.map_err(| e| errors.push(e)).ok() )
        .collect::<Vec<[f64;3]>>();

        if errors.len() != 0 {

            Err(errors)
        
        } else {

            Ok(Twu91(lmn).into())

        }
     
    }

}

impl AlphaModel for Twu91 {

    fn alpha(&self,t:f64, tc:&[f64]) -> Vec<f64>{

        tc.iter().enumerate().map(|(i, tc)| {

            let tr = t / tc;
            let [l, m, n] = self.0[i];
            tr.powf(n * (m - 1.)) * (l * (1. - tr.powf(n * m))).exp()

        })
        .collect()
    }

    fn dalpha_dt(&self,t:f64,tc:&[f64]) -> Vec<f64> {

        let alpha = self.alpha(t, tc);

        tc.iter().enumerate().map(|(i, tc)| {

            let tr = t / tc;
            let [l, m, n] = self.0[i];
            alpha[i] * n * (m - 1. - l * m * tr.powf(n * m)) / t

        })
        .collect()
    }
    fn to_string(&self) -> String {

        format!("[L,M,N]={:?}", self.0)
    
    }

}

#[cfg(test)]
mod tests {
    use approx::assert_relative_eq;

    use crate::models::cubic::{alpha::{AlphaModel, AlphaOption, AlphaRecord, Alphaa, CubicOptions}, models::{PR78}};



    #[test]
    fn options(){

        let alpha = AlphaOption::Twu91;
        let model = PR78.into();
        let options = CubicOptions::new(alpha, model);
        
        let r = AlphaRecord::Twu91 { l:0.28726, m:0.83405 , n:2.01991 }; 
        let alphaa = Alphaa::build(&[""], vec![r], &options).unwrap();

        let a = alphaa.alpha(198.15, &[507.60])[0];
        let s = alphaa.to_string();
        
        println!("{}",s);
        assert_relative_eq!(1.7223476412263239, a, epsilon=1e-9);



    }
}
// 