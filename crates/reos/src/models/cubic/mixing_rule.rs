//! Return w, the vector of mixture parameters (B, D, d1, d2)
//! and the derivatives in respect to state variables (n, T, V)


use enum_dispatch::enum_dispatch;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

use crate::models::cubic::{DwDni, DwDt, W, combining_rule::CombiningRuleModel, models::CubicModel, parameters::CubicParameters};


#[enum_dispatch]
pub trait MixingRuleModel: Clone + std::fmt::Debug{

    fn apply(&self, t:f64, v:f64, x:&Array1<f64>, parameters:&CubicParameters) -> W ;

    fn dw_dt(&self, t:f64, v:f64, x:&Array1<f64>, w:&W, parameters:&CubicParameters) -> DwDt ;

    fn dw_dni(&self, t:f64, v:f64, x:&Array1<f64>, w:&W, parameters:&CubicParameters) -> DwDni ;

    fn to_string(&self)->String;
}

#[derive(Clone,Debug, PartialEq)]
pub struct Quadratic;



impl MixingRuleModel for  Quadratic {

    fn to_string(&self) -> String {
        "Quadratic".to_string()
    }
    
    fn apply(&self, t:f64, _:f64, x:&Array1<f64>, parameters:&CubicParameters) -> W {
        
        let n = parameters.tc.len();
        let [aij_mat, bij_mat] = [&parameters.aij, &parameters.bij];
        // let epsilon = parameters.model.eps();  
        // let sigma = parameters.model.sig();  
        // let bin = &parameters.binary;
        let kij_mat = &parameters.kij;
        let alpha = &parameters.alpha.alpha(t, &parameters.tc);

        // let [ac, bc, _] = [&parameters.a,&parameters.b,&parameters.c];
        let [mut d,mut b, _] = [0.,0.,0.];
        let combr = &parameters.combr;

        for i in 0..n {
            
            // c += x[i] * vt[i];

            for j in 0..n {

                // let bi_bj = [bc[i], bc[j]];
                // let ai_aj = [aij[(i,i)] * alpha[i], ac[j] * alpha[j]];
                
                let alpha_ij = combr.alpha_ij(alpha[i], alpha[j]);

                let kij = kij_mat[(i,j)].a + kij_mat[(i,j)].b * t;
                // let [aij, bij] = combr.apply(ai_aj, bi_bj, kij);
                

                b += x[i] * x[j] * bij_mat[(i,j)];
                d += x[i] * x[j] * aij_mat[(i,j)] * alpha_ij * (1.0 - kij);
                

            }
        }

        let d1 = parameters.model.eps();
        let d2 = parameters.model.sig();
        W{b, d, d1, d2}

    }

    fn dw_dt(&self, t:f64, _:f64, x:&Array1<f64>, _:&W, parameters:&CubicParameters) -> DwDt{

        let kij_mat = &parameters.kij;
        let alpha = &parameters.alpha.alpha(t, &parameters.tc);
        let dalpha_dt = parameters.alpha.dalpha_dt(t, &parameters.tc);

        // let [ac, _, _] = [&parameters.a,&parameters.b,&parameters.c];
        let [ai] = [&parameters.aij.diag()];

        let n = alpha.len();
        let mut dd = 0.;
        let combr = &parameters.combr;


        for i in 0..n {
            for j in 0..n{

                let dkij = kij_mat[(i,j)].b;
                let kij = kij_mat[(i,j)].a + dkij * t;
                 
                let ai_aj = [ai[i] * alpha[i] , ai[j] * alpha[j]];
                let dai_daj = [ai[i] * dalpha_dt[i] , ai[j] * dalpha_dt[j]];

                let daij_dt = combr.dt(ai_aj, dai_daj, kij, dkij);
                // let sqrt = (ai * aj).sqrt();

                // // dd += 0.5 * (1. - kij) * ( dai * aj + daj * ai) / sqrt - sqrt * dkij;
                // let daij_dt = 0.5 * (1. - kij) * ( dai * aj + daj * ai) / sqrt - sqrt * dkij;

                dd += x[i] * x[j] * daij_dt

            }
        }   

        DwDt { d: dd }

    }

    fn dw_dni(&self, t:f64, _:f64, x:&Array1<f64>, w:&W, parameters:&CubicParameters) -> DwDni {

        let n = x.len();
        let kij_mat = &parameters.kij;
        let alpha = &parameters.alpha.alpha(t, &parameters.tc);
        // let [ac, bc, _] = [&parameters.a,&parameters.b,&parameters.c];
        let [aij_mat, bij_mat] = [&parameters.aij, &parameters.bij];

        let combr = &parameters.combr;
        let mut db = Vec::with_capacity(n);
        let mut dd = Vec::with_capacity(n);
        
        for i in 0..n {

            let mut sum_db = 0.; 
            let mut sum_dd = 0.; 

            for j in 0..n {
                
                // let bij = 0.5 * (bc[i] + bc[j]);
                let alpha_ij = combr.alpha_ij(alpha[i], alpha[j]);

                let kij = kij_mat[(i,j)].a + kij_mat[(i,j)].b * t;

                // let [ai, aj] = [ac[i] * alpha[i], ac[j] * alpha[j]];
                // let sqrt = (ai * aj).sqrt();

                // let aij = (1. - kij) * sqrt;

                sum_dd += x[j] * aij_mat[(i, j)] * alpha_ij * (1.0 - kij);
                sum_db += x[j] * bij_mat[(i, j)];

            }

            db.push(2. * sum_db - w.b);
            dd.push(2. * sum_dd)
        }

        DwDni { b: db, d: dd }

    }
}

#[enum_dispatch(MixingRuleModel)]
#[derive(Clone,Debug, PartialEq)]
pub enum MixingRule{
    Quadratic,
}

#[derive(Clone,Debug,Serialize,Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]

pub enum MixingRuleOption {
    Quadratic,
}

impl Default for MixingRuleOption {
    
    fn default() -> Self {
        Self::Quadratic
    }
}

impl From<MixingRuleOption> for MixingRule {
    
    fn from(value: MixingRuleOption) -> Self {
        
        match value {

            MixingRuleOption::Quadratic => Quadratic.into(),

        }    
    }
}


#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn assert_model_parse() {

        let quad:MixingRule = MixingRuleOption::Quadratic.into();

        assert_eq!(quad, MixingRule::Quadratic(Quadratic));
        
    }
}