//! Return w, the vector of mixture parameters (B, D, d1, d2)
//! and the derivatives in respect to state variables (n, T, V)

use std::{fmt::Debug, str::FromStr};

use enum_dispatch::enum_dispatch;
use ndarray::Array1;

use crate::models::cubic::{DwDni, DwDt, W, combining_rule::CombiningRuleModel, models::CubicModel, options::OptionsParseError, parameters::CubicParameters};


#[enum_dispatch]
pub trait MixingRuleModel: Clone + Debug{

    fn apply(&self, t:f64, v:f64, x:&Array1<f64>, parameters:&CubicParameters) -> W ;

    fn dw_dt(&self, t:f64, v:f64, x:&Array1<f64>, w:&W, parameters:&CubicParameters) -> DwDt ;

    fn dw_dni(&self, t:f64, v:f64, x:&Array1<f64>, w:&W, parameters:&CubicParameters) -> DwDni ;

    fn to_string(&self)->String;
}

#[derive(Clone,Debug,)]
pub struct Quadratic;



impl MixingRuleModel for  Quadratic {

    fn to_string(&self) -> String {
        "Quadratic".to_string()
    }
    
    fn apply(&self, t:f64, _:f64, x:&Array1<f64>, parameters:&CubicParameters) -> W {
        let n = parameters.a.len();
        // let epsilon = parameters.options.model.eps();  
        // let sigma = parameters.options.model.sig();  
        let bin = &parameters.binary;
        let alpha = &parameters.options.alpha.alpha(t, &parameters.tc);

        let [ac, bc, _] = [&parameters.a,&parameters.b,&parameters.c];
        let [mut d,mut b, _] = [0.,0.,0.];
        let combr = &parameters.options.combr;

        for i in 0..n {
            
            // c += x[i] * vt[i];

            for j in 0..n {

                let bi_bj = [bc[i], bc[j]];
                let ai_aj = [ac[i] * alpha[i], ac[j] * alpha[j]];
                
                let kij = bin[(i,j)].kij + bin[(i,j)].lij * t;

                let [aij, bij] = combr.apply(ai_aj, bi_bj, kij);
                

                b += x[i] * x[j] * bij;
                d += x[i] * x[j] * aij;
                

            }
        }

        let d1 = parameters.options.model.eps();
        let d2 = parameters.options.model.sig();
        W{b, d, d1, d2}

    }

    fn dw_dt(&self, t:f64, _:f64, x:&Array1<f64>, _:&W, parameters:&CubicParameters) -> DwDt{

        let bin = &parameters.binary;
        let alpha = &parameters.options.alpha.alpha(t, &parameters.tc);
        let dalpha_dt = parameters.options.alpha.dalpha_dt(t, &parameters.tc);

        let [ac, _, _] = [&parameters.a,&parameters.b,&parameters.c];
        let n = alpha.len();
        let mut dd = 0.;
        let combr = &parameters.options.combr;


        for i in 0..n {
            for j in 0..n{

                let dkij = bin[(i,j)].lij;
                let kij = bin[(i,j)].kij + dkij * t;
                 
                let ai_aj = [ac[i] * alpha[i] , ac[j] * alpha[j]];
                let dai_daj = [ac[i] * dalpha_dt[i] , ac[j] * dalpha_dt[j]];

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
        let bin = &parameters.binary;
        let alpha = &parameters.options.alpha.alpha(t, &parameters.tc);
        let [ac, bc, _] = [&parameters.a,&parameters.b,&parameters.c];

        let mut db = Vec::with_capacity(n);
        let mut dd = Vec::with_capacity(n);
        
        for i in 0..n {

            let mut sum_db = 0.; 
            let mut sum_dd = 0.; 

            for j in 0..n {
                
                let bij = 0.5 * (bc[i] + bc[j]);

                let kij = bin[(i,j)].kij + bin[(i,j)].lij * t;

                let [ai, aj] = [ac[i] * alpha[i], ac[j] * alpha[j]];
                let sqrt = (ai * aj).sqrt();

                let aij = (1. - kij) * sqrt;

                sum_dd += x[j] * aij;
                sum_db += x[j] * bij;

            }

            db.push(2. * sum_db - w.b);
            dd.push(2. * sum_dd)
        }

        DwDni { b: db, d: dd }

    }
}

#[enum_dispatch(MixingRuleModel)]
#[derive(Clone,Debug)]
pub enum MixingRule{
    Quadratic,
}

impl Default for MixingRule {

    fn default() -> Self {
        Self::Quadratic(Quadratic)
    }
}


impl FromStr for MixingRule {
    type Err = OptionsParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {

        match s.to_lowercase().as_str() {
            "quadratic" => Ok(Self::Quadratic(Quadratic)),
            "" => Ok(Self::default()),
            _ => Err(OptionsParseError(format!("{} is not a mixing rule implemented",s)))
        }
    }
}