pub mod cubic;
use std::sync::Arc;

use serde::{Deserialize, Serialize};

use crate::{models::IDEAL_GAS_CONST, parameters::Parameters, residual::Residual, state::State};



pub type ScalingRecord = [f64;6];

// // Vec<f64> , [c1, c2, c3, c4, c5, c6, ..., Cm] ?


// impl Parameters<ScalingRecord, (), ()> for ScalingParameters {

//     fn from_raw(pure:Vec<ScalingRecord>, _: Vec<crate::parameters::BinaryParameter<()>>, _: Option<crate::parameters::Properties>, _: ()) -> Self {
        
//         // let n = pure.len();
//         // let mut a1:Vec<f64> = Vec::with_capacity(n); 
//         // let mut a2:Vec<f64> = Vec::with_capacity(n); 
//         // let mut b1:Vec<f64> = Vec::with_capacity(n); 
//         // let mut b2:Vec<f64> = Vec::with_capacity(n); 
//         // let mut c : Vec<f64> = Vec::with_capacity(n); 
//         // let mut d : Vec<f64> = Vec::with_capacity(n);
//         // for i in 0..n {
//         //     a1.push(pure[i].a1); 
//         //     a2.push(pure[i].a2); 
//         //     b1.push(pure[i].b1); 
//         //     b2.push(pure[i].b2); 
//         //     c .push(pure[i].c ); 
//         //     d .push(pure[i].d ); 
//         // }
//         // Self { a1, a2, b1, b2, c, d }
//         Self (pure.into_iter().next().unwrap())
        

//     }
// }
// pub struct ScalingParameters(ScalingRecord);
    // pub a1:Vec<f64>, 
    // pub a2:Vec<f64>, 
    // pub b1:Vec<f64>, 
    // pub b2:Vec<f64>, 
    // pub c: Vec<f64>, 
    // pub d: Vec<f64>,



pub struct Scaling<R:Residual> {
    pub parameters: ScalingRecord,
    pub entropy: f64,
    pub reference:f64,
    pub state: Arc<State<R>>,
    
}


impl<R: Residual> Scaling<R> {


    pub fn new(state: &Arc<State<R>>, parameters: ScalingRecord) -> Self {
        
        let entropy = state.entropy_isov() / IDEAL_GAS_CONST;
        let reference = super::rosenfeld_viscosity(state.t, state.d, state.molar_weight()[0]);

        Self { state: Arc::clone(state), entropy, reference, parameters }
    
    }


    pub fn xscaling(&self, sc:f64) -> f64 {
        
        let s = self.entropy;
        let r = s / sc;
        - r  - r.ln()
        
    }


    pub fn yscaling(&self, sc:f64) -> f64{
        
        let x = self.xscaling(sc);
        let s = self.entropy;
        let parameters = &self.parameters;
        let lhs = (parameters[0] + parameters[1] * s) / (1. + (parameters[4] * x).exp()) 
                    + (parameters[2] + parameters[3] * s) / (1. + (- parameters[4] * x).exp());

        let rhs = parameters[5] / sc;

        lhs * x + rhs
    }
    
    pub fn viscosity(&self, sc:f64) -> f64 {

        let mw = self.state.molar_weight()[0];
        let y = self.yscaling(sc);
        let reff = super::rosenfeld_viscosity(self.state.t, self.state.d, mw);
        
        reff * y.exp()

    }


}


impl<R:Residual> std::fmt::Display for Scaling<R> {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "ScalingDehlouz(\n\t\tstate={},\n\t\tη_ref={} Pa.s,\n\t\tS/R={},\n\t\tparameters={:?}",self.state,self.reference,self.entropy,self.parameters)
    }
}
// pub fn viscosity(residual: Residua)


