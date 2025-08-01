use ndarray::Array1;

use crate::{models::{associative::Associative, cpa::CPA}, state::eos::EosResult};

pub trait Residual{

    ///Adimensional
    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>;
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64>;
    fn bmix(&self,x:&Array1<f64>)->f64;
    fn components(&self)->usize;
}








