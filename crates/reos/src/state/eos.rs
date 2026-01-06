
use core::f64;

use ndarray::Array1;
// use ndarray_linalg::error::LinalgError;
use thiserror::Error;

use crate::{models::{ IDEAL_GAS_CONST as R}, residual::Residual};

pub type EosResult<T> = Result<T, EosError>;



pub struct EquationOfState<R>{
    pub residual:R,
    // pub properties:Option<PureProperties>
}

impl<R:Residual> EquationOfState<R> {
    
    pub fn from_residual(r:R)->Self{
        Self{
            residual:r,

        }
    }
}
impl <R:Residual> EquationOfState<R> {
    
    pub fn ideal_gas_pressure(&self,t: f64,d: f64)->f64{
        d * R * t
    }

    pub fn pressure(&self,t: f64,d: f64,x: &Array1<f64>)->f64{

        let r_pres = self.residual.r_pressure(t, d, x);
        let r_pig = d;
        R * t * ( r_pres + r_pig)

    }

    pub fn helmholtz(&self,t: f64,d: f64,x: &Array1<f64>)->f64 {

        R * t * self.residual.r_helmholtz(t, d, x)
    }

    pub fn compressibility(&self,t:f64,d:f64,x:&Array1<f64>)->f64{

        self.pressure(t, d, x) / self.ideal_gas_pressure(t, d)

    }

    pub fn lnphi(&self,t:f64,d:f64,x:&Array1<f64>)->Array1<f64>{

        self.residual.r_chemical_potential(t, d, x) - self.compressibility(t, d, x).ln()
        
    }
}

#[derive(Error,Debug)]
pub enum EosError {
    #[error("`{0}` Not Converged.")]
    NotConverged(String),
    #[error("Phase must be 'liquid' or 'vapor'.")]
    PhaseError,


}