
use ndarray::Array1;
use thiserror::Error;

use crate::{models::{ IDEAL_GAS_CONST}, residual::Residual};

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
    
    pub fn ideal_gas_pressure(&self,t: f64,rho: f64)->f64{
        rho*IDEAL_GAS_CONST*t
    }
    pub fn pressure(&self,t: f64,rho: f64,x: &ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>)->EosResult<f64>{
        let p_res_dimenssionless = self.residual.r_pressure(t, rho, x)?;
        let p_ig_dimenssionless = rho;
        // Ok(IDEAL_GAS_CONST*t*self.residual.r_pressure(t, rho, x)? + Self::ideal_gas_pressure(t, rho))
        Ok( IDEAL_GAS_CONST*t*( p_res_dimenssionless + p_ig_dimenssionless))
    }
    pub fn compressibility(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64>{

        Ok(self.pressure(t, rho, x)?/self.ideal_gas_pressure(t, rho))

    }

    pub fn lnphi(&self,t:f64,rho:f64,x:&Array1<f64>)->Result<Array1<f64>,EosError>{
        Ok(
            self.residual.r_chemical_potential(t, rho, x)? - self.compressibility(t, rho, x)?.ln()
        )
    }
}

#[derive(Error,Debug)]
pub enum EosError {
    #[error("`{0}` Not Converged.")]
    NotConverged(String),
    #[error("Phase must be 'liquid' or 'vapor'.")]
    PhaseError,
    #[error("Rule must be: 'cr1','ecr','mcr1' or 'exp'")]
    RuleError
}