
use ndarray::Array1;
use thiserror::Error;

use crate::{models::IDEAL_GAS_CONST, residual::Residual};


// Define: 
// 1) trait Residual: 
// 
// 2) struct EquationOfState: (Equação de estado genérica
// de uma Estado Termodinâmico, possibilitando assim cálculos de propriedades)
// (genérica pois possui um objeto que necessariamente implemente trait/interface 'Residual';


pub type EosResult<T> = Result<T, EosError>;


pub struct EquationOfState<R>{
    pub residual:R,
}

impl<R:Residual> EquationOfState<R> {
    
    pub fn from_residual(r:R)->Self{
        Self{residual:r}
    }
}
impl <R:Residual> EquationOfState<R> {
    
    fn ideal_gas_pressure(t: f64,rho: f64)->f64{
        rho*IDEAL_GAS_CONST*t
    }
    pub fn pressure(&self,t: f64,rho: f64,x: &ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>)->EosResult<f64>{
        Ok(self.residual.pressure(t, rho, x)? + Self::ideal_gas_pressure(t, rho))
    }
    pub fn compressibility(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64>{

        // P_residual + P_gas_ideal
        // dbg!(self.pressure(t, rho, x)?);
        Ok(self.pressure(t, rho, x)?/(rho*IDEAL_GAS_CONST*t))
    }

    /// Potencial Químico Adimensional - lnZ
    pub fn lnphi(&self,t:f64,rho:f64,x:&Array1<f64>)->Result<Array1<f64>,EosError>{
        Ok(
            self.residual.residual_chemical_potential(t, rho, x)? - self.compressibility(t, rho, x)?.ln()
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