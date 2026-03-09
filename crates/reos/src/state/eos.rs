
use core::f64;

use ndarray::Array1;
// use ndarray_linalg::error::LinalgError;
use thiserror::Error;

use crate::{models::R_GAS, parameters::Properties, residual::Residual};

pub type EosResult<T> = Result<T, EosError>;


/// Equation of State API
pub struct EquationOfState<R>{
    residual:R,
}

impl<R:Residual> EquationOfState<R> {
    
    pub fn from_residual(r:R)->Self{
        Self{
            residual:r,

        }
    }

    pub fn residual(&self)->&R{
        &self.residual
    }

}

impl<T: Residual> From<T> for EquationOfState<T> {
    
    fn from(residual: T) -> Self {
        Self { residual }
    }
}



impl <R:Residual> EquationOfState<R> {
    
    /// Ideal gas pressure in `[Pa]`
    pub fn ideal_gas_pressure(&self,t: f64,d: f64)->f64{
        d * R_GAS * t
    }

    /// Compressibility factor, defined as:
    /// 
    /// `
    /// Z(T,ρ,x) = Zⁱᵍ + Zʳ = 1 - V∂F/∂V
    /// `
    pub fn compressibility(&self,t:f64, d:f64, x:&Array1<f64>)->f64{

        let df_dv= self.residual.df_dv(t, d, x);
        let z_res =  - df_dv / d;
        1.0 + z_res 

    }

    /// Pressure `[Pa]`, defined as:
    /// 
    /// `
    /// P(T,ρ,x) = ZρRT
    /// `
    pub fn pressure(&self,t: f64, d: f64,x: &Array1<f64>)->f64{

        let z = self.compressibility(t, d, x);
        z * d * R_GAS * t


    }
    
    /// Residual Helmholtz free energy `[J / mol]`, defined as:
    /// 
    /// `
    /// Aʳ(T,ρ,x) = RT ⋅ F
    /// ` 
    pub fn helmholtz(&self,t: f64,d: f64,x: &Array1<f64>)->f64 {
        R_GAS * t * self.residual.helmholtz(t, d, x)
    }

    /// Residual TV Entropy `[J / mol / K]`, defined as:
    /// 
    /// `
    /// Sʳ(T,ρ,x) = R(- F - T∂F/∂T)
    /// `
    pub fn tv_entropy(&self,t: f64, d: f64,x: &Array1<f64>)->f64 {

        let df_dt = self.residual.df_dt(t, d, x);
        let f = self.residual.helmholtz(t, d, x);
        let s = - f - t * df_dt;

        R_GAS * s
    }

    /// Residual TP Entropy `[J / mol / K]`, defined as:
    /// 
    /// `
    /// Sʳ(T,P,x) =  Sʳ(T,ρ,x) + Rln(Z)
    /// `
    pub fn tp_entropy(&self,t: f64,d: f64, x: &Array1<f64>)->f64 {

        // let s_isov = self.residual.r_entropy(t, d, x);
        let z = self.compressibility(t, d, x);
        let s_tv = self.tv_entropy(t, d, x);
        s_tv + R_GAS * z.ln()
        
        
    }


    /// Natural logarithm of the fugacity coefficient, defined as:
    /// 
    /// `
    /// ln(ϕᵢ) = ∂F/∂nᵢ - ln(Z)
    /// `
    pub fn lnphi(&self, t:f64, d:f64, x:&Array1<f64>)->Array1<f64>{

        self.residual.df_dn(t, d, x) - self.compressibility(t, d, x).ln()
        
    }

    /// Residual TP Chemical potential in `[J / mol]`, defined as:
    /// 
    /// `
    /// μᵢʳ(T,P,x) = RTln(ϕᵢ)
    /// `
    pub fn chem_pot(&self,t:f64, d:f64, x:&Array1<f64>)->Array1<f64> {
        
        let lnphi = self.lnphi(t, d, x);
        R_GAS * t * lnphi 

    }

    /// Residual Gibbs energy in `[J / mol]`, defined as:
    /// 
    /// `
    /// Gʳ(T,P,x) = Aʳ + RTZʳ - RTln(Z)
    /// `
    pub fn gibbs(&self, t: f64, d: f64, x: &Array1<f64>) -> f64 {

        let f = self.residual.helmholtz(t, d, x);
        let z = self.compressibility(t, d, x);
        let g = f + z - 1.0 - z.ln();
        R_GAS * t * g

    }

    pub fn max_density(&self, x:&Array1<f64>)->f64{

        self.residual.max_density(x)
    }
    pub fn molar_weight(&self)->&Array1<f64>{

        self.residual.molar_weight()

    }

    pub fn get_properties(&self)-> &Properties{

        self.residual.get_properties()
        
    }
}

#[derive(Error,Debug)]
pub enum EosError {
    #[error("{0}")]
    NotConverged(String),
    #[error("Phase must be 'liquid', 'vapor' or 'stable'")]
    PhaseError,


}