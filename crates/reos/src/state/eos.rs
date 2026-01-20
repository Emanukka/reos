
use core::f64;

use ndarray::Array1;
// use ndarray_linalg::error::LinalgError;
use thiserror::Error;

use crate::{models::IDEAL_GAS_CONST as R, parameters::Properties, residual::Residual};

pub type EosResult<T> = Result<T, EosError>;


/// Equation of State API
pub struct EquationOfState<R>{
    residual:R,
    // pub properties:Option<PureProperties>
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
    
    /// Ideal gas pressure in Pa
    pub fn ideal_gas_pressure(&self,t: f64,d: f64)->f64{
        d * R * t
    }

    /// Pressure in Pa
    pub fn pressure(&self,t: f64,d: f64,x: &Array1<f64>)->f64{

        let r_pres = self.residual.r_pressure(t, d, x);
        let r_pig = d;
        R * t * ( r_pres + r_pig)

    }
    /// Residual Isovolumetric Helmholtz free energy in J / mol
    pub fn helmholtz_isov(&self,t: f64,d: f64,x: &Array1<f64>)->f64 {

        R * t * self.residual.r_helmholtz(t, d, x)
    }

    /// Residual Isovolumetric Entropy in J / mol / K
    pub fn entropy_isov(&self,t: f64,d: f64,x: &Array1<f64>)->f64 {

        R * self.residual.r_entropy(t, d, x)
    }

    /// Residual Entropy in J / mol / K
    pub fn entropy(&self,t: f64,d: f64,x: &Array1<f64>)->f64 {

        let isov = self.residual.r_entropy(t, d, x);
        let z = 1.0 + self.residual.compressibility(t, d, x);
        let s = isov + z.ln();
        
        R * s
        
    }
    /// Compressibility factor
    pub fn compressibility(&self,t:f64,d:f64,x:&Array1<f64>)->f64{

        self.pressure(t, d, x) / self.ideal_gas_pressure(t, d)

    }

    /// Natural logarithm of the fugacity coefficient
    pub fn lnphi(&self,t:f64,d:f64,x:&Array1<f64>)->Array1<f64>{

        self.residual.r_chemical_potential(t, d, x) - self.compressibility(t, d, x).ln()
        
    }

    /// Residual Isovolumetric Chemical potential in J / mol
    pub fn chem_pot_isov(&self,t:f64, d:f64, x:&Array1<f64>)->Array1<f64> {
        
        R * t * self.residual.r_chemical_potential(t, d, x) 

    }

    /// Residual Chemical potential in J / mol
    pub fn chem_pot(&self,t:f64, d:f64, x:&Array1<f64>)->Array1<f64> {
        
        let lnphi = self.lnphi(t, d, x);
        R * t * lnphi 

    }

    /// Residual Gibbs energy in J / mol
    pub fn gibbs(&self, t: f64, d: f64, x: &Array1<f64>) -> f64 {

        let a_res_isov = self.residual.r_helmholtz(t, d, x);
        let z_res = self.residual.compressibility(t, d, x);
        let z = 1.0 + z_res;

        let g = a_res_isov + z_res - z.ln();

        R * t * g
        // R * t * self.residual.r_gibbs(t, d, x)
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
    #[error("Phase must be 'liquid' or 'vapor'.")]
    PhaseError,


}