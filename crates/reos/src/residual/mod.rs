use approx::assert_relative_eq;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

use crate::state::eos::EosResult;

#[derive(Default,Serialize,Deserialize)]
pub struct ResidualDerivedProperties{
    pub a:f64,
    pub dadv:f64,
    pub dadni:Vec<f64>,
}

impl ResidualDerivedProperties {
    
    pub fn comparison(&self,other:ResidualDerivedProperties,tol:Option<f64>){

        let tol = tol.unwrap_or(1e-8);
        assert_relative_eq!(self.a,other.a,epsilon = tol);
        assert_relative_eq!(self.dadv,other.dadv,epsilon = tol);
        
        let n = self.dadni.len();
        for i in 0..n{
            assert_relative_eq!(self.dadni[i],other.dadni[i],epsilon = tol);
        }
    }
}

///A trait that provides a interface for all types 
///that have a ```non-ideal``` contribution to the thermodynamic state. 
///All properties are in the reduced form (less the pressure), and they
/// are at the natural coordinates ```NVT```  
/// Example:
/// ```text
///    s = S/R
///    a = A/R/T
///    h = H/R/T
///    P = d(A/RT)/dV
/// ```
/// 
/// 
pub trait Residual{

    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>;
    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64>;
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64>;
    fn bmix(&self,x:&Array1<f64>)->f64;

    fn components(&self)->usize;

    fn all_derived_properties(&self,t:f64,rho:f64,x:&Array1<f64>)->ResidualDerivedProperties{

        let mut r_properties = ResidualDerivedProperties::default();
        r_properties.a = self.r_helmholtz(t, rho, x).unwrap();
        r_properties.dadv = self.r_pressure(t, rho, x).unwrap();
        r_properties.dadni = self.r_chemical_potential(t, rho, x).unwrap().to_vec();

        r_properties
        // arr

    }
}








