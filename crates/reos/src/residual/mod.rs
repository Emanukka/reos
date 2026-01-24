use approx::assert_relative_eq;
use ndarray::Array1;
use serde::{Deserialize, Serialize};

use crate::parameters::Properties;


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

/// API for computing dimensionless residual isovolumetric thermodynamic properties
pub trait Residual{

    
    fn get_properties(&self)->&Properties;
    
    fn molar_weight(&self)->&Array1<f64>;
    
    fn r_entropy(&self,t:f64, d:f64, x:&Array1<f64>)->f64;

    fn r_chemical_potential(&self,t:f64,d: f64,x:&Array1<f64>)->Array1<f64>;

    fn r_helmholtz(&self,t:f64,d: f64,x:&Array1<f64>)->f64;

    fn r_pressure(&self,t:f64,d: f64,x:&Array1<f64>)->f64;

    fn compressibility(&self,t:f64,d:f64,x:&Array1<f64>)->f64{
        let r_pres = self.r_pressure(t, d, x);
        r_pres / d
    }


    fn max_density(&self,x:&Array1<f64>)->f64;

    fn components(&self)->usize;

    fn all_derived_properties(&self,t:f64,d: f64,x:&Array1<f64>)->ResidualDerivedProperties{

        let mut r_properties = ResidualDerivedProperties::default();
        r_properties.a = self.r_helmholtz(t, d, x);
        r_properties.dadv = self.r_pressure(t, d, x);
        r_properties.dadni = self.r_chemical_potential(t, d, x).to_vec();

        r_properties

    }
}


