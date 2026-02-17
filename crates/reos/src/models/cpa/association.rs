use ndarray::Array1;
use crate::models::associative::Associative;
use crate::models::associative::parameters::{ AssociativeParameters};
use crate::models::cpa::rdf::{Rdf, RdfModel};
use crate::residual::Residual;


#[derive(Clone)]
pub struct AssociativeCPA<R:RdfModel>{
    pub assoc:Associative,
    pub rdf:Rdf<R>,

}


impl<R:RdfModel> AssociativeCPA<R> {
    
    pub fn from_parameters(
        parameters:AssociativeParameters,
        rdf:Rdf<R>
        )->Self{
        Self
        {
            assoc: Associative::from_parameters(parameters),
            rdf,
        }
    }
}


impl<R:RdfModel> Residual for AssociativeCPA<R> {
    
    fn get_properties(&self)->&crate::parameters::Properties {
        panic!()
    }
    fn molar_weight(&self)->&Array1<f64> {
        panic!()
    }
    fn components(&self)->usize {
        // self.assoc.parameters.gamma.nrows()
        panic!()
    }
    fn max_density(&self,_x:&Array1<f64>)->f64 {
        panic!()
    }
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>) -> f64 {

        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let unbonded = self.assoc.unbonded_sites_fraction(x, k);
        let h = self.assoc.h(x, &unbonded);

        let dlng_drho = self.rdf.dlngdrho(rho, x);

        self.assoc.r_pressure(h, rho, dlng_drho)
    }

    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->Array1<f64>{

        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let unbonded = &self.assoc.unbonded_sites_fraction(x, k);
        let h = self.assoc.h(x, &unbonded);

        let ndlng_dni = &self.rdf.ndlngdni(rho, x);

        self.assoc.r_chemical_potential(h, ndlng_dni, unbonded)
    }
    
    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>) -> f64 {

        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let unbonded = &self.assoc.unbonded_sites_fraction(x, k);

        self.assoc.r_helmholtz(x, unbonded)
    }

    fn r_entropy(&self,t:f64, rho:f64, x:&Array1<f64>)->f64 {

        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let k = &self.assoc.association_constants(t, rho, x, &volf);

        self.assoc.r_entropy(t, x ,k)
        
    }
    
}

// #[cfg(test)]
// mod tests{
//     use std::sync::Arc;

//     use approx::assert_relative_eq;
//     use ndarray::{Array1, array};

//     use crate::{arr_eq, models::{cpa::{SCPA, parameters::readyto::{from_records, acetic1a, acoh_octane, co2, methanol3b, octane, water4c, water4c_acetic1a, water4c_co2}}, cubic::models::SRK}, state::{E, S, density_solver::DensityInitialization}};


// }