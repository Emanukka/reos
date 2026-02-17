pub mod parameters;
pub mod rdf;
pub mod association;

use crate::models::cpa::association::AssociativeCPA;
use crate::models::cpa::parameters::{CPABinaryRecord, CPAParameters, CPAPureRecord};
use crate::models::cpa::rdf::{Kontogeorgis, RdfModel};
use crate::models::cubic::options::CubicOptions;
use crate::models::cubic::{Cubic};

use crate::parameters::Parameters;
use crate::parameters::records::{BinaryRecord, PureRecord};
use crate::residual::Residual;
use core::f64;
use ndarray::Array1;
#[cfg(test)]
mod tests;


pub struct CPA <R:RdfModel>{
    pub cubic: Cubic,
    pub assoc: AssociativeCPA<R>,
}

pub type Pure = PureRecord<CPAPureRecord>;
pub type Binary = BinaryRecord<CPABinaryRecord>;
// type Options = CubicOptions;

impl<R:RdfModel> CPA<R> {

    pub fn from_parameters(parameters:CPAParameters)->Self{
        let cubic= Cubic::from_parameters(parameters.cubic);

        let b_components = &cubic.parameters.b;

        let rdf = R::new(Array1::from_vec(b_components.clone()));
        let assoc = AssociativeCPA::from_parameters(parameters.assoc,rdf);

        Self
        {
            cubic,
            assoc,
        }
    }



}

impl<R:RdfModel> Residual for CPA<R> {
    
    fn get_properties(&self)->&crate::parameters::Properties {
        &self.cubic.parameters.properties
    }
    fn molar_weight(&self)->&Array1<f64> {
        self.cubic.molar_weight()
    }
    fn components(&self)->usize {
        self.cubic.components()
    }
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->f64 {

        self.cubic.r_pressure(t, rho, x) + self.assoc.r_pressure(t, rho, x)

    }
    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->Array1<f64> {

        self.cubic.r_chemical_potential(t, rho, x) + self.assoc.r_chemical_potential(t, rho, x)
    }
    fn max_density(&self,x:&Array1<f64>)->f64 {
        self.cubic.max_density(x)
    }

    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->f64{

        self.cubic.r_helmholtz(t, rho, x) + self.assoc.r_helmholtz(t, rho, x)
        
    }

    fn r_entropy(&self,t:f64,rho:f64,x:&Array1<f64>)->f64{

        self.cubic.r_entropy(t, rho, x) + self.assoc.r_entropy(t, rho, x)
        
    }
}



impl<R: RdfModel> CPA<R> {
    
    pub fn unbonded_sites(&self,t:f64, d:f64, x:&Array1<f64>) -> Array1<f64> {
        // self.model.assoc.unbonded_sites(&self.state)
        let assoc = &self.assoc;
        let volf = &assoc.rdf.bij * assoc.rdf.g(d, x);
        let k = assoc.assoc.association_constants(t, d, x, &volf);
        let u = assoc.assoc.unbonded_sites_fraction(x, &k);
        u
    }

    pub fn association_constants(&self,t:f64, d:f64, x:&Array1<f64>) -> ndarray::Array2<f64> {
        let assoc = &self.assoc;
        let volf = &assoc.rdf.bij * assoc.rdf.g(d, x);
        let k = assoc.assoc.association_constants(t, d, x, &volf);
        k
    }
}
pub type SCPA = CPA<Kontogeorgis>;
