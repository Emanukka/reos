pub mod parameters;
pub mod rdf;
pub mod association;
#[cfg(test)]
pub mod bench;

use crate::models::cpa::association::AssociativeCPA;
use crate::models::cpa::parameters::{CPABinaryRecord, CPAParameters, CPAPureRecord};
use crate::models::cpa::rdf::{Kontogeorgis, RdfModel};
use crate::models::cubic::{Cubic};

use crate::parameters::Parameters;
use crate::parameters::records::{BinaryRecord, PureRecord};
use crate::residual::Residual;
use core::f64;
use ndarray::Array1;



pub struct CPA <R:RdfModel>{
    pub cubic: Cubic,
    pub assoc: AssociativeCPA<R>,
}

type Pure = CPAPureRecord;
type Binary = CPABinaryRecord;

impl<R:RdfModel> CPA<R> {

    pub fn from_parameters(parameters:CPAParameters)->Self{
        let cubic= Cubic::from_parameters(parameters.cubic);

        let b_components = &cubic.parameters.b;

        let rdf = R::new(b_components.flatten().to_owned());
        let assoc = AssociativeCPA::from_parameters(parameters.assoc,rdf);

        Self
        {
            cubic,
            assoc,
        }
    }

    pub fn from_records(pure_records: Vec<PureRecord<Pure>>, binary_records: Vec<BinaryRecord<Binary>>, cubic_model: crate::models::cubic::models::CubicModels)->Self{
        
        let parameters = CPAParameters::new(pure_records, binary_records, cubic_model);
        Self::from_parameters(parameters)
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

#[cfg(test)]
mod tests {

    use approx::assert_relative_eq;
    use ndarray::array;
    use crate::{models::cpa::{CPA, parameters::CPAParameters, rdf::Kontogeorgis}, parameters::Parameters, residual::Residual};
    use super::parameters::readyto::water4c;
    use crate::models::cubic::models::{SRK};
    // use super::SCPA;

    
    fn water()->CPA<Kontogeorgis>{
        let water = water4c();
        let p = CPAParameters::new(vec![water], vec![], SRK.into());
        let model = CPA::from_parameters(p);
        model

    }
    
    #[test]
    fn test_scpa_helmholtz() {
        
        let model = water();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];
        let val = model.r_helmholtz(t, d, &x);

        assert_relative_eq!(val, -0.058144295861 + -1.597921023379 , epsilon = 1e-9)

    }

    #[test]
    fn test_scpa_entropy() {
        
        let model = water();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = model.r_entropy(t, d, &x);

        assert_relative_eq!(val, -0.041951593945 + -4.713659269705, epsilon = 1e-9)

    }

    #[test]
    fn test_scpa_chem_pot() {
        
        let model = water();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = model.r_chemical_potential(t, d, &x);

        assert_relative_eq!(val[0], -0.115660251059 + -2.54386196979185 , epsilon = 1e-10)


    }
    
    #[test]
    fn test_scpa_pressure() {
        
        let model = water();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = model.r_pressure(t, d, &x);

        assert_relative_eq!(val, -57.5159551979349 + -945.9409464127781, epsilon = 1e-8)

    }
}