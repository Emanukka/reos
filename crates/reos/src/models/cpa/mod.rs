pub mod parameters;
pub mod rdf;
pub mod association;
#[cfg(test)]
pub mod bench;

use crate::models::cpa::association::AssociativeCPA;
use crate::models::cpa::parameters::{CPABinaryRecord, CPAParameters, CPAPureRecord};
use crate::models::cpa::rdf::{CS, Kontogeorgis, RdfModel};
use crate::models::cubic::{Cubic, CubicModel, SRK};


use crate::parameters::Parameters;
use crate::parameters::records::{BinaryRecord, PureRecord};
use crate::residual::Residual;
use crate::state::State;
use core::{f64, str};
use ndarray::{Array1, Array2};



// #[derive(Clone)]
pub struct CPA <C:CubicModel,R:RdfModel>{
    pub cubic: Cubic<C>,
    pub assoc: AssociativeCPA<R>,
}

impl<C:CubicModel,R:RdfModel> CPA<C,R> {

    pub fn from_parameters(parameters:CPAParameters<C>)->Self{
        let cubic=Cubic::from_parameters(parameters.cubic);

        let b_components = &cubic.parameters.b;

        let rdf = R::new(b_components.flatten().to_owned());
        let assoc = AssociativeCPA::from_parameters(parameters.assoc,rdf);

        Self
        {
            cubic,
            assoc,
        }
    }

    pub fn from_records(pure_records:Vec<PureRecord<CPAPureRecord>>,binary_records:Vec<BinaryRecord<CPABinaryRecord>>)->Self{

        let parameters = CPAParameters::new(pure_records, binary_records);
        Self::from_parameters(parameters)

    }
    // let l_t = &parameters.assoc.lambda_t;
    // let g_t = &parameters.assoc.gamma_t;
    // let b_sites = l_t.dot(&g_t.dot(b_components));
    // let bij = 0.5 * (&b_sites + &b_sites.t());


}

impl<C:CubicModel,R:RdfModel> Residual for CPA<C,R> {
    
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


impl<C:CubicModel,R:RdfModel> State<CPA<C,R>> {
    
    pub fn non_bonded_sites(&self)->Array1<f64> {

        // let (t,rho,x)=(self.t,self.d,&self.x);
        // let kmat = 
        // self.eos.residual.assoc.assoc.x_tan(x, kmat).unwrap_or_default()

        unimplemented!()
    }

}


pub type SCPA = CPA<SRK,Kontogeorgis>;


