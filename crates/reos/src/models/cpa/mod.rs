pub mod parameters;
pub mod rdf;
pub mod association;
pub mod tests;
use crate::models::cpa::association::AssociativeCPA;
use crate::models::cpa::parameters::CPAParameters;
use crate::models::cpa::rdf::{CS, Kontogeorgis, RdfModel};
use crate::models::cubic::{Cubic, CubicModel, SRK};


use crate::residual::Residual;
use crate::state::eos::EosResult;
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


        let bik = 0.5 * (b_components + &b_components.t());

        let rdf = R::new(b_components.flatten().to_owned(),bik);
        let assoc = AssociativeCPA::from_parameters(parameters.assoc,rdf);

        Self
        {
            cubic,
            assoc,
        }
    }

    // let l_t = &parameters.assoc.lambda_t;
    // let g_t = &parameters.assoc.gamma_t;
    // let b_sites = l_t.dot(&g_t.dot(b_components));
    // let bij = 0.5 * (&b_sites + &b_sites.t());


}

impl<C:CubicModel,R:RdfModel> Residual for CPA<C,R> {
    
    fn components(&self)->usize {
        self.cubic.components()
    }
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        Ok(
        self.cubic.r_pressure(t, rho, x)?
        +self.assoc.r_pressure(t, rho, x)?
        )

    }
    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        Ok(
            self.cubic.r_chemical_potential(t, rho, x)?
            +self.assoc.r_chemical_potential(t, rho, x)?
        )
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.cubic.bmix(x)
    }

    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        Ok(
        self.cubic.r_helmholtz(t, rho, x)?+self.assoc.r_helmholtz(t, rho, x)?)
        
    }
}


impl<C:CubicModel,R:RdfModel> State<CPA<C,R>> {
    
    pub fn non_bonded_sites(&self)->Array1<f64> {

        let (t,rho,x)=(self.t,self.rho,&self.x);
        // let kmat = 
        // self.eos.residual.assoc.assoc.x_tan(x, kmat).unwrap_or_default()

        todo!()
    }

}


pub type SCPA = CPA<SRK,Kontogeorgis>;


