use std::marker::PhantomData;

use ndarray::{Array1, Array2};
use crate::models::associative::Associative;
use crate::models::associative::parameters::{ AssociativeParameters};
use crate::models::cpa::rdf::{Rdf, RdfModel};
use crate::residual::Residual;
use crate::state::S;
use crate::state::eos::{EosError, EosResult};


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

    fn components(&self)->usize {
        // self.assoc.parameters.gamma.nrows()
        todo!()
    }
    fn bmix(&self,_x:&Array1<f64>)->f64 {
        todo!()
    }
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let kmat = self.assoc.association_constants(t, rho, x, volf);
        let unbonded= self.assoc.x_tan(x,&kmat)?;
        let dlng_drho = self.rdf.dlngdrho(rho, x);

        Ok(
        self.assoc.r_pressure(rho, dlng_drho, x, &unbonded)
        )
    }

    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let kmat = self.assoc.association_constants(t, rho, x, volf);
        let unbonded= &self.assoc.x_tan(x,&kmat)?;

        let ndlng_dni = &self.rdf.ndlngdni(rho, x);
        Ok(
        self.assoc.r_chemical_potential(x, ndlng_dni, unbonded)
        )
    }
    
    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        let volf = &self.rdf.bij * self.rdf.g(rho, x);
        let kmat = self.assoc.association_constants(t, rho, x, volf);
        let unbonded= &self.assoc.x_tan(x,&kmat)?;

        Ok(
        self.assoc.r_helmholtz(x, unbonded)
        )
    }
    
}


