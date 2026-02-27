use ndarray::Array1;
use crate::models::associative::Associative;
use crate::models::associative::parameters::{ AssociativeParameters};
use crate::models::cpa::rdf::{RDFcpa};
use crate::residual::Residual;


#[derive(Clone)]
pub struct AssociativeCPA{
    pub assoc:Associative,
    pub rdf:RDFcpa,

}

impl AssociativeCPA {
    
    pub fn from_parameters(parameters:AssociativeParameters,rdf:RDFcpa) -> Self{
        Self{
            assoc: Associative{parameters},
            rdf,
        }
    }
}


impl Residual for AssociativeCPA {
    
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
    fn df_dv(&self,t:f64,rho:f64,x:&Array1<f64>) -> f64 {

        let volf = self.rdf.g(rho, x) * &self.rdf.b;
    // let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let unbonded = self.assoc.unbonded_sites_fraction(x, k);
        let h = self.assoc.h(x, &unbonded);

        let dlng_drho = self.rdf.dlngdrho(rho, x);

        self.assoc.df_dv(h, rho, dlng_drho)
    }

    fn df_dn(&self,t:f64,rho:f64,x:&Array1<f64>)->Array1<f64>{

        let volf = self.rdf.g(rho, x) * &self.rdf.b;
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let unbonded = &self.assoc.unbonded_sites_fraction(x, k);
        let h = self.assoc.h(x, &unbonded);

        let ndlng_dni = &self.rdf.ndlngdni(rho, x);

        self.assoc.df_dn(h, ndlng_dni, unbonded)
    }
    
    fn helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>) -> f64 {

        let volf = self.rdf.g(rho, x) * &self.rdf.b;
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let unbonded = &self.assoc.unbonded_sites_fraction(x, k);

        self.assoc.helmholtz(x, unbonded)
    }

    fn df_dt(&self,t:f64, rho:f64, x:&Array1<f64>)->f64 {

        let volf = self.rdf.g(rho, x) * &self.rdf.b;
        let k = &self.assoc.association_constants(t, rho, x, &volf);
        let dk_dt = &self.assoc.dk_dt(t, rho, x, &volf);

        // dbg!(dk_dt);        

        // dbg!(k / rho );
        self.assoc.df_dt(t, x , k, dk_dt)
        
    }
    
}

