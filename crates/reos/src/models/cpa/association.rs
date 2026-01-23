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

#[cfg(test)]
mod tests{
    use std::sync::Arc;

    use approx::assert_relative_eq;
    use ndarray::{Array1, array};

    use crate::{arr_eq, models::{cpa::{SCPA, parameters::readyto::{acetic1a, acoh_octane, co2, methanol3b, octane, water4c, water4c_acetic1a, water4c_co2}}, cubic::models::SRK}, state::{E, S, density_solver::DensityInitialization}};


    #[test]
    fn test_assoc_consts_4c(){

        let t=298.15;
        let rho= 1000.0;

        let x= &array![1.0];
        let comp1 = water4c();
        let assoc = SCPA::from_records(vec![comp1],vec![], SRK.into()).assoc;
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;

        let k = assoc.assoc.association_constants(t, rho, x,&volf);
        let reff = array![[0.           , 0.83518048731],[0.83518048731, 0.           ]];
        
        k.iter().zip(reff.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });

    }

    #[test]
    fn test_assoc_consts_3b(){

        let t=298.15;
        let rho= 1000.0;

        let x= &array![1.0];
        let comp1 = methanol3b();
        let assoc = SCPA::from_records(vec![comp1],vec![], SRK.into()).assoc;
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;

        let k = assoc.assoc.association_constants(t, rho, x,&volf);
        let reff = array![[0.           , 0.76195179893041],
                                                     [0.76195179893041, 0.           ]];  

        k.iter().zip(reff.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });

    }


    #[test]
    fn test_assoc_consts_1a_inert(){
        let t=298.15;
        let rho= 7_500.;
        let x= &array![0.25,0.75];

        let comp1 = acetic1a();
        let comp2 = octane();
        let b = acoh_octane();
        let assoc = SCPA::from_records(vec![comp1,comp2],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k_rust = assoc.assoc.association_constants(t, rho, x,&volf);

        let ok = arr_eq!(k_rust, array![ 1980.9483616 ], 1e-7);

        assert!(ok)
    }
    #[test]
    fn unbonded_acetic_octane(){
        let t=298.15;
        let rho= 7_500.;
        let x= &array![0.25,0.75];

        let comp1 = acetic1a();
        let comp2 = octane();
        let b = acoh_octane();
        let assoc = SCPA::from_records(vec![comp1,comp2],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = assoc.assoc.association_constants(t, rho, x,&volf);

        let m = &assoc.assoc.sites_mole_frac(x);
        let tan = assoc.assoc.x_tan(m, &k).unwrap();
        // let michelsen = assoc.assoc.x_michelsen(m, &k).unwrap();
        let analytic = assoc.assoc.unbonded_sites_fraction(x,&k);
        
        let ok = arr_eq!(tan, analytic, 1e-10);


        assert!(ok)
    }

    #[test]
    fn phi_acetic_octane(){
        let t=298.15;
        let rho= 7_500.;
        let x= &array![0.25,0.75];

        let comp1 = acetic1a();
        let comp2 = octane();
        let b = acoh_octane();
        // let assoc = SCPA::from_records(vec![comp1,comp2],vec![b]).assoc;
        
        let cpa = SCPA::from_records(vec![comp1,comp2],vec![b], SRK.into());
        
        let s = S::new_trx(Arc::new(E::from_residual(cpa)), t, rho, x.clone());
        let phi_rust = s.lnphi().exp();
        
        let ok = arr_eq!(phi_rust,array![0.0002995 , 0.00071522],1e-6);

        assert!(ok)
    }
    #[test]
    fn association_constants_water_co2() {

        let t=298.15;
        let rho= 55_000.0;
        let x= &array![0.5,0.5];

        let comp1 = water4c();
        let comp2 = co2();
        let b = water4c_co2();
        let assoc = SCPA::from_records(vec![comp1,comp2],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k_rust = assoc.assoc.association_constants(t, rho, x,&volf);

        let ok = arr_eq!(k_rust, array![ 0. , 25.04896564,  3.21025658, 25.04896564,  0., 0., 3.21025658,  0.,  0.], 1e-7);

        assert!(ok)
    }

    #[test]
    fn unbonded_water_co2() {
        let t= 298.15;
        let rho= 55_000.0;
        let x= &Array1::from_vec(vec![0.5,0.5]);
        
        let water = water4c();
        let co2 = co2();
        let b = water4c_co2();
        let assoc = SCPA::from_records(vec![water,co2],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        
        let kmat = assoc.assoc.association_constants(t, rho, x, &volf);
        let m = &assoc.assoc.sites_mole_frac(x);
        let unbonded = assoc.assoc.x_tan(m, &kmat).unwrap();
        // let unbonded = assoc.assoc.x_michelsen(m, &kmat).unwrap();

        // let ok = arr_eq!(unbonded, array![0.03874146, 0.20484502, 0.66778819],1e-5);
        // let ok = arr_eq!(unbonded, array![0.03874146, 0.20484502, 0.66778819],1e-5);

        dbg!(unbonded);
        // assert!(ok);
    }

    #[test]
    fn association_constants_water_acetic() {

        let t=298.15;
        let rho= 21_500.0;
        let x= &array![0.25,0.75];

        let comp1 = water4c();
        let comp2 = acetic1a();
        let b = water4c_acetic1a();
        let assoc = SCPA::from_records(vec![comp1,comp2],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k_rust = assoc.assoc.association_constants(t, rho, x, &volf);

        let reff =  array![0.00000000e+00, 1.84368158e+00, 3.00115845e+02,
                                              1.84368158e+00, 0.00000000e+00, 3.00115845e+02, 
                                              3.00115845e+02, 3.00115845e+02, 4.88530782e+04] ;

        let ok = arr_eq!(k_rust,reff,1e-6);
        assert!(ok)
    }

    #[test]
    fn unbonded_water_acetic() {
        let t=298.15;
        let rho= 21_500.0;
        let x= &array![0.25,0.75];

        let comp1 = water4c();
        let comp2 = acetic1a();
        let b = water4c_acetic1a();
        let assoc = SCPA::from_records(vec![comp1,comp2],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = assoc.assoc.association_constants(t, rho, x, &volf);

        let unb = assoc.assoc.unbonded_sites_fraction(x,&k);

        let reff = array![0.1598388 , 0.1598388 , 0.00241471];

        let ok = arr_eq!(unb,reff,1e-6);

        assert!(ok)
    }

    #[test]
    fn permutation_water_acetic() {

        let t=298.15;
        let rho= 21_500.0;
        let x= &array![0.75, 0.25];

        let comp1 = water4c();
        let comp2 = acetic1a();
        let b = water4c_acetic1a();
        let assoc = SCPA::from_records(vec![comp2,comp1],vec![b], SRK.into()).assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;

        let k = assoc.assoc.association_constants(t, rho, x, &volf);
        let unb = assoc.assoc.unbonded_sites_fraction(x,&k);

        let reff = array![0.00241471, 0.1598388 , 0.1598388 ];

        let ok = arr_eq!(unb, reff, 1e-6);

        assert!(ok)
    }
    #[test]
    fn phi_water_co2(){

        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);
        let water = water4c();
        let co2 = co2();
        let b = water4c_co2();
        let cpa = SCPA::from_records(vec![water,co2],vec![b], SRK.into());
        
        let s = S::new_tpx(Arc::new(E::from_residual(cpa)), t, p, x.clone(), Some(DensityInitialization::Vapor)).unwrap();
        let phi_rust = s.lnphi().exp();
        
        let ok = arr_eq!(phi_rust,array![2.14385745e-04, 5.65853284e-01],1e-6);

        assert!(ok)


    }


}