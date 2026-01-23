pub mod eos;
pub mod density_solver;

use ndarray::{Array1, array};
use std::sync::Arc;

use crate::parameters::Properties;
use crate::residual::Residual;
use crate::state::eos::{EosError, EquationOfState};
use crate::state::density_solver::{DensityInitialization, density, liquid_solver, vapor_solver};
pub type E<R> = EquationOfState<R>;
pub type S<R> = State<R>;

pub type StateResult<R> = Result<S<R>, EosError>;



#[derive(Clone)]
pub struct State<R:Residual>{

    pub eos:Arc<EquationOfState<R>>,
    pub t:f64,
    pub p:f64,
    pub d:f64,
    pub x:Array1<f64>,
}



impl<R:Residual> State<R> {

    pub fn new_trx(eos:Arc<E<R>>, t:f64, d:f64, x:Array1<f64>)-> Self{

        let p = eos.pressure(t, d, &x);

        Self{
            eos:Arc::clone(&eos),
            t,
            p,
            d,
            x,
        }


    }
    pub fn new_tpx(eos:Arc<E<R>>, t:f64, p:f64, x:Array1<f64>, phase:Option<DensityInitialization>)->StateResult<R>{

        let phase = match phase {
            Some(ph) => ph,
            None => DensityInitialization::Stable
        };

        match phase {

            DensityInitialization::Liquid=>{

               liquid_solver(eos, t, p, x)

            }

            DensityInitialization::Vapor=>{

                vapor_solver(eos, t, p, x)

            }

            DensityInitialization::Guess(old_density)=>{

                let bm = 1./eos.max_density(&x);
                let guess = old_density * bm;
                density(eos, t, p, x, guess)

            }

            DensityInitialization::Stable => {  

                let vapor = vapor_solver(Arc::clone(&eos), t, p, x.clone());
                let liquid = liquid_solver(Arc::clone(&eos), t, p, x);

                match (vapor, liquid) {

                    (Ok(vapor), Ok(liquid)) => {
                        
                        if vapor.gibbs() < liquid.gibbs() {
                            Ok(vapor)
                        } else {
                            Ok(liquid)
                        }
                    }

                    (Ok(state), Err(_)) | (Err(_), Ok(state)) => Ok(state),
                    (Err(e1), Err(_e2)) => Err(e1),
                    
                }
            }
            
        }
        
    }

    pub fn new_tp(eos:Arc<E<R>>, t:f64, p:f64, phase:Option<DensityInitialization>)->StateResult<R>{

        State::new_tpx(eos, t, p, array![1.0], phase)
    }

    pub fn new_tr(eos:Arc<E<R>>, t:f64, d:f64)->StateResult<R>{

        Ok(State::new_trx(eos, t, d, array![1.0]))
    }

}

impl<R:Residual> State<R> {

    pub fn pressure(&self)->f64{
        self.p
    }

    pub fn composition(&self)->&[f64]{
        self.x.as_slice().unwrap()
    }

    pub fn temperature(&self)->f64{
        self.t
    }

    pub fn density(&self)->f64{
        self.d
    }

    pub fn lnphi(&self)->Array1<f64>{
        
        self.eos.lnphi(self.t, self.d, &self.x)
    }

    pub fn helmholtz(&self)->f64{
        self.eos.helmholtz_isov(self.t, self.d, &self.x)
    }

    pub fn entropy(&self)->f64{
        self.eos.entropy(self.t, self.d, &self.x)
    }

    pub fn chem_pot(&self)->Array1<f64>{
        self.eos.chem_pot(self.t, self.d, &self.x)
    }

    pub fn gibbs(&self) ->f64{
        self.eos.gibbs(self.t, self.d, &self.x)
    }

    pub fn max_density(&self)->f64{

        self.eos.max_density(&self.x)
    }


    pub fn compressibility(&self)->f64{

        self.eos.compressibility(self.t, self.d, &self.x)

    }

    pub fn molar_weight(&self) -> &Array1<f64> {

        self.eos.molar_weight()

    }
    pub fn get_properties(&self)-> &Properties{

        self.eos.get_properties()
        
    }
    pub fn mass_density(&self)->f64{

        let vmw = self.eos.molar_weight();
        let x = &self.x;

        let mw = vmw.dot(x);

        self.d * mw

    }
    
}


#[cfg(test)]
mod tests {

    use approx::assert_relative_eq;
    use ndarray::array;

    use crate::{models::{cpa::{CPA, parameters::readyto::*, rdf::Kontogeorgis}, cubic::models::SRK}, parameters::Parameters};
    use super::*;
    
    #[test]
    fn liq_stable() {

        let water = water4c();

        let p = CPAParameters::new(vec![water], vec![], SRK.into());
        let r = CPA::<Kontogeorgis>::from_parameters(p);

        let eos = EquationOfState::from_residual(r).into();

        let t = 298.15;
        let p = 1e5;
        let x = array![1.0];
        
        let stable = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), None).unwrap();
        let vapor = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), Some(DensityInitialization::Vapor)).unwrap();
        let liquid = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), Some(DensityInitialization::Liquid)).unwrap();
        

        assert_relative_eq!(stable.d, liquid.d, epsilon = 1e-10);
        assert!( vapor.gibbs() > liquid.gibbs() );


    }

    #[test]
    fn vap_stable() {

        let water = water4c();

        let p = CPAParameters::new(vec![water], vec![], SRK.into());
        let r = CPA::<Kontogeorgis>::from_parameters(p);

        let eos = EquationOfState::from_residual(r).into();

        let t = 298.15 + 100.;
        let p = 1e5;
        let x = array![1.0];
        
        let stable = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), None).unwrap();
        let vapor = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), Some(DensityInitialization::Vapor)).unwrap();
        let liquid = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), Some(DensityInitialization::Liquid)).unwrap();


        assert_relative_eq!(stable.d, vapor.d, epsilon = 1e-10);
        assert!( vapor.gibbs() < liquid.gibbs() );

    }

    #[test]
    fn invalid_pressure() {
            let water = water4c();

        let p = CPAParameters::new(vec![water], vec![], SRK.into());
        let r = CPA::<Kontogeorgis>::from_parameters(p);

        let eos = EquationOfState::from_residual(r).into();

        let t = 298.15;
        let p = 0.0;
        let x = array![1.0];
        
        let s = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), Some(DensityInitialization::Vapor));

        assert!(s.is_err());
        // println!("Density = {}", s.d);
        
    }
    #[test]
    fn invalid_frac() {
            let water = water4c();

        let p = CPAParameters::new(vec![water], vec![], SRK.into());
        let r = CPA::<Kontogeorgis>::from_parameters(p);

        let eos = EquationOfState::from_residual(r).into();

        let t = 298.15;
        let p = 1e5;
        let x = array![-12.0];
        
        let s = State::new_tpx(Arc::clone(&eos), t, p, x.clone(), Some(DensityInitialization::Vapor));

        assert!(s.is_err());
        // println!("Density = {}", s.d);
        
    }

    
}