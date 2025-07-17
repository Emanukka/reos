// mod's
pub mod eos;
pub mod density_solver;




use ndarray::{Array1, Array2};
use std::sync::Arc;

use crate::models::IDEAL_GAS_CONST;
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult, EquationOfState};
use crate::state::density_solver::{DensityInitialization,density};
pub type E<R>=EquationOfState<R>;
pub type S<R>=State<R>;

pub type StateResult<R>=Result<S<R>,EosError>;


#[derive(Clone)]
pub struct State<R:Residual>{

    pub eos:Arc<EquationOfState<R>>,
    // Kelvin
    pub t:f64,
    pub p:f64,
    pub rho:f64,
    pub x:Array1<f64>,
    // pub cache: Option<Cache>
}


impl<R:Residual> State<R> {

    pub fn new_trx(eos:&Arc<E<R>>,t:f64,rho:f64,x:Array1<f64>)->StateResult<R>{

        let p =eos.pressure(t, rho, &x);

        match p {
            Ok(pres)=>Ok(
                Self{
                    eos:Arc::clone(&eos),
                    t,
                    p:pres,
                    rho,
                    x,
                }
            ),
            Err(e)=>Err(e)
        }

    }
    pub fn new_tpx(eos:&Arc<E<R>>,t:f64,p:f64,x:Array1<f64>,phase:DensityInitialization)->StateResult<R>{

        match phase {
            DensityInitialization::Liquid=>{

                let guess=0.99;
                // let guess=0.9;
                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }

            }

            DensityInitialization::Vapor=>{

                let bm=eos.residual.bmix(&x);
                let guess = bm/(bm + (IDEAL_GAS_CONST*t)/p);

                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }
            }

            DensityInitialization::Guess(old_density)=>{
                let bm=eos.residual.bmix(&x);
                let guess=old_density*bm;

                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }
            }
            
        }
        
    }




}


impl<R:Residual> State<R> {
    
    pub fn lnphi(&self)->EosResult<Array1<f64>>{
        self.eos.lnphi(self.t, self.rho, &self.x)
    }
    pub fn pressure(&self)->EosResult<f64>{
        self.eos.pressure(self.t, self.rho, &self.x)
    }
    pub fn bmix(&self)->f64{
        self.eos.residual.bmix(&self.x)
    }
}



#[cfg(test)]
mod tests{
    use std::sync::Arc;

    use approx::{assert_relative_eq, assert_relative_ne, relative_eq};
    use ndarray::{array, Array1};

    use crate::{models::{associative::{}, cpa::{self, CPA}, cubic::{Cubic}}, parameters::Parameters, state::{eos::EquationOfState, State,Residual}};


    // comparar valores com PYTHON

    // #[test]

    fn water_sat_100(){

        // // let t=300.0;
        // let vp = Array1::linspace(1.0, 600., 100);
        // // let mut yw = Array1::zeros(100);
        // // let eos = co2_and_water(*p, t, 0.5, 0.5);
        // // let eos = co2_and_water();
        // // let comps = vec!["co2","water"];
        // let comps = vec!["co2","ch4","water"];
        // let set_1= "konteo_1_mod.json";

        // let cubic = CubicParameters::from_file(set_1, comps.clone());
        // let assoc= ASCParameters::from_file(set_1, comps.clone());

        // let eos = CPAeos::<SRKeos>::new(cubic, assoc);
        // // println!("{}",&cubic);

        // let t = 273.15 +3.0;
        // let p =325e5;
        // let x=vec![0.03,0.97,1e-6];
        // let x=Array1::from_vec(x);
        // let state=State::new_tpx(Arc::new(eos),  t,p, x, Vapor).unwrap();

        // println!("{}",state.water_saturation().unwrap())

    }
}