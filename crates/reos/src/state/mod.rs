// mod's
pub mod eos;
pub mod density_solver;

use ndarray::{Array1};
use std::sync::Arc;

use crate::models::IDEAL_GAS_CONST;
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult, EquationOfState};
use crate::state::density_solver::{DensityInitialization,density};
pub type E<R> = EquationOfState<R>;
pub type S<R> = State<R>;

pub type StateResult<R> = Result<S<R>,EosError>;


#[derive(Clone)]
pub struct State<R:Residual>{

    pub eos:Arc<EquationOfState<R>>,
    // Kelvin
    pub t:f64,
    pub p:f64,
    pub d:f64,
    pub x:Array1<f64>,
}



impl<R:Residual> State<R> {

    pub fn new_trx(eos:Arc<E<R>>,t:f64,d:f64,x:Array1<f64>)-> Self{

        let p = eos.pressure(t, d, &x);

        Self{
            eos:Arc::clone(&eos),
            t,
            p,
            d,
            x,
        }


    }
    pub fn new_tpx(eos:Arc<E<R>>,t:f64,p:f64,x:Array1<f64>,phase:DensityInitialization)->StateResult<R>{

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

                let bm = 1./eos.residual.max_density(&x);
                let guess = bm/(bm + (IDEAL_GAS_CONST*t)/p);

                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }
            }

            DensityInitialization::Guess(old_density)=>{

                let bm = 1./eos.residual.max_density(&x);
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


//properties methods
impl<R:Residual> State<R> {
    
    pub fn lnphi(&self)->Array1<f64>{
        
        self.eos.lnphi(self.t, self.d, &self.x)
        // self.eos.lnphi(self.t, self.rho, &self.x)
    }
    pub fn pressure(&self)->f64{
        // self.eos.pressure(self.t, self.rho, &self.x)
        // self.eos.pressure(self.t, self.rho, &self.x)
        self.p
    }

    pub fn max_density(&self)->f64{

        self.eos.residual.max_density(&self.x)
    }
    pub fn composition(&self)->&[f64]{
        self.x.as_slice().unwrap()
    }
    pub fn temperature(&self)->f64{
        self.t
    }
    
}
