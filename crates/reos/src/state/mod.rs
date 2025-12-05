// mod's
pub mod eos;
pub mod density_solver;




use ndarray::{Array1, Array2};
use std::fmt::{Display, Write};
use std::sync::Arc;

use crate::models::cpa::CPA;
use crate::models::IDEAL_GAS_CONST;
use crate::models::cpa::rdf::ElliotRDF;
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
pub fn fmt_array2<W: Write>(
    f: &mut W,
    name: &str,
    mat: &Array2<f64>,
    decimals: usize,
) -> std::fmt::Result {
    writeln!(f, "{} = [", name)?;
    for row in mat.rows() {
        write!(f, "    [")?;
        for (j, val) in row.iter().enumerate() {
            if j > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:.*}", decimals, val)?;
        }
        writeln!(f, "],")?;
    }
    writeln!(f, "]")?;
    Ok(())
}
// 
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
    pub fn composition(&self)->Array1<f64>{
        self.x.clone()
    }
    pub fn temperature(&self)->f64{
        self.t
    }

}
