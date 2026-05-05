use std::sync::Arc;

use ndarray::Array1;

use crate::{models::R_GAS, state::{E, Residual, State, StateResult, eos::EosError}};


#[derive(Clone, Copy)]

pub enum Phase{
    Vapor,
    Liquid,
    Unkown
}

impl std::fmt::Display for Phase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Phase::Vapor=> write!(f,"vapor"),
            Phase::Liquid=> write!(f,"liquid"),
            Phase::Unkown=> write!(f,"unkown"),
        }
    }
}
#[derive(Clone, Copy)]
pub enum DensityInitialization {
    Vapor,
    Liquid,
    Stable,
    Guess(f64)
}

impl Into<Phase> for DensityInitialization {
    fn into(self) -> Phase {
        match self {
            DensityInitialization::Vapor=>{
                Phase::Vapor
            }
            DensityInitialization::Liquid=>{
                Phase::Liquid
            }
            DensityInitialization::Guess(_)=>{
                Phase::Vapor
            }

            _ => {
                Phase::Unkown
            }
        }
    }
}
impl DensityInitialization {
    pub fn from_str(s:&str)->Result<Self,EosError>{

        if s.to_lowercase() == "vapor"{
            Ok(DensityInitialization::Vapor)
        }else if s.to_lowercase()=="liquid" {
            Ok(DensityInitialization::Liquid)
        }else if s.to_lowercase()=="stable" {
            Ok(DensityInitialization::Stable)
        }else {
            Err(EosError::PhaseError)
        }
    }

    // pub fn inverse(&self)->DensityInitialization{

    //     match *self {
    //         DensityInitialization::Vapor=>{
    //             Self::Liquid
    //         }
    //         DensityInitialization::Liquid=>{
    //             Self::Vapor
    //         }
    //         Self::Guess(_)=>{
    //             Self::Vapor
    //         }

    //     }
    // }
}


pub fn density<R:Residual>(eos:Arc<E<R>>,t:f64,p:f64,x:Array1<f64>,guess:f64)-> StateResult<R>{           

    const EPS:f64 = 1e-5;
    const DEN:f64 = 2.0 * EPS;

    let rhomax = eos.max_density(&x);
    
    let f = |s:f64| { 

        let rho = s * rhomax;
        let p_iter = eos.pressure(t,rho,&x);
        
        // eprintln!("p_iter:{}", p_iter);
        (1.0 - s) * (p_iter - p) 

    };

    let dfds = |s:f64|{

        let s_for = s + EPS;
        let s_back = s - EPS;
        
        let rho_for = s_for * rhomax;
        let rho_back = s_back * rhomax;
        
        let p_for = eos.pressure(t, rho_for, &x);
        let p_back = eos.pressure(t, rho_back, &x);
        
        let forward = (1.0 - s_for ) * (p_for - p) ;
        let backward = (1.0 - s_back ) * (p_back - p) ;
        
        (forward - backward) / DEN

    };
    
    let tol = 1e-8;
    let it_max = 100;
    
    let mut f0 = 1.;
    let mut df_at_s0 = 1.;
    let mut s_min = 0.0;
    let mut s_max = 1.0;
    
    let mut it = 0;
    let mut res: f64 = 1.0;
    let mut s1 = guess;

    while (res.abs() > tol) & (it < it_max) {
        
        // debug
        let s0 = s1;
        f0 = f(s0);
        df_at_s0 = dfds(s0);
        s1 = s0 - f0 / df_at_s0;
        
        it += 1;

        let bis = bissec(s0, s1, f0, &mut s_max, &mut s_min, &mut res);

        // eprintln!("Iteration: {}, s1: {}, s_min: {}, s_max: {}, s0: {}, f0: {}, df_at_s0: {}, res: {}", it, s1, s_min, s_max, s0, f0, df_at_s0, res);

        match bis{
            Bissection::In=>{continue;}
            Bissection::Out(x)=>{
                s1 = x
            }
        }

    }
    

    if it == it_max{

        return Err(EosError::NotConverged("density solver: max iterations".to_string()))

    }
    
    else if res.is_nan(){

        return Err(EosError::NotConverged(("density solver: nan at it=".to_string()+&it.to_string()).to_string()));

    }
    else {

        let density = rhomax * s1;
        Ok(State::new_trx(eos, t, density, x))
    }
}



pub fn bissec(
    s0:f64,
    s1:f64,
    f0:f64,
    s_max:&mut f64,
    s_min:&mut f64,
    res:&mut f64,

)->Bissection{

    if f0 > 0.0 {
        
        *s_max = s0;
    
    } else if f0 < 0.0 {
        
        *s_min = s0;
    
    }
    
    *res = ((s1 - s0) / s0).abs();
    
    if (s1 >= *s_min) & (s1 <= *s_max) {
    
        Bissection::In
    
    } else {
        
        let s1 = 0.5 * (*s_max + *s_min);

        Bissection::Out(s1)

    }
}

pub enum Bissection{
    Out(f64),
    In
}

pub fn vapor_solver<R:Residual>(eos:Arc<E<R>>, t:f64, p:f64, x:Array1<f64>) -> StateResult<R> {

    let bm = 1. / eos.max_density(&x);
    let guess = bm / (bm + (R_GAS * t) / p);
    density(eos, t, p, x, guess)

}

pub fn liquid_solver<R:Residual>(eos:Arc<E<R>>, t:f64, p:f64, x:Array1<f64>) -> StateResult<R> {

    let guess = 0.99;
    density(eos, t, p, x, guess)

}