use std::{f64::NAN, sync::Arc};

use ndarray::Array1;

use crate::state::{eos::EosError,Residual, State, StateResult, E};


#[derive(Clone, Copy)]

pub enum Phase{
    Vapor,
    Liquid
}
// guess vai ser um tipo especifico (vapor ou liquido)

#[derive(Clone, Copy)]
pub enum DensityInitialization {
    Vapor,
    Liquid,
    Guess(f64)
}

impl DensityInitialization {
    pub fn from_str(s:&str)->Result<Self,EosError>{

        if s.to_lowercase()=="vapor"{
            Ok(DensityInitialization::Vapor)
        }else if s.to_lowercase()=="liquid" {
            Ok(DensityInitialization::Liquid)
        }else {
            Err(EosError::PhaseError)
        }

    }

    pub fn inverse(&self)->DensityInitialization{

        match *self {
            DensityInitialization::Vapor=>{
                Self::Liquid
            }
            DensityInitialization::Liquid=>{
                Self::Vapor
            }
            Self::Guess(_)=>{
                Self::Vapor
            }

        }
    }
}


pub fn density<R:Residual>
    (
    eos:&Arc<E<R>>,
    t:f64,
    p:f64,
    x:Array1<f64>,
    guess:f64
    )-> StateResult<R>
{           


        let rhomax = 1./eos.residual.bmix(&x);
        let eps = 1e-5;

        let f = |s:f64| { 
        // let f = |s:f64| { 
            let rho = s*rhomax;
            let p_iter = eos.pressure(t,rho,&x).unwrap_or(NAN);

            (1.0-s)*(p_iter - p) 

        };

        let dfds = |s:f64|{

            let s_mais = s+eps;
            let s_menos = s-eps;

            let rho_mais = s_mais*rhomax;
            let rho_menos = s_menos*rhomax;

            let pfwd=eos.pressure(t, rho_mais, &x).unwrap_or(NAN);
            let pbwd=eos.pressure(t, rho_menos, &x).unwrap_or(NAN);

            let fwd = (1.0 - (s+eps) )*(pfwd - p) ;
            let bwd = (1.0 - (s-eps) )*(pbwd - p) ;

            (fwd - bwd)/(2.*eps)

        };



        let mut f0: f64 = 1.;

        let mut s_min =0.;
        let mut s_max =1.;

        let tol = 1e-7;

        let it_max = 100;

        let mut it =0;

        let mut df_at_s0:f64;
        let mut res: f64 = 1.;

        let mut s1= guess;

        while (res.abs()>tol) & (it<it_max) {

            let s0: f64 = s1;
            // println!("s={}",s0);

            // (f0,df_at_s0) = f_and_dfds(s0);

            f0 = f(s0);
        // dbg!("aqui");

            df_at_s0=dfds(s0);
            s1 = s0 - f0/df_at_s0;
            
            it+=1;

            let bis=bissec(s0, s1, f0,&mut s_max,&mut s_min, &mut res);
            match bis{
                Bissection::In=>{continue;}

                Bissection::Out(x)=>{
                    s1=x
                }
            }

        }

        // println!("it d={it}");
        if it==it_max{
            return Err(EosError::NotConverged("density".to_string()))
        }
        
        else if res.is_nan(){
            return Err(EosError::NotConverged(("density (NaN value) at it=".to_string()+&it.to_string()).to_string()));

        }

        else {
            let density= rhomax*s1;
            // println!("iterações solver={}",it);

            Ok(State::new_trx(eos, t, density, x).unwrap())
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

    if f0>0.{
        *s_max = s0;
    }
    else if f0<0.{
        *s_min = s0;
    }
    
    *res = ((s1-s0)/s0).abs();
    
    if (s1>=*s_min) & (s1<=*s_max){
        Bissection::In
    }
    else {
        let s1 = (*s_max+*s_min)/2.;
        Bissection::Out(s1)

    }
}

pub enum Bissection{
    Out(f64),
    In
}