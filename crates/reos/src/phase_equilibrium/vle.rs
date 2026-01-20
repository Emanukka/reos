use ndarray::{ Array1};

use crate::{phase_equilibrium::{Antoine, AntoineRecord, PhaseEquilibrium}, residual::Residual, state::{density_solver::DensityInitialization, eos::EosError, StateResult, S}};

use crate::phase_equilibrium::LogBase;
impl<R:Residual> PhaseEquilibrium<R>{
    
    pub fn vapor(&self,t:f64,p:f64,y:Array1<f64>)
    ->StateResult<R>{
        // S::new_tpx(self.0.clone(), t, p,y, DensityInitialization::Vapor)
        S::new_tpx(self.eos.clone(), t, p, y, Some(DensityInitialization::Vapor))
    }
    pub fn liquid(&self,t:f64,p:f64,x:Array1<f64>)
    ->StateResult<R>{
        S::new_tpx(self.eos.clone(), t, p, x, Some(DensityInitialization::Liquid))
    }

    pub fn set_antoine(&mut self,antoine:Antoine){
        self.correl=Some(antoine)
    }

    pub fn pyraoult(&self,t:f64,x:&Array1<f64>)->Result<(f64,Array1<f64>),EosError>{
        
        // let psat=Array1::from_elem(x.len(), 1e5);[]
        let mut psat=Array1::from_elem(x.len(), 1e5);

        match &self.correl {
            Some(a)=>{psat=a.psat(t)}
            None=>{}
        }
        let p=psat.dot(x);
        let y=(psat*x)/p;
        Ok((p,y))
    }

    //bbpy:
    //bbty: 
    //orvpx
    pub fn bbpy(&self,
        t:f64,
        x:Array1<f64>,
        tol_p:Option<f64>,
        top_y:Option<f64>,
    )->Result<(f64,Array1<f64>),EosError>{


        let (mut p0,mut y0)=self.pyraoult(t, &x)?;
        // let mut p0=1e5;
        // let mut y0=Array1::from_vec(vec![0.5,0.5]);
        // dbg!(p0,&y0);
        let mut p1=p0*1.0;
        let mut y1=y0*1.0;
        
        // let liquid=self.liquid(t, p1, x.clone())?;
        // let vapor =self.vapor(t, p1, y1.clone())?;


        let mut res_p: f64=10.0; // Pressure loop (j)

        let tol_p=tol_p.unwrap_or(1e-7);
        let tol_y=top_y.unwrap_or(1e-6);

        let max_j=100;
        let max_i=100;

        let mut j=0;

        while (res_p.abs()>tol_p) & (j<max_j){

            
            let mut i =0;
            let mut res_y: f64=10.0; // y loop (i)

            p0=p1;
            // println!("LIQUID");

            let liquid = self.liquid(t, p0, x.clone())?;

            let phi_l = liquid.lnphi().exp();
            // dbg!(j);
            // dbg!(p0);
            // Reduzir n de iterações ou toly
            // dbg!(&y1,&phi_l,&x);
            while (res_y.abs()>tol_y) &(i<max_i) {

                
                y0=y1.clone();
                let ynorm = &y1/y1.sum();

                // println!("VAPOR");
                let vapor =self.vapor(t, p0, ynorm)?;

                let phi_v=vapor.lnphi().exp();

                y1=(&phi_l/phi_v)*&x;

                res_y= ((&y1-&y0)).mapv(|u|u*u).sum().sqrt();
                i+=1;
                
            }

            // dbg!(i);
            // dbg!(res_y);
            let sum_y=y1.sum();
            res_p=sum_y-1.0;
            // println!("{}",p0);
            p1=p0*sum_y;
            j+=1;

        }

        // println!("{}",j);
        if j<max_j{

            Ok(
                (p1,y1)
            )
        }else {
            Err(EosError::NotConverged("Bubble Point, max it".to_string()))
        }
    }

}

#[allow(unused)]

fn antoine_water_acetic_acid()->Antoine{

    let ant1 = AntoineRecord::record(6.20963,2354.731,7.559,LogBase::Log10 );
    let ant2 = AntoineRecord::record(4.68206,1642.54,-39.764,LogBase::Log10 );
    Antoine::from_records(vec![ant1,ant2])
}

#[cfg(test)]
pub mod tests{
    use std::sync::Arc;

    use approx::assert_relative_eq;
    use ndarray::Array1;

    use crate::{ models::{cpa::{CPA, SCPA, parameters::readyto::{acetic1a, water4c, water4c_acetic1a}, rdf::Kontogeorgis}, cubic::models::SRK}, phase_equilibrium::PhaseEquilibrium, state::E};
    #[allow(unused_imports)]
    // use crate::{models::{ cpa::{CPA,parameters::{octane_acoh,acoh_octane, methanol_2b, methanol_3b, water_acetic_acid}}, cubic::Cubic}, phase_equilibrium::{vle::antoine_water_acetic_acid, Antoine, AntoineRecord, LogBase, PhaseEquilibrium}, state::E};

    
    fn water_acetic_acid()-> CPA<Kontogeorgis>{
        let water = water4c();
        let acetic = acetic1a();
        let b = water4c_acetic1a();
        let cpa = SCPA::from_records(vec![water,acetic],vec![b], SRK.into());
        cpa
    }

    #[test]
    pub fn cmp_bbpy_water_acoh(){

        println!("bbpy-water-acoh");

        let eos = E::from_residual(water_acetic_acid());
        let peq=PhaseEquilibrium::new(
            &Arc::new(eos),
            None);
            // Some(antoine_water_acetic_acid()));
        //State Variables
        let t = 300.0;
        let x=Array1::from_vec(vec![0.5,0.5]);
        
        let (pb,y)=peq.bbpy(t, x,Some(1e-7),Some(1e-7)).unwrap();
        
        let cmp=[0.59739271, 0.40260708];

        assert_relative_eq!(pb,3127.2493944115,epsilon=1e-4);
        
        assert_relative_eq!(y[0],cmp[0],epsilon = 1e-6);
        assert_relative_eq!(y[1],cmp[1],epsilon = 1e-6);


    }
 
    // #[test]
    // pub fn bench(){

    //     // println!("bbpy-water-acoh");
    //     let n = 100;
    //     let eos = E::from_residual(water_acetic_acid());
    //     let peq=PhaseEquilibrium::new(
    //         &Arc::new(eos),
    //         None);
    //         // Some(antoine_water_acetic_acid()));
    //     //State Variables
    //     let t = 300.0;
    //     let x=Array1::from_vec(vec![0.5,0.5]);
        
    //     for _ in 0..n{
    //         let (_,_)=peq.bbpy(t, x.clone(),Some(1e-7),Some(1e-7)).unwrap();
            
    //     }

    // }
 
    // #[test]
    // pub fn cmp_bbpy_acoh_octane(){

    //     println!("bbpy-acoh-octane");
    //     let eos = acoh_octane();
    //     let peq=PhaseEquilibrium::new(
    //         &Arc::new(eos),
    //         None);
    //     //State Variables
    //     let t=343.15;
    //     let x1=0.6;
    //     let x=Array1::from_vec(vec![x1,1.0-x1]);
        
    //     let (pb,y)=peq.bbpy(t, x,Some(1e-6),Some(1e-6)).unwrap();
        
    //     let cmp=[0.6542769760724502, 0.34572302386625076];

    //     assert_relative_eq!(pb,28987.67376094271,epsilon=1e-2);
        
    //     assert_relative_eq!(y[0],cmp[0],epsilon=1e-6);
    //     assert_relative_eq!(y[1],cmp[1],epsilon=1e-6);


    // }
    // #[test]
    // pub fn cmp_bbpy_octane_acoh(){


    //     println!("bbpy-octane-acoh");
    //     let eos = octane_acoh();
    //     let peq=PhaseEquilibrium::new(
    //         &Arc::new(eos),
    //         None);
    //     //State Variables
    //     let t=343.15;

    //     let x1=0.4;
    //     let x=Array1::from_vec(vec![x1,1.0-x1]);
        
    //     let (pb,y)=peq.bbpy(t, x,Some(1e-6),Some(1e-6)).unwrap();
        
    //     let cmp=[ 0.34572302386625076,0.6542769760724502];
        
    //     // println!("P bol ={}",pb);
    //     // println!("y={}",y);

    //     assert_relative_eq!(pb,28987.67376094271,epsilon=1e-2);
        
    //     assert_relative_eq!(y[0],cmp[0],epsilon=1e-6);
    //     assert_relative_eq!(y[1],cmp[1],epsilon=1e-6);


    // }

    // #[test]
    // pub fn vle(){

        
    //     cmp_bbpy_acoh_octane();

    //     cmp_bbpy_octane_acoh();
        
    //     cmp_bbpy_water_acoh();

    // }
}
