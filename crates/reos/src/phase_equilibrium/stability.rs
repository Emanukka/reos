use ndarray::{ Array1};

use crate::{phase_equilibrium::PhaseEquilibrium, residual::Residual, state::{density_solver::DensityInitialization, eos::{EosError, EosResult}, State}, tools::newton};


impl <R:Residual> State<R> {

    pub fn tpd(
        &self,
        x:Array1<f64>,
        xphase:DensityInitialization)
        ->EosResult<f64>{
        
        let hy = self.lnphi()?+self.x.ln();
        let lnx=x.ln();
        let daughter_phase = State::new_tpx(&self.eos, self.t, self.p, x, xphase )?;
        let hx=daughter_phase.lnphi()?+&lnx;
        Ok(
        (lnx.exp()*(hx-hy)).sum()
        )
        

    }
    pub fn min_tpd(
        &self,
        xphase:DensityInitialization,
        xguess:Array1<f64>,
        tol:Option<f64>,
        it_max:Option<i32>)
        ->EosResult<(f64,State<R>)>{
        let tol =tol.unwrap_or(1e-8);
        let it_max=it_max.unwrap_or(200);
        let mut it=0;
        let mut res=10.0;

        // Cte
        let p=self.p;
        let z=&self.x;
        let t=self.t;
        let hy = self.lnphi()?+z.ln();

        //Var
        let mut daughter_phase=State::new_tpx(&self.eos, t, p, xguess.clone(), xphase)?;

        let mut old_dens=daughter_phase.rho;
        let mut vw_old:Array1<f64>;
        let mut vw =xguess.clone();
        let mut x=xguess;

        // xPhase guess
        let mut lnphix:Array1<f64>=Array1::zeros(x.len());

        while (res>tol) & (it<it_max) {

            vw_old = vw;
            // incipient state
            daughter_phase = State::new_tpx(&self.eos, t, p, x, DensityInitialization::Guess(old_dens) )?;
            old_dens=daughter_phase.rho;
            lnphix=daughter_phase.lnphi()?;

            vw = (&hy-&lnphix).exp();
            
            x = &vw/vw.sum();

            res = (&vw - vw_old).pow2().sum().sqrt();
            it+=1
        }
        let hx = x.ln()+ lnphix;
        let dg = (vw*(hx-hy)).sum();

        Ok((dg,daughter_phase))

    }
    

}

pub enum TemperatureOrPressure {
    Temperature(f64),
    Pressure(f64)
}
impl<R:Residual> PhaseEquilibrium<R>{

    // pub fn root_min_tpd(
    //     &self,
    //     temperature_or_pressure:TemperatureOrPressure,
    //     x0:f64,
    //     zphase:DensityInitialization,
    //     z:Array1<f64>,
    //     xguess:Array1<f64>,
    //     )->[State<R>;2]{

    //     let mut p:f64;
    //     let mut t: f64;
    //     let xphase= zphase.inverse();
    //     match temperature_or_pressure {
            
    //         TemperatureOrPressure::Temperature(t)=>{

    //             let delta_g=|p|->f64{
                
    //             let state=State::new_tpx(&self.eos, t, p, z.clone(), zphase);

    //             state.unwrap().min_tpd(xphase, xguess, None, None).unwrap().0
            
    //             };
    //             match newton(delta_g, x0, Option(1e-6), Option(100)){
            
    //             Ok(x)=>{
    //                 let result=x.x;
    //                 let state=State::new_tpx(eos, t, result, z, zphase).unwrap();
    //                 state.min_tpd(xphase, xguess, tol, it_max).unwrap()
    //             }
    //             Err(

    //             )
    //     }    
    //         }
    //         TemperatureOrPressure::Pressure(p)=>{
    //             p=p;
    //             t=x0
    //         }
    //     }

        




    // }

    pub fn tpd(
        &self,
        t:f64,
        p:f64,
        z:Array1<f64>,
        xphase:DensityInitialization,
        xguess:Array1<f64>,
        tol:Option<f64>,
        it_max:Option<i32>
        )->Result<(f64,Array1<f64>),EosError>{
        let tol =tol.unwrap_or(1e-8);
        let it_max=it_max.unwrap_or(200);
        let mut it=0;
        let mut res=10.0;

        // Cte
        let mother_phase=State::new_tpx(&self.eos, t, p, z.clone(), xphase.inverse())?;
        let mut daughter_phase=State::new_tpx(&self.eos, t, p, xguess.clone(), xphase)?;

        let mut old_dens=daughter_phase.rho;
        let hy = mother_phase.lnphi()?+z.ln();
        let mut vw_old:Array1<f64>;
        let mut vw =xguess.clone();
        let mut x=xguess;

        // xPhase guess
        let mut lnphix:Array1<f64>=Array1::zeros(x.len());

        while (res>tol) & (it<it_max) {

            vw_old = vw;
            // incipient state
            daughter_phase = State::new_tpx(&self.eos, t, p, x, DensityInitialization::Guess(old_dens) )?;
            old_dens=daughter_phase.rho;
            lnphix=daughter_phase.lnphi()?;

            vw = (&hy-&lnphix).exp();
            
            x = &vw/vw.sum();

            res = (&vw - vw_old).pow2().sum().sqrt();
            it+=1
        }
        let hx = x.ln()+ lnphix;
        let dg = (vw*(hx-hy)).sum();

        //Retorna TPDresult (DG,StateMae,StateFilha)
        
        Ok((dg,x))

    }

}
#[cfg(test)]
pub mod tests{
    use std::sync::Arc;

    use ndarray::{array, Array1};
 
    use crate::{models::cpa::parameters::water_acetic_acid, phase_equilibrium::PhaseEquilibrium, state::{density_solver::DensityInitialization::{Vapor,Liquid}, State}};

    #[test]
    fn verify_tpd_close_to_bbpoint(){

        let eos = water_acetic_acid();
        let eos=Arc::new(eos);

        let peq=PhaseEquilibrium::new(
            &eos,
            None);
            // Some(antoine_water_acetic_acid()));
        //State Variables
        let t=313.15;
        let x=Array1::from_vec(vec![0.5,0.5]);
        
        let (pb,y)=peq.bbpy(t, x.clone(),Some(1e-10),Some(1e-10)).unwrap();
        
        println!("{}",&y);
        let p_mais =pb*1.2;
        let p_menos =pb*0.2;

        println!("Pbolha= {} bar",pb/1e5);
        println!("Pmais= {} bar",p_mais/1e5);
        println!("Pmenos= {} bar",p_menos/1e5);

        let zphase=State::new_tpx(&eos, t, pb, x.clone(), Liquid).unwrap();

        // let dg_zero= peq
        // .tpd(
        //     t, 
        //     pb,
        //     x.clone(),
        //     Vapor, 
        //     &y*1., 
        //     Some(1e-10),
        //     None).unwrap();
        let dg_zero=zphase.min_tpd(Vapor, array![0.5,0.5], None, None).unwrap();
        println!("ΔG em pBolha={}",dg_zero.0);
        println!("y TPD={}",dg_zero.1.x);

        let dg_mais= peq
        .tpd(
            t, 
            p_mais,
            x.clone(),
            Vapor , 
            &y*1., 
            Some(1e-10),
            None).unwrap();
        
        println!("ΔG em P>Pbolha ={}",dg_mais.0);

        let dg_menos= peq
        .tpd(
            t, 
            p_menos,
            x.clone(),
            Vapor , 
            &y*1., 
            Some(1e-10),
            None).unwrap();
        
        println!("ΔG em P<Pbolha ={}",dg_menos.0);

    }

}

