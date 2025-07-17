use std::sync::Arc;

use ndarray::{Array, Array1};

use crate::{phase_equilibrium::{PhaseEquilibrium}, residual::Residual, state::{density_solver::DensityInitialization, eos::EosError, State, StateResult, E, S}, tools::newton};

impl<R:Residual> PhaseEquilibrium<R>{

    pub fn tpd(
        &self,
        t:f64,
        p:f64,
        z:Array1<f64>,
        xphase:DensityInitialization,
        xguess:Array1<f64>,
        tol:Option<f64>,
        it_max:Option<i32>
        )->Result<f64,EosError>{
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
        Ok(dg)

    }

}
pub mod tests{
    use std::sync::Arc;

    use ndarray::Array1;

    use crate::{models::cpa::water_acetic_acid, phase_equilibrium::PhaseEquilibrium};
    use crate::state::density_solver::DensityInitialization::{Liquid,Vapor};

    #[test]

    fn tpd_1(){

        let eos = water_acetic_acid();
        let peq=PhaseEquilibrium::new(
            Arc::new(eos),
            None);
            // Some(antoine_water_acetic_acid()));
        //State Variables
        let t=313.15;
        let x=Array1::from_vec(vec![0.5,0.5]);
        
        let (pb,y)=peq.bbpy(t, x.clone(),Some(1e-10),Some(1e-10)).unwrap();
        
        // println!("{}",&y);
        let p_mais =pb*1.2;
        let p_menos =pb*0.2;

        println!("Pbolha= {} bar",pb/1e5);
        println!("Pmais= {} bar",p_mais/1e5);
        println!("Pmenos= {} bar",p_menos/1e5);
        let dg_zero= peq
        .tpd(
            t, 
            pb,
            x.clone(),
            Vapor, 
            &y*1., 
            Some(1e-10),
            None).unwrap();
        
        println!("ΔG em pBolha={}",dg_zero);

        let dg_mais= peq
        .tpd(
            t, 
            p_mais,
            x.clone(),
            Vapor , 
            &y*1., 
            Some(1e-10),
            None).unwrap();
        
        println!("ΔG em P>Pbolha ={}",dg_mais);

        let dg_menos= peq
        .tpd(
            t, 
            p_menos,
            x.clone(),
            Vapor , 
            &y*1., 
            Some(1e-10),
            None).unwrap();
        
        println!("ΔG em P<Pbolha ={}",dg_menos);

    }

}
// impl <E:EquationOfState> State<E> {

//     pub fn tpd(&self,xphase:DensityInitialization,xguess:Array1<f64>)->Result<f64,EosError>{

//         let tol =1e-8;
//         let it_max=200;
//         let mut it=0;
//         let mut res=10.0;

//         // Cte
//         let hy = self.lnphi()?+self.x.ln();

//         let mut vw_old:Array1<f64>;
//         let mut vw =xguess.clone();
//         let mut x=xguess;

//         // xPhase guess
//         let mut lnphix:Array1<f64>=Array1::zeros(x.len());
//         let mut old_density = State::
//         new_tpx(self.eos.clone(), self.t, self.p, x.clone(), xphase)?.rho;

//         while (res>tol) & (it<it_max) {

//             vw_old = vw;

//             // incipient state 
//             let incipient_state= State::
//             new_tpx(self.eos.clone(), self.t, self.p, x, DensityInitialization::Guess(old_density))?;            
//             old_density=incipient_state.rho;

//             lnphix=incipient_state.lnphi()?;

//             vw = ( &hy-&lnphix).exp();
            
//             x = &vw/vw.sum();

//             res = (&vw - vw_old).pow2().sum().sqrt();
//             it+=1
//         }
//         // dbg!(&vx);
        
//         let hx = x.ln()+ lnphix;
//         let dg = (vw*(hx-hy)).sum();

//         Ok(dg)
//     }
    
// }

