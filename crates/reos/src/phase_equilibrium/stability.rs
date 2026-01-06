use std::{iter::Zip, sync::Arc};

use ndarray::{ Array1};
use serde::de;

use crate::{phase_equilibrium::PhaseEquilibrium, residual::Residual, state::{E, State, density_solver::DensityInitialization, eos::{EosError, EosResult}}, tools::newton};


impl <R:Residual> State<R> {

    pub fn tpd(&self, other: &State<R>) -> f64{
        
        let hy = self.lnphi() + self.x.ln();
        let lnx= other.x.clone().ln();
        let hx = other.lnphi() + &lnx;

        (lnx.exp()*(hx-hy)).sum()
    }

    pub fn min_tpd(
        &self,
        xphase:DensityInitialization,
        x:Array1<f64>,
        tol:Option<f64>,
        it_max:Option<i32>) -> EosResult<(f64,State<R>)>{

        let tol = tol.unwrap_or(1e-10);
        let max= it_max.unwrap_or(200);

        let ln_phiy = self.lnphi();
        let z = &self.x;
        let hy = Array1::from_shape_fn(z.len(), |i| ln_phiy[i] + z[i].ln() );

        let guess = State::new_tpx(Arc::clone(&self.eos), self.t, self.p, x.clone(), xphase)?;
       
        let mut x = x;
        let mut dens = guess.d;

        for it in 0..max {

            if self.min_tpd_step(&mut dens, &mut x, &hy, tol)? {
                
                break;
            } 

            if it == max - 1 {

                return Err(EosError::NotConverged("TPD min".to_string()))
            }
        }

        let x_phase = State::new_trx(Arc::clone(&self.eos), self.t, dens, x.clone());
        let ln_phix = x_phase.lnphi(); 
        let mut dg = 0.0;

        for i in 0..x.len() {
            
            let hx = ln_phix[i] + x[i].ln();
            dg += x[i] * (hx - hy[i])
        }

        Ok((dg, x_phase))

    }
    
    fn min_tpd_step(
        &self,
        dens: &mut f64,
        x:    &mut Array1<f64>, 
        hy:   &Array1<f64>, 
        tol:  f64
        )->Result<bool, EosError> {

            let n = x.len();
            let s = State::new_tpx(Arc::clone(&self.eos), self.t, self.p, x.clone(), DensityInitialization::Guess(*dens) )?;

            *dens = s.d;
            
            let ln_phi = s.lnphi();
            let ww = Array1::from_shape_fn(n, |i| (hy[i] - ln_phi[i]).exp());
            
            let sum = ww.sum();

            let mut err = 0.0;

            x.iter_mut().zip(&ww).for_each(|(x,&w)| {

                let new = w / sum;
                err += (*x - new).powi(2);
                *x = new;
            });
            
            Ok( err.sqrt() < tol)

    }
}



pub enum TemperatureOrPressure {
    Temperature(f64),
    Pressure(f64)
}
impl<R:Residual> PhaseEquilibrium<R>{



}
#[cfg(test)]
pub mod tests{
    use std::sync::Arc;

    use approx::assert_relative_eq;
    use ndarray::{array, Array1};
 
    use crate::{ models::cpa::{SCPA, parameters::readyto::{acetic1a, water4c, water4c_acetic1a}}, phase_equilibrium::PhaseEquilibrium, state::{E, State, density_solver::DensityInitialization::{Liquid, Vapor}}};

    fn water_acetic_acid() -> E<SCPA>{
        
        let w = water4c();
        let a = acetic1a();
        let b = water4c_acetic1a();

        let cpa = SCPA::from_records(vec![w, a], vec![b]);

        E::from_residual(cpa)


    }
    #[test]
    fn verify_tpd_close_to_bbpoint(){

        let t=313.15;
        let x=Array1::from_vec(vec![0.5,0.5]);

        let eos = water_acetic_acid();
        let eos=Arc::new(eos);
        let peq=PhaseEquilibrium::new(&eos, None);

        
        let (pb,_) = peq.bbpy(t, x.clone(),Some(1e-10),Some(1e-10)).unwrap();
        

        let zphase=State::new_tpx(eos.clone(), t, pb, x.clone(), Liquid).unwrap();
        let (dg_zero, _) = zphase.min_tpd(Vapor, array![0.5,0.5], None, None).unwrap();

        assert_relative_eq!(dg_zero, 0.0, epsilon = 1e-10);

        let zphase = State::new_tpx(eos.clone(), t, pb * 1.2, x.clone(), Liquid).unwrap();
        let (dg, _) = zphase.min_tpd(Vapor, array![0.5,0.5], None, None).unwrap();
        
        assert!( dg > 0.0);
       
        let zphase = State::new_tpx(eos.clone(), t, pb * 0.8, x.clone(), Liquid).unwrap();
        let (dg, _) = zphase.min_tpd(Vapor, array![0.5,0.5], None, None).unwrap();
       
        assert!( dg < 0.0);


    }

}

