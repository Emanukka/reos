use std::f64::NAN;

use ndarray::Array1;

use crate::tools::newton;
mod stability;

use super::{cpa::EquationOfState, density_solver::{density, DensityInitialization}, equation_of_state::EosError, state::State};

const XW: f64 = 0.999999;
const YWGUESS:f64 = 2000e-6;

pub struct WaterSaturarion<E:EquationOfState>{

    saturation:f64,
    state:State<E>,
}

impl<E: EquationOfState> std::fmt::Display for WaterSaturarion<E> {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let state = &self.state;

        writeln!(f, "WSAT Calculation:")?;
        writeln!(f, "  Temperature       = {:>10.5} K", state.t)?;
        writeln!(f, "  Pressure          = {:>10.5} bar", state.p / 1e5)?;
        writeln!(f, "  Water Saturation  = {:>10.6} ppm", self.saturation * 1e6)?;
    
        Ok(())

    }
}

impl <E:EquationOfState> WaterSaturarion<E>{

    pub fn saturation(dry_gas:State<E>)->Result<Self,EosError>{

        let ncomp = dry_gas.x.len();
        let water_index = ncomp-1;
        let y_dry_gas=dry_gas.x.clone();

        let eos = dry_gas.eos.clone();

        let t=dry_gas.t;
        let p=dry_gas.p;

        let res = |log10_yw0:f64|{
    
            let yw0 = 10.0_f64.powf(log10_yw0);
            
            let mut y_wet_0 = (&y_dry_gas/y_dry_gas.sum())*(1.-yw0);

            y_wet_0[water_index] = yw0;
    
            let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;

            vx[water_index]= XW;

            let vapor = State::new_tpx(eos.clone(), t, p, y_wet_0, DensityInitialization::Vapor).unwrap(); 

            let liquid = State::new_tpx(eos.clone(), t, p, vx, DensityInitialization::Liquid).unwrap(); 

            let lnphi_water_vap=vapor.lnphi().unwrap_or(Array1::from_elem(ncomp, NAN))[water_index];
            let lnphi_water_liq = liquid.lnphi().unwrap_or(Array1::from_elem(ncomp, NAN))[water_index];
            

            yw0.ln()+lnphi_water_vap-lnphi_water_liq
    
        };
        
        // testar com e sem chute de antoine
        // modificar pra receber chute de antoine
        let log10yw0 = YWGUESS.log10();
        let result = newton(res, log10yw0, None, None);
        

        if let Ok(res)=result{

            let yguess =10.0_f64.powf(res.x);
            // dbg!(yguess);

            // 
            // tpd_root()
            // att. valor de y=> modifica o state_y
            let log10_ywguess_from_phiphi = yguess.log10();
            // let ncomp = eos.asc.parameters.ncomp;
            // let water_index = ncomp-1;



            // residual function, used inside newton root find
            let res = |log10_yw0:f64|{
                
                let yw0 = 10.0_f64.powf(log10_yw0);

                let mut vy_wet = (&y_dry_gas/ y_dry_gas.sum())*(1.-yw0);

                vy_wet[water_index] = yw0;

                let vap_state= 
                State::new_tpx(eos.clone(), t, p, vy_wet, DensityInitialization::Vapor).unwrap();

                // dbg!(&vap_state);

                let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;
                vx[water_index]= XW;
    
                let tpd_result= vap_state.tpd(DensityInitialization::Liquid, vx).unwrap();

                // dbg!(tpd_result);
                tpd_result
        
            };

            let newton_result = newton(res, log10_ywguess_from_phiphi, None, None);


            match newton_result{
                
                Ok(ok)=>{
                    let sat= 10.0_f64.powf(ok.x);
                    let y = (&y_dry_gas/ y_dry_gas.sum())*(1.-sat);

                    Ok(
                    Self{
                        saturation:sat,
                        state: State::new_tpx(eos, t, p, y, DensityInitialization::Vapor)?
                    }
                    )

                }
                Err(_)
                =>{
                    Err(EosError::NotConverged("tpd root find".to_string()))
                }
            }            

        }else {
            
            return Err(EosError::NotConverged("phi-phi guess root find".to_string()))

        }

    }

}


impl<E:EquationOfState> State<E> {

    pub fn water_saturation(&self)->Result<f64,EosError>{


        let ncomp = self.x.len();
        let water_index = self.x.len()-1;

        let mut wet = self.x.clone();

        wet[water_index]=0.0;

        let y_dry_gas = (1.0/(wet.sum()))*wet;

        dbg!(&y_dry_gas);

        let eos = self.eos.clone();

        let t=self.t;
        let p=self.p;

        let res = |log10_yw0:f64|{
    
            let yw0 = 10.0_f64.powf(log10_yw0);
            
            let mut y_wet_0 = (&y_dry_gas/y_dry_gas.sum())*(1.-yw0);

            y_wet_0[water_index] = yw0;
    
            let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;

            vx[water_index]= XW;

            let vapor = State::new_tpx(eos.clone(), t, p, y_wet_0, DensityInitialization::Vapor).unwrap(); 

            let liquid = State::new_tpx(eos.clone(), t, p, vx, DensityInitialization::Liquid).unwrap(); 

            let lnphi_water_vap=vapor.lnphi().unwrap_or(Array1::from_elem(ncomp, NAN))[water_index];
            let lnphi_water_liq = liquid.lnphi().unwrap_or(Array1::from_elem(ncomp, NAN))[water_index];
            

            yw0.ln()+lnphi_water_vap-lnphi_water_liq
    
        };
        
        // testar com e sem chute de antoine
        // modificar pra receber chute de antoine
        let log10yw0 = YWGUESS.log10();
        let result = newton(res, log10yw0, None, None);
        

        if let Ok(res)=result{

            let yguess =10.0_f64.powf(res.x);
            // dbg!(yguess);

            // 
            // tpd_root()
            // att. valor de y=> modifica o state_y
            let log10_ywguess_from_phiphi = yguess.log10();
            // let ncomp = eos.asc.parameters.ncomp;
            // let water_index = ncomp-1;



            // residual function, used inside newton root find
            let res = |log10_yw0:f64|{
                
                let yw0 = 10.0_f64.powf(log10_yw0);

                let mut vy_wet = (&y_dry_gas/ y_dry_gas.sum())*(1.-yw0);

                vy_wet[water_index] = yw0;

                let vap_state= 
                State::new_tpx(eos.clone(), t, p, vy_wet, DensityInitialization::Vapor).unwrap();

                // dbg!(&vap_state);

                let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;
                vx[water_index]= XW;
    
                let tpd_result= vap_state.tpd(DensityInitialization::Liquid, vx).unwrap();

                // dbg!(tpd_result);
                tpd_result
        
            };

            let newton_result = newton(res, log10_ywguess_from_phiphi, None, None);


            match newton_result{
                
                Ok(ok)=>{
                    let sat= 10.0_f64.powf(ok.x);
                    // let y = (&y_dry_gas/ y_dry_gas.sum())*(1.-sat);

                    Ok(sat)

                }
                Err(_)
                =>{
                    Err(EosError::NotConverged("tpd root find".to_string()))
                }
            }            

        }else {
            
            return Err(EosError::NotConverged("phi-phi guess root find".to_string()))

        }

    }
    
}

