use std::{path::Components, sync::Arc};

use reeos::{models::cpa::CPA, parameters::{association::{AssociationPureRecord, AssociationRule}, cubic::CubicPureRecord, CompRecord}, phase_equilibrium::PhaseEquilibrium, state::{density_solver::DensityInitialization, eos::{EosError, EquationOfState}, State, E}, tools::newton, Array1};

pub mod data;
const XW: f64 = 0.999999;
const YWGUESS:f64 = 2000e-6;

pub struct WaterSaturation{

    pub eos:Arc<E<CPA>>
}

impl WaterSaturation {
    
    fn new(eos:Arc<E<CPA>>)->Self{
        Self{eos}
    }
}

    
impl WaterSaturation {

pub fn phi_phi_guess(&self,t:f64,p:f64,y_dry_gas:&Array1<f64>,ncomp:usize)->Result<f64,EosError>{

    let water_index = ncomp-1;
    let res = |log10_yw0:f64|{



        let yw0 = 10.0_f64.powf(log10_yw0);
        

        let mut vy_wet = (y_dry_gas/y_dry_gas.sum())*(1.-yw0);
        vy_wet[water_index] = yw0;

        let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;
        vx[water_index]= XW;


        // fase vapor (mistura)
        let vapor= State::new_tpx(&self.eos, t, p, vy_wet.clone(), DensityInitialization::Vapor).unwrap();
        // let rho_vap = 1./eos.volume(t, None,p, &vy_wet, "Vapor").x; 
        let phi_water_vap = vapor.lnphi().unwrap().exp()[water_index];


        // fase liq (~agua pura)
        let liquid= State::new_tpx(&self.eos, t, p, vx.clone(), DensityInitialization::Liquid).unwrap();

        let phi_water_liq =liquid.lnphi().unwrap().exp()[water_index];

        ((yw0*phi_water_vap)/(phi_water_liq)).ln()
        // yw0.ln() + phi_water_vap.ln() - phi_water_liq.ln()

    };
    
    let log10yw0 = YWGUESS.log10();

    let result = newton(res, log10yw0, None, None);
    // PhiPhiResult{x:10.0_f64.powf(x.x)}

    match result {

        Ok(x)=> {
            Ok(10.0_f64.powf(x.x))
        },
        Err(e)=> Err(EosError::NotConverged("newton rhapson at phi_phi guess didnt converge".to_string()))
        
    } 
    }

pub fn watcon(&self,t:f64,p:f64,y_dry_gas:&Array1<f64>)->Result<f64,EosError>{
    
        let ncomp=y_dry_gas.len();
        let phi_phi_value=self.phi_phi_guess( t, p, y_dry_gas,ncomp)?;

        let log10_ywguess_from_phiphi = phi_phi_value.log10();
        let water_index = ncomp-1;
        let peq=PhaseEquilibrium::new(self.eos.clone(), None);
        // residual function, used inside newton root find
        let res = |log10_yw0:f64|{
            
            let yw0 = 10.0_f64.powf(log10_yw0);

            // dbg!(yw0);
            let mut vy_wet = (y_dry_gas/ y_dry_gas.sum())*(1.-yw0);
            vy_wet[water_index] = yw0;
            // dbg!(yw0);

            let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;
            vx[water_index]= XW;

            // let tpd_result=tpd(eos,t,p,XW,vy_wet,None,None);
            let tpd_result=peq.tpd(
                t, 
                p, 
                vy_wet.clone(),
                DensityInitialization::Liquid,
                vx.clone(), 
                None, 
                None).unwrap();

            // dbg!(tpd_result);
            tpd_result

        };
        // println!("Newton Watcon");
        // println!("res={log10_ywguess_from_phiphi},watcon");
        let newton_result = newton(res, log10_ywguess_from_phiphi, None, None);

        match newton_result {

            Ok(x)=> {Ok(10.0_f64.powf(x.x))}

            Err(e)=>{

                return Err(EosError::NotConverged("error at newton-watcon".to_string()));
                // println!("ERRO:{e}");
                }

            }
    }
}



pub fn co2_water()->E<CPA>{
        //1:Water, 2:CO2
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new(0.35079, 0.0272e-3, 0.7602, 304.12);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
            0.0145e-3
        );
        let a2=AssociationPureRecord::solvate(
            [0,1,0],
            0.0272e-3
        );


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c2,c1],
            vec![a2,a1]);

        cpa.cubic.set_binary(0,1, Some(0.000877),-0.15508 );
        cpa.assoc.set_binary(0, 1, AssociationRule::MCR1, None, Some(0.1836));
        
        //Create new State
        E::from_residual(cpa)
}
mod tests{
    use std::sync::Arc;

    use reeos::Array1;

    use crate::{co2_water, WaterSaturation};

    
    #[test]
    fn t1(){

        let eos=Arc::new(co2_water());

        let water_sat=WaterSaturation::new(eos.clone());

        let t=298.15;
        let p=500e5;
        let y_dry_gas=Array1::<f64>::from_vec(vec![1.0,0.0]);

        let yw=water_sat.watcon(t, p, &y_dry_gas).unwrap();
        
        // println!("P= {} yw= {}")
    }

}