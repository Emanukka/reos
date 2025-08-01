use core::f64;
use ndarray::Array1;
use crate::models::associative::Associative;
use crate::models::cubic::Cubic;
use crate::parameters::association::{ASCParameters, AssociationPureRecord};
use crate::parameters::cubic::{CubicParameters, CubicPureRecord};
use crate::residual::Residual;
use crate::state::eos::EosResult;


pub struct CPA{
    pub cubic: Cubic,
    pub assoc:Associative,

}

impl CPA {
    pub fn new(c:Cubic,a:Associative)->Self{
        Self{cubic:c,assoc:a}
    }
    pub fn from_records(c:Vec<CubicPureRecord>,a:Vec<AssociationPureRecord>)->Self{
        
        CPA::new(
            Cubic::new(CubicParameters::from_records(c),super::cubic::CubicModel::SRK), 
            Associative::new(ASCParameters::from_records(a)))
    }
}
impl Residual for CPA {
    
    fn components(&self)->usize {
        self.cubic.components()
    }
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        Ok(
        self.cubic.pressure(t, rho, x)?
        +self.assoc.pressure(t, rho, x)?
        )

    }
    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        Ok(
            self.cubic.residual_chemical_potential(t, rho, x)?
            +self.assoc.residual_chemical_potential(t, rho, x)?
        )
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.cubic.bmix(x)
    }

}


#[cfg(test)]
mod tests {
    use std::sync::Arc;

    pub use approx::assert_relative_eq;
    
    // use nalgebra::{DMatrix, DVector};
    use ndarray::{array, Array1};

    use crate::{ parameters::association::{methanol_2b, methanol_3b, water_acetic_acid, water_co2}, state::{density_solver::DensityInitialization, State}};

    // use ndarray_linalg::{lapack::solve, solve, Solve};

    #[test]
    fn cmp_phi_water_co2_vapor(){


        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);

        let eos =water_co2();
        let s=State::new_tpx(&Arc::new(eos), t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        
        // let rho=1e3;
        // // let D=eos.residual.assoc.calc_delta_mat(t, rho, &x);
        // // let X=eos.residual.assoc.calc_non_assoc_sites_mat(1e3, &x, None, &D);
        let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        
        let rho=s.rho;
        let xasc=s.eos.residual.assoc.X_tan(t, rho, &x).unwrap();
        // let kmat=s.eos.residual.assoc.association_constants(t, rho, &x, s.eos.residual.assoc.g_func(rho, &x));
        // println!("rh");
        println!("rho={}",rho);
        println!("X=\n{}",xasc);
        // println!("K=\n{}",kmat);
        // println!("phi={}",phi);
        let phi_cmp=array![2.14385745e-04, 5.65853284e-01];

        assert_relative_eq!(phi[0],phi_cmp[0],epsilon=1e-8);
        assert_relative_eq!(phi[1],phi_cmp[1],epsilon=1e-8);

    }



    #[test]
    fn cmp_metoh_3b_2b_xassoc(){

        //State Variables
        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![1.]);

        let eos_1 = Arc::new(methanol_2b());
        let eos_2 = Arc::new(methanol_3b());

        let s1=State::new_tpx(&eos_1, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        let s2=State::new_tpx(&eos_2, t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        let phi_1=s1.lnphi().unwrap().exp();
        let phi_2=s2.lnphi().unwrap().exp();

        let cmp_1= [0.0006861214969843528];
        let cmp_2= [0.0007471714606619553];


        let rho_1=s1.rho;
        let rho_2=s2.rho;
        let x_1=eos_1.residual.assoc.X_tan(t, rho_1, &x).unwrap();
        let x_2=eos_2.residual.assoc.X_tan(t, rho_2, &x).unwrap();


        let xa_2=x_2[0];
        let xb_2=x_2[1];
        // println!("cmp_metoh_3b_2b_xassoc");
        // println!("rh");
        println!("X2B=\n{}",eos_1.residual.assoc.X_tan(t, s1.rho, &x).unwrap());
        println!("X3B=\n{}",eos_2.residual.assoc.X_tan(t, s2.rho, &x).unwrap());
        assert_relative_eq!(2.0*xa_2-1.0,xb_2,epsilon=1e-8);
        assert_relative_eq!(x_1[0],x_1[1],epsilon=1e-8);
        assert_relative_eq!(phi_1[0],cmp_1[0],epsilon=1e-9);
        assert_relative_eq!(phi_2[0],cmp_2[0],epsilon=1e-9);

    }

    #[test]
    fn cmp_phi_water_acoh(){


        //State Variables
        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);

        let eos = water_acetic_acid();
        // let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();   //0.09s
        let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();   //

        let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();


        let cmp=array![0.00010095530780761838, 8.66809157609047e-05];
        // dbg!(s.rho);
        assert_relative_eq!(phi[0],cmp[0],epsilon=1e-10);
        assert_relative_eq!(phi[1],cmp[1],epsilon=1e-10);

    }
    #[test]
    fn test_val(){


        //State Variables
        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![1e-6,1.0]);

        let eos = water_acetic_acid();
        // let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();   //0.09s
        let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();   //

        let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();

        let cmp=array![0.00009551080457744488, 0.00007903286580072037];
        // dbg!(s.rho);
        assert_relative_eq!(phi[0],cmp[0],epsilon=1e-10);
        assert_relative_eq!(phi[1],cmp[1],epsilon=1e-10);

    }

}