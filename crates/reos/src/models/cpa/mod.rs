
pub mod associative;
pub mod sites;
pub mod parameters;

use crate::models::cpa::parameters::{ASCParameters, AssociationPureRecord};
use crate::models::cubic::parameters::{CubicParameters, CubicPureRecord};
use crate::models::{cubic::Cubic,cpa::associative::Associative};
use crate::residual::Residual;
use crate::state::eos::EosResult;
use crate::state::State;
use core::f64;
use ndarray::Array1;




#[derive(Clone)]
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

    fn residual_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        Ok(
        self.cubic.residual_helmholtz(t, rho, x)?+self.assoc.residual_helmholtz(t, rho, x)?)
        
    }
}


impl State<CPA> {
    
    pub fn non_bonded_sites(&self)->Array1<f64> {

        let (t,rho,x)=(self.t,self.rho,&self.x);
        self.eos.residual.assoc.X_tan(t, rho, x).unwrap()
    }

}

#[cfg(test)]
mod tests {
    use std::sync::Arc;

    pub use approx::assert_relative_eq;
    
    // use nalgebra::{DMatrix, DVector};
    use ndarray::{array, Array1};

    use crate::{ models::cpa::parameters::{acetic_acid_water, methanol_2b, methanol_3b, water_acetic_acid, water_co2, water_octane_acetic_acid}, state::{S, State, density_solver::DensityInitialization}};

    // use ndarray_linalg::{lapack::solve, solve, Solve};
    // #[test]
    // pub fn test_phi_4c_inert_1a(){
        
    //     println!("---WATER & OCTANE & ACETIC ACID---\n");
    //     let eos = water_octane_acetic_acid().into();
    //     let p=500e5;
    //     let t=298.15;

    //     let x=array![0.2,1e-12,0.8];

    //     let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        
    //     let phi=state.lnphi().unwrap().exp();

    //     let phi_4c_1a_with_inert=array![phi[0],phi[2]];

    //     let phi_4c_1a=get_phi_4c_1a();
    //     dbg!(&phi_4c_1a);
    //     dbg!(&phi_4c_1a_with_inert)


    //     let dif=&phi_4c_1a-phi_4c_1a_with_inert;
    //     let err_norm=dif.mapv(|x|x.powi(2)).sum().sqrt();
        
    //     assert_relative_eq!(err_norm,0.0,epsilon=1e-12);
    // }

    #[test]
    fn test_phi_water_co2(){


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

    pub fn get_phi_4c_1a()->Array1<f64>{
        
        println!("---WATER & ACETIC ACID---\n");
        let eos = water_acetic_acid().into();
        let p=500e5;
        let t=298.15;

        let x=array![0.2,0.8];
        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        // let state=S::new_trx(eos, t, rho, x);
        
        // let p=state.eos.residual.assoc.parameters.clone();

        let phi=state.lnphi().unwrap().exp();

        phi
    }

    pub fn get_phi_1a_4c()->Array1<f64>{
        
        println!("---ACETIC ACID & WATER---\n");
        let eos = acetic_acid_water().into();
        let p=500e5;
        let t=298.15;
        let x=array![0.8,0.2];

        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        // println!{"{}",format!("{}",state)};
        let phi_1a_4c=state.lnphi().unwrap().exp();

        phi_1a_4c


        
    }

    #[test]
    pub fn test_permutation_between_1a_4c(){

        let phi_1a4c=get_phi_1a_4c();
        let phi_4c1a_inv=get_phi_4c_1a().slice(ndarray::s![..;-1]).to_owned();

        let dif=&phi_1a4c-phi_4c1a_inv;
        let err_norm=dif.mapv(|x|x.powi(2)).sum().sqrt();
        assert_relative_eq!(err_norm,0.0,epsilon=1e-12);

    }

    #[test]
    fn test_phi_metoh_2b_3b(){

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
        println!("X2B=\n{}",x_1);
        println!("X3B=\n{}",x_2);
        assert_relative_eq!(2.0*xa_2-1.0,xb_2,epsilon=1e-8);
        assert_relative_eq!(x_1[0],x_1[1],epsilon=1e-8);
        assert_relative_eq!(phi_1[0],cmp_1[0],epsilon=1e-9);
        assert_relative_eq!(phi_2[0],cmp_2[0],epsilon=1e-9);


        let grad1=eos_1.residual.assoc.grad(t, rho_1, &x, &x_1);
        let grad2=eos_2.residual.assoc.grad(t, rho_2, &x, &x_2);

        // println!("dX1={}\n",grad1);
        // println!("dX2={}",grad2);



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


        // #[test]
    pub fn dbg_assoc_p(){

        test_permutation_between_1a_4c(); // also test ECR 
        println!("\n");
        test_phi_metoh_2b_3b();
        println!("\n");

        test_phi_water_co2();
        println!("\n");

        // cargo test dbg_assoc_p -- --nocapture >src/parameters/dbg/dbg_assoc_p.txt

    }
}