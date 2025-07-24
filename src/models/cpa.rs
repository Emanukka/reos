use core::f64;
use ndarray::Array1;
use crate::models::associative::Associative;
use crate::models::cubic::Cubic;
use crate::parameters::association::{ASCParameters, AssociationPureRecord, AssociationRule};
use crate::parameters::cubic::{CubicParameters, CubicPureRecord};
use crate::residual::Residual;
use crate::state::eos::EosResult;
use crate::state::E;


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



mod tests {
    use std::sync::Arc;

    use approx::assert_relative_eq;
    use nalgebra::{DMatrix, DVector};
    use ndarray::{array, Array1, Array2};

    use crate::{ parameters::{association::{methanol_2b, methanol_3b, water_acetic_acid, water_co2, ASCParameters, AssociationPureRecord, AssociationRule}, cubic::CubicParameters, AssociativeTerm, Parameters}, state::{density_solver::DensityInitialization, eos::EquationOfState, State}};

    // use ndarray_linalg::{lapack::solve, solve, Solve};

    #[test]
    fn cpa_t2(){

        let comps = vec!["water","co2"];
        let ncomp = comps.len();

        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);

        let eos =water_co2();
        let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();

        // let rho=1e3;
        // // let D=eos.residual.assoc.calc_delta_mat(t, rho, &x);
        // // let X=eos.residual.assoc.calc_non_assoc_sites_mat(1e3, &x, None, &D);
        let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        
        // println!("X ok:{}",X.unwrap());
        // let volume=1./s.rho;

        let phi_cmp=array![2.14385745e-04, 5.65853284e-01];

        assert_relative_eq!(phi[0],phi_cmp[0],epsilon=1e-8);
        assert_relative_eq!(phi[1],phi_cmp[1],epsilon=1e-8);
        // dbg!(phi);

    }

    #[test]
    fn cpa_t3(){


        //State Variables
        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);

        let eos = water_acetic_acid();
        let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();

        let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        let cmp=array![0.00010095530780761838, 8.66809157609047e-05];
        // dbg!(s.rho);
        assert_relative_eq!(phi[0],cmp[0],epsilon=1e-10);
        assert_relative_eq!(phi[1],cmp[1],epsilon=1e-10);

    }

    #[test]
    fn cpa_t5(){

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


        let xa_2=x_2[(0,0)];
        let xb_2=x_2[(1,0)];

        assert_relative_eq!(2.0*xa_2-1.0,xb_2,epsilon=1e-8);
        assert_relative_eq!(x_1[(0,0)],x_1[(1,0)],epsilon=1e-8);
        assert_relative_eq!(phi_1[0],cmp_1[0],epsilon=1e-9);
        assert_relative_eq!(phi_2[0],cmp_2[0],epsilon=1e-9);

    }

    #[test]
    fn cpa_t6(){


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
    fn cpa_t7(){


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
    #[test]
    fn mat1(){


        // let a=Array2::from_shape_vec((3,2), vec![1,2,3,4,5,6]).unwrap();
        // println!("A={}",a);
        // let a_flat=a.flatten_with_order(ndarray::Order::ColumnMajor);
        // println!("A vetor [ordem Fortran (coluna)]={}",a_flat);
        let a: Array2<f64>=Array2::from_shape_vec((2,2), vec![1.,2.,3.,4.]).unwrap();
        let b: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=array![1.0,2.];


    }
//     fn convert_ndarray_to_nalgebra(nd: &Array2<f64>) -> DMatrix<f64> {
//     let (nrows, ncols) = nd.dim();
//     DMatrix::from_row_slice(nrows, ncols, nd.as_slice().unwrap())
// }
    // #[test]
    // fn mat2(){


    //     let a=Array2::from_shape_vec((3,2), vec![1.,2.,3.,4.,5.,6.]).unwrap();

    //     // let b=convert_ndarray_to_nalgebra(&a);

    //     println!("A={}",a);
    //     println!("a={}",b);
    //     // let a_flat=a.flatten_with_order(ndarray::Order::ColumnMajor);
    //     // println!("A vetor [ordem Fortran (coluna)]={}",a_flat);

    // }



}