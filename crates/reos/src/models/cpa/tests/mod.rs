use crate::{models::{associative::parameters::AssociationPureRecord, cpa::{SCPAsrkCR1, SCPAsrkECR, parameters::{CPAParameters, CPAPureRecord}}, cubic::parameters::CubicPureRecord}, parameters::Parameters, state::E};


pub fn water_acetic_acid()->E<SCPAsrkECR>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
        );
        let a2=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
        );


        let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
        let mut parameters = CPAParameters::from_records(records);
        parameters.cubic.set_kij(0, 1, -0.222);
        let cpa = SCPAsrkECR::from_parameters(parameters);
        E::from_residual(cpa)

} 
// pub fn water_octane_acetic_acid()->EOS{
//             //Records
//         //1:Water, 2:Acetic Acid
//         let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
//         let c2=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);
//         let c3=CubicPureRecord::new_set1(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

//         let a1=AssociationPureRecord::associative(
//             166.55e2, 
//             0.0692, 
//             [2,2,0],
//             0.0145e-3
//         );
//         let a2=AssociationPureRecord::associative(
//             403.23e2, 
//             4.5e-3, 
//             [0,0,1],
//             0.0468e-3
//         );

//         let a3=AssociationPureRecord::inert(0.1424e-3);


//         //CPA eos
//         // let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1,c2,c3]));
//         // let a = Associative::from_parameters(ASCParameters::<Konteogeorgis>::from_records(vec![a1,a2,a3]));
//         // let mut cpa = CPA::from_models(c, a);

//         //Set binary parameters 
//         cpa.cubic.parameters.set_kij(0, 1, -0.222);
//         cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
//         //Create new State
//         E::from_residual(cpa)

//                 // cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
//         let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2),CPAPureRecord::new(c3, a3)];
//         let mut parameters = CPAParameters::from_records(records);
//         parameters.cubic.set_kij(0, 1, -0.222);
//         let cpa = SCPAsrkECR::from_parameters(parameters);
//         //Create new State
//         E::from_residual(cpa)

// } 
pub fn acetic_acid_water()->E<SCPAsrkECR>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
        );
        let a2=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
        );

        let records = vec![CPAPureRecord::new(c2, a2),CPAPureRecord::new(c1, a1)];
        let mut parameters = CPAParameters::from_records(records);
        parameters.cubic.set_kij(0, 1, -0.222);
        let cpa = SCPAsrkECR::from_parameters(parameters);
        E::from_residual(cpa)

} 

pub fn water_co2()->E<SCPAsrkCR1>{
            //Records
        let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
        );
        let a2=AssociationPureRecord::solvate(
            [0,1,0]);

        let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
        let mut parameters = CPAParameters::from_records(records);
        // parameters.cubic.set_kij(0, 1, -0.222);

        parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
        parameters.assoc.set_binary_from_owners(0, 1, None, Some(0.1836));

        let cpa = SCPAsrkCR1::from_parameters(parameters);
        //Create new State
        E::from_residual(cpa)
} 
pub fn co2_water()->E<SCPAsrkCR1>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
        );
        let a2=AssociationPureRecord::solvate(
            [0,1,0],
        );

        let records = vec![CPAPureRecord::new(c2, a2),CPAPureRecord::new(c1, a1)];
        let mut parameters = CPAParameters::from_records(records);

        parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
        parameters.assoc.set_binary_from_owners(0, 1, None, Some(0.1836));

        let cpa = SCPAsrkCR1::from_parameters(parameters);
        //Create new State
        E::from_residual(cpa)

} 


pub fn methanol_2b()->E<SCPAsrkCR1>{
            //Records
        //1:metoh, 2:oct
        let c1=CubicPureRecord::new_set1(0.40531, 0.0000309, 0.4310, 513.);

        let a1=AssociationPureRecord::associative(
            24591.0, 
            0.01610, 
            [1,1,0],
        );
        let records = vec![CPAPureRecord::new(c1, a1)];
        let parameters = CPAParameters::from_records(records);


        let cpa = SCPAsrkCR1::from_parameters(parameters);
        //Create new State
        E::from_residual(cpa)

        
} 
pub fn methanol_3b()->E<SCPAsrkCR1>{
            //Records
        //1:metoh, 2:oct
        let c1=
        CubicPureRecord::new_set1(
            4.5897e-1, 
            0.0334e-3, 
            1.0068,
            513.);

        let a1=AssociationPureRecord::associative(
            160.70e2, 
            34.4e-3, 
            [2,1,0],
        );

        let records = vec![CPAPureRecord::new(c1, a1)];
        let parameters = CPAParameters::from_records(records);


        let cpa = SCPAsrkCR1::from_parameters(parameters);
        E::from_residual(cpa)
} 


pub fn acoh_octane()->E<SCPAsrkCR1>{
            //Records
        //1:Acetic Acid, 2:Octane
        let c1=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c2=CubicPureRecord::new_set1(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

        let a1=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
        );
        let a2=AssociationPureRecord::inert();


        let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
        let parameters = CPAParameters::from_records(records);


        let mut cpa = SCPAsrkCR1::from_parameters(parameters);
        cpa.cubic.parameters.set_kij(0, 1, 0.064);

        E::from_residual(cpa)

        
} 
pub fn octane_acoh()->E<SCPAsrkCR1>{
            //Records
        //1:Acetic Acid, 2:Octane
        let c1=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c2=CubicPureRecord::new_set1(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

        let a1=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
        );
        let a2=AssociationPureRecord::inert();


        let records = vec![CPAPureRecord::new(c2, a2),CPAPureRecord::new(c1, a1)];
        let parameters = CPAParameters::from_records(records);


        let mut cpa = SCPAsrkCR1::from_parameters(parameters);

        cpa.cubic.parameters.set_kij(0, 1, 0.064);

        E::from_residual(cpa)

        
} 


#[cfg(test)]
mod tests {
    use std::sync::Arc;

    pub use approx::assert_relative_eq;
    use ndarray::{Array, Array1, Dimension, array};

    use crate::{models::{associative::parameters::CombiningRule, cpa::{association::Associative, rdf::RdfModel, tests::water_co2}}, state::{S, State, density_solver::DensityInitialization}};
    use super::*;
    // use crate::{ models::{IDEAL_GAS_CONST, associative::parameters::CombiningRule, cpa::{CR1, parameters::{acetic_acid_water, methanol_2b, methanol_3b, water_acetic_acid, water_co2}}}, state::{S, State, density_solver::DensityInitialization}};

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
    // macro_rules! arr_eq {
    //     ($a:expr, $b:expr, $tol:expr) => {{
    //         if $a.len() != $b.len() {
    //             false
    //         } else {
    //             $a.iter().zip($b.iter()).all(|(x, y)| {
    //               let err = (x - y).abs();
    //               if err <= $tol{
    //                 true
    //               }
    //               else {
    //                 println!("left = {x}, right = {y}");
    //                 false
    //               }
    //             })
    //           }}
    //       }
    // }
    
    fn arr_eq<D:Dimension>(a:Array<f64, D>,b:Array<f64,D>,tol:f64)->bool{

      if a.shape() != b.shape() {
          return false;
      }
      a.iter()
      .zip(b.iter())
      .all(|(x,y)|{
        let err = (x-y).abs();
        if err <= tol{
          true
        }
        else {
          println!("left = {x}, right = {y}");
          false
        }
      })
    }

    fn dbg_assoc<T:CombiningRule,R:RdfModel>(
      assoc:&Associative<T,R>,
      t:f64,
      rho:f64,
      x:&Array1<f64>,
      p_ok:Array1<f64>,
      k_ok:Array1<f64>,
      l_ok:Array1<f64>,
      g_ok:Array1<f64>,
      x_ok:Array1<f64>,
      tol:f64,
      )->bool{
        
      let m = assoc.m(x);

      let lambda = assoc.parameters.lambda.flatten().to_owned();
      let gamma = assoc.parameters.gamma.flatten().to_owned();

      let p      = assoc.parameters.p.flatten().to_owned();
      let k      = assoc.k(t, rho, x, &m).flatten().to_owned();
      
      let x_tan = assoc.x_tan(t, rho, x).unwrap();
      let ok_p = arr_eq(p, p_ok, tol);
      let ok_k = arr_eq(k, k_ok, tol);
      let ok_l = arr_eq(lambda, l_ok, tol);
      let ok_g = arr_eq(gamma, g_ok, tol);
      let ok_x = arr_eq(x_tan, x_ok, tol);

      assert!(ok_p, "matrix P error");
      assert!(ok_k, "matrix K error");
      assert!(ok_l, "matrix lambda error");
      assert!(ok_g, "matrix gamma error");
      assert!(ok_x, "x_tan error");

      true
    }


    #[test]
    fn dbg_water_co2_matrices() {
      let rho= 55_000.0;
      let t=298.15;
      let x= &Array1::from_vec(vec![0.5,0.5]);
      let assoc = &water_co2().residual.assoc;
      
      let p_ok = array![0., 1., 1., 1., 0., 0., 1., 0., 0.];
      let k_ok = array![ 0. , 25.04896564,  3.21025658, 25.04896564,  0., 0., 3.21025658,  0.,  0.];
      let l_ok = array![1., 1., 0., 0., 0., 1.];
      let g_ok = array![1., 0., 0., 1.];
      let x_ok = array![0.03874146, 0.20484502, 0.66778819];
      let tol = 1e-5;
      let ok = dbg_assoc(assoc, t, rho, x, p_ok, k_ok, l_ok, g_ok, x_ok, tol);

      assert!(ok);

    }

    #[test]
    fn test_phi_water_co2(){

        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);
        let eos = water_co2();
        let s=State::new_tpx(&Arc::new(eos), t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        let phi=s.lnphi().unwrap().exp();
        let phi_cmp=array![2.14385745e-04, 5.65853284e-01];
        
        assert_relative_eq!(phi[0],phi_cmp[0],epsilon=1e-8);
        assert_relative_eq!(phi[1],phi_cmp[1],epsilon=1e-8);

    }

    pub fn get_phi_4c_1a()->Array1<f64>{
        
        let eos = water_acetic_acid().into();
        let p=500e5;
        let t=298.15;

        let x=array![0.2,0.8];
        let state = S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        let phi=state.lnphi().unwrap().exp();

        phi
    }

    pub fn get_phi_1a_4c()->Array1<f64>{
        
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
    pub fn permutation_between_water_acoh(){

        let phi_1a4c= get_phi_1a_4c();
        let phi_4c1a_inv= get_phi_4c_1a().slice(ndarray::s![..;-1]).to_owned();

        let dif=&phi_1a4c-phi_4c1a_inv;
        let err_norm=dif.mapv(|x|x.powi(2)).sum().sqrt();
        assert_relative_eq!(err_norm,0.0,epsilon=1e-12);

    }

    #[test]
    fn phi_metoh_2b_3b(){

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
        let x_1=eos_1.residual.assoc.x_tan(t, rho_1, &x).unwrap();
        let x_2=eos_2.residual.assoc.x_tan(t, rho_2, &x).unwrap();


        let xa_2=x_2[0];
        let xb_2=x_2[1];
        // println!("cmp_metoh_3b_2b_xassoc");
        // println!("rh");
        // println!("X2B=\n{}",x_1);
        // println!("X3B=\n{}",x_2);
        assert_relative_eq!(2.0*xa_2-1.0,xb_2,epsilon=1e-8);
        assert_relative_eq!(x_1[0],x_1[1],epsilon=1e-8);
        assert_relative_eq!(phi_1[0],cmp_1[0],epsilon=1e-9);
        assert_relative_eq!(phi_2[0],cmp_2[0],epsilon=1e-9);


        // let grad1=eos_1.residual.assoc.grad(t, rho_1, &x, &x_1);
        // let grad2=eos_2.residual.assoc.grad(t, rho_2, &x, &x_2);

        // println!("dX1={}\n",grad1);
        // println!("dX2={}",grad2);



    }

    #[test]
    fn phi_water_acoh(){


        //State Variables
        let p=500e5;
        let t=298.15;
        let x=Array1::from_vec(vec![0.5,0.5]);

        let eos = water_acetic_acid();
        let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();   //

        let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        let cmp=array![0.00010095530780761838, 8.66809157609047e-05];
        assert_relative_eq!(phi[0],cmp[0],epsilon=1e-10);
        assert_relative_eq!(phi[1],cmp[1],epsilon=1e-10);

    }
}