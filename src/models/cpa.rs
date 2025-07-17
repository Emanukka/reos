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


pub fn water_acetic_acid()->E<CPA>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
            0.0145e-3
        );
        let a2=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
            0.0468e-3
        );


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c1,c2],
            vec![a1,a2]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.0),-0.222 );
        cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
        //Create new State
        E::from_residual(cpa)
} 

pub fn acoh_octane()->E<CPA>{
            //Records
        //1:Acetic Acid, 2:Octane
        let c1=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c2=CubicPureRecord::new(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

        let a1=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
            0.0468e-3
        );
        let a2=AssociationPureRecord::inert(0.1424e-3);


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c1,c2],
            vec![a1,a2]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.0), 0.064 );
        
        //Create new State
        E::from_residual(cpa)
        
} 

mod tests {
    use std::sync::Arc;

    use approx::assert_relative_eq;
    use ndarray::{array, Array1, Array2};

    use crate::{ parameters::{association::{ASCParameters, AssociationPureRecord, AssociationRule}, cubic::CubicParameters, AssociativeTerm, Parameters}, state::{density_solver::DensityInitialization, eos::EquationOfState, State}};


    // #[test]
    // fn cpa_t2(){

    //     let comps = vec!["water","co2"];
    //     let ncomp = comps.len();

    //     let water=AssociationPureRecord::associative(
    //         166.55e2, 
    //         0.0692, 
    //         [2,2,0],
    //         0.0145e-3
    //     );
    //     let co2=AssociationPureRecord::solvate(
    //         [0,1,0],
    //         0.0272e-3
    //     );

    //     let mut pasc=ASCParameters::from_records(vec![water,co2]);

    //     // println!("Antes de setar parametro binário:{}",&pasc);

        
    //     pasc.set_binary_interaction(0, 1, AssociationRule::MCR1, None, Some(0.1836));
    //     // println!("Depois de setar parametro binário:{}",&pasc);

        
    //     let set_1= "konteo_1_mod.json";

    //     let p=500e5;
        
    //     let t=298.15;
    //     let x=Array1::from_vec(vec![0.5,0.5]);
    //     let cubic = CubicParameters::from_file(set_1, comps.clone());

        
    //     let c=Cubic::new(cubic, crate::models::cubic::CubicModel::SRK);
    //     let a =Associative::new(pasc);
    //     let SMAT=&a.parameters.site_multiplicity;
    //     let cpa = CPA::new(c, a);
        
    //     let eos =EquationOfState::from_residual(cpa);
    //     let s=State::new_tpx(&Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();

    //     // let rho=1e3;
    //     // // let D=eos.residual.assoc.calc_delta_mat(t, rho, &x);
    //     // // let X=eos.residual.assoc.calc_non_assoc_sites_mat(1e3, &x, None, &D);
    //     let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        
    //     // println!("X ok:{}",X.unwrap());
    //     // let volume=1./s.rho;

    //     let phi_cmp=array![2.14385745e-04, 5.65853284e-01];

    //     assert_relative_eq!(phi[0],phi_cmp[0],epsilon=1e-8);
    //     assert_relative_eq!(phi[1],phi_cmp[1],epsilon=1e-8);
    //     // dbg!(phi);

    // }
    // #[test]
    // fn cpa_t3(){

    //     let comps = vec!["water","acetic"];
    //     let ncomp = comps.len();

    //     let water=AssociationPureRecord::associative(
    //         166.55e2, 
    //         0.0692, 
    //         [2,2,0],
    //         0.0145e-3
    //     );
    //     let acetic=AssociationPureRecord::associative(
    //         403.23e2, 
    //         4.5e-3, 
    //         [0,0,1],
    //         0.0468e-3
    //     );

    //     let mut pasc=ASCParameters::from_records(vec![water,acetic]);

    //     println!("Antes de setar parametro binário:{}",&pasc);

    //     // pasc.set_binary_interaction(0, 1, AssociationRule::CR1, None, None);
    //     // println!("Depois de setar parametro binário:{}",&pasc);
    //     // println!("Depois de setar parametro binário:{}",&);

        
    //     let set_1= "konteo_1_mod.json";

    //     let p=500e5;
        
    //     let t=298.15;
    //     let x=Array1::from_vec(vec![0.5,0.5]);
    //     let cubic = CubicParameters::from_file(set_1, comps.clone());

    //     println!("Depois de setar parametro binário:{}",&cubic);
        
    //     let c=Cubic::new(cubic, crate::models::cubic::CubicModel::SRK);
    //     let a =Associative::new(pasc);
    //     let cpa = CPA::new(c, a);
        
    //     let eos =Arc::new(EquationOfState::from_residual(cpa));
    //     let s=State::new_tpx(&eos.clone(), t, p, x, DensityInitialization::Vapor).unwrap();

    //     // let rho=1e3;
    //     // let D=eos.residual.assoc.calc_delta_mat(t, rho, &x);
    //     // let X=eos.residual.assoc.calc_non_assoc_sites_mat(1e3, &x, None, &D);
    //     let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        
    //     // println!("{}",phi);
    //     // // println!("X ok:{}",X.unwrap());
    //     // // let volume=1./s.rho;

    //     let phi_cmp=array![9.884693365193102e-05, 7.676491784564717e-05];

    //     assert_relative_eq!(phi[0],phi_cmp[0],epsilon=1e-10);
    //     assert_relative_eq!(phi[1],phi_cmp[1],epsilon=1e-10);

    // }
    // #[test]
    // fn cpa_t4(){

    //     let comps = vec!["water","acetic"];
    //     let ncomp = comps.len();

    //     let water=AssociationPureRecord::associative(
    //         166.55e2, 
    //         0.0692, 
    //         [2,2,0],
    //         0.0145e-3
    //     );
    //     let acetic=AssociationPureRecord::associative(
    //         403.23e2, 
    //         4.5e-3, 
    //         [0,0,1],
    //         0.0468e-3
    //     );

    //     let pasc=ASCParameters::from_records(vec![water,acetic]);
 
    //     // println!("Depois de setar parametro binário:{}",&pasc);
        
    //     let set_1= "konteo_1_mod.json";

    //     let p=500e5;
        
    //     let t=298.15;
    //     let x=Array1::from_vec(vec![0.5,0.5]);
    //     let cubic = CubicParameters::from_file(set_1, comps.clone());

    //     // println!("Depois de setar parametro binário:{}",&cubic);
        
    //     let c=Cubic::new(cubic, crate::models::cubic::CubicModel::SRK);
    //     let a =Associative::new(pasc);
    //     let cpa = CPA::new(c, a);
        
    //     let eos =Arc::new(EquationOfState::from_residual(cpa));
    //     let s=State::new_tpx(&eos, t, p, x, DensityInitialization::Vapor).unwrap();

    //     // let rho=1e3;
    //     // let D=eos.residual.assoc.calc_delta_mat(t, rho, &x);
    //     // let X=eos.residual.assoc.calc_non_assoc_sites_mat(1e3, &x, None, &D);
    //     let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
        

    //     let phi_cmp=array![9.884693365193102e-05, 7.676491784564717e-05];

    //     assert_relative_eq!(phi[0],phi_cmp[0],epsilon=1e-10);
    //     assert_relative_eq!(phi[1],phi_cmp[1],epsilon=1e-10);

    // }
    // #[test]
    // fn cpa_t5(){


    //     //State Variables
    //     let p=500e5;
    //     let t=298.15;
    //     let x=Array1::from_vec(vec![0.5,0.5]);

    //     let eos = water_acetic_acid();
    //     // let s=State::new_tpx(Arc::new(eos), t, p, x, DensityInitialization::Vapor).unwrap();

    //     // let phi: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>>=s.lnphi().unwrap().exp();
    //     // let cmp=array![0.00010095530780761838, 8.66809157609047e-05];
    //     // dbg!(s.rho);
    //     // assert_relative_eq!(phi[0],cmp[0],epsilon=1e-10);
    //     // assert_relative_eq!(phi[1],cmp[1],epsilon=1e-10);

    // }
    // #[test]
    // fn cpa_t6(){


    //     //State Variables
    //     let p=500e5;
    //     let t=298.15;
    //     let x=Array1::from_vec(vec![0.5,0.5]);

    //     let eos = Arc::new(acoh_octane());
    //     let s=State::new_tpx(&eos, t, p, x, DensityInitialization::Vapor).unwrap();
    //     s.lnphi();

    // }


}