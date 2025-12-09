
use crate::models::associative::Associative;
use crate::models::associative::parameters::{ASCParameters, AssociationPureRecord, AssociationRule, AssociationType};
use crate::models::cpa::CPA;
use crate::models::cpa::rdf::{CarnahanStarlingRDF, ElliotRDF, RDF};

use crate::models::cubic::{Cubic, CubicModel, SRK};
use crate::models::cubic::parameters::{CubicParameters, CubicPureRecord};
use crate::parameters::{Parameters, PureRecord};
use crate::state::E;




// pub type CPAPureRecord = (CubicPureRecord, AssociationPureRecord);

pub struct CPAPureRecord{
    //Cubic
    a0:f64,
    b:f64,
    c1:f64,
    tc:f64,
    //Associative
    epsilon:f64,
    beta:f64,
    na:usize,
    nb:usize,
    nc:usize,
    typ:AssociationType

}
impl PureRecord for CPAPureRecord {
    
}
pub struct CPAParameters<C:CubicModel,R:RDF> ( pub CubicParameters<C>, pub ASCParameters<R>);


impl<C: CubicModel, R: RDF> Parameters<CPAPureRecord> for CPAParameters<C,R> {

    fn from_records(records:Vec<CPAPureRecord>)->Self {
        
        let mut cubic_records:Vec<CubicPureRecord>=Vec::new();
        let mut asc_records:Vec<AssociationPureRecord>=Vec::new();

        for record in records.iter(){
            cubic_records.push(CubicPureRecord::new(
                record.a0, 
                record.b, 
                record.c1, 
                record.tc));
            
            asc_records.push(AssociationPureRecord::new(
                record.epsilon,
                record.beta,
                record.na,
                record.nb,
                record.nc,
                record.typ.clone(),
                record.b)
            );
        }

        let cubic_params=CubicParameters::<C>::new(C::model(), cubic_records);
        let asc_params=ASCParameters::<R>::new(R::model(), asc_records);

        CPAParameters(cubic_params, asc_params)

    }
}
type EOS = E<CPA<SRK,ElliotRDF>>;
//ready-to parameters


pub fn water_acetic_acid()->EOS{
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

        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1,c2]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a1,a2]));
        let mut cpa = CPA::from_models(c, a);

        // let mut cpa=CPA::<ElliotRDF>::srk_from_records(
        //     vec![c1,c2],
        //     vec![a1,a2]);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij(0, 1, -0.222);
        cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
        //Create new State
        E::from_residual(cpa)

} 
pub fn water_octane_acetic_acid()->EOS{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c3=CubicPureRecord::new(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

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

        let a3=AssociationPureRecord::inert(0.1424e-3);


        //CPA eos
        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1,c2,c3]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a1,a2,a3]));
        let mut cpa = CPA::from_models(c, a);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij(0, 1, -0.222);
        cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
        //Create new State
        E::from_residual(cpa)

} 
pub fn acetic_acid_water()->EOS{
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


        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c2,c1]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a2,a1]));
        let mut cpa = CPA::from_models(c, a);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij(0, 1, -0.222);
        cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
        //Create new State
        E::from_residual(cpa)

} 

pub fn water_co2()->EOS{
            //Records
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
            0.0272e-3);


        //CPA eos
        
        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1,c2]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a1,a2]));
        let mut cpa = CPA::from_models(c, a);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij_temperature_dependent(0, 1,-0.15508,0.000877);

        //Set binary parameters 
        // cpa.cubic.set_binary(0,1, Some(0.000877),-0.15508 );
        cpa.assoc.set_binary(0, 1, AssociationRule::MCR1, None, Some(0.1836));

        //Create new State
        E::from_residual(cpa)

} 
pub fn co2_water()->EOS{
            //Records
        //1:Water, 2:Acetic Acid
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
            0.0272e-3);


        //CPA eos
        
        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c2,c1]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a2,a1]));
        let mut cpa = CPA::from_models(c, a);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij_temperature_dependent(0, 1,-0.15508,0.000877);

        //Set binary parameters 
        // cpa.cubic.set_binary(0,1, Some(0.000877),-0.15508 );
        cpa.assoc.set_binary(0, 1, AssociationRule::MCR1, None, Some(0.1836));

        //Create new State
        E::from_residual(cpa)

} 


pub fn methanol_2b()->EOS{
            //Records
        //1:metoh, 2:oct
        let c1=CubicPureRecord::new(0.40531, 0.0000309, 0.4310, 513.);

        let a1=AssociationPureRecord::associative(
            24591.0, 
            0.01610, 
            [1,1,0],
            0.0000309
        );

        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a1]));
        let cpa = CPA::from_models(c, a);

        
        //Create new State
        E::from_residual(cpa)
} 
pub fn methanol_3b()->EOS{
            //Records
        //1:metoh, 2:oct
        let c1=
        CubicPureRecord::
        new(
            4.5897e-1, 
            0.0334e-3, 
            1.0068,
            513.);

        let a1=AssociationPureRecord::associative(
            160.70e2, 
            34.4e-3, 
            [2,1,0],
            0.0334e-3
        );

        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a1]));
        let cpa = CPA::from_models(c, a);

        
        //Create new State
        E::from_residual(cpa)
} 


pub fn acoh_octane()->EOS{
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


        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1,c2]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a1,a2]));
        let mut cpa = CPA::from_models(c, a);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij(0, 1,  0.064 );

        //Set binary parameters 
        // cpa.cubic.set_binary(0,1, Some(0.0), 0.064 );
        
        //Create new State
        E::from_residual(cpa)

        
} 
pub fn octane_acoh()->EOS{
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
        let c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c2,c1]));
        let a = Associative::from_parameters(ASCParameters::<ElliotRDF>::from_records(vec![a2,a1]));
        let mut cpa = CPA::from_models(c, a);

        //Set binary parameters 
        cpa.cubic.parameters.set_kij(0, 1,  0.064 );
        //Set binary parameters 
        // cpa.cubic.set_binary(0,1, Some(0.0), 0.064 );
        
        //Create new State
        E::from_residual(cpa)

        
} 

#[cfg(test)]
pub mod tests{
    // use std::sync::Arc;

    // use approx::assert_relative_eq;
    // use ndarray::{Array1, Array2, arr1, array, linalg::Dot};

    // use crate::{models::cpa::{CPA, parameters::{acetic_acid_water, acoh_octane, co2_water, methanol_2b, methanol_3b, octane_acoh, water_acetic_acid, water_co2, water_octane_acetic_acid}}, state::{S, density_solver::DensityInitialization, eos::EquationOfState}};
    

    // #[test]
    // fn show_sites(){
    //     let eos=water_co2();
    //     let sites=&eos.residual.assoc.parameters.f;

    //     println!("Sites = {}" ,sites);
    //     let eos=water_octane_acetic_acid();
    //     let sites=&eos.residual.assoc.parameters.f;

    //     println!("Sites = {}" ,sites);
    // }
}