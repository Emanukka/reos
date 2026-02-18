mod pure;
mod association;

pub mod recipes{

    use crate::{models::{associative::{parameters::{AssociationBinaryRecord, AssociationPureRecord}, sites::CombiningRule}, cpa::{Binary, CPA, Pure, SCPA, parameters::{CPABinaryRecord, CPAParameters, CPAPureRecord}}, cubic::{models::SRK, options::CubicOptions, parameters::{CubicBinaryRecord, CubicPureRecord}}}, parameters::{BinaryRecord, Parameters, PureRecord}};

    
    pub fn scpa(pure_records: Vec<Pure>, binary_records: Vec<Binary>)->Result<SCPA, Box<dyn std::error::Error>>{
        
        let parameters = CPAParameters::new(pure_records, binary_records, CubicOptions::classic_soave(SRK.into()))?;
        Ok(CPA::from_parameters(parameters))
    }


    pub fn water4c()->Pure{

        let c = CubicPureRecord::regressed_soave(0.12277, 0.0145e-3, 647.14, 0.6736, None);

        let a = AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0]);

        let m = CPAPureRecord::new(c, a);
            
        PureRecord::new(0.0, "water".to_string(), m)

    }

    pub fn acetic1a()->Pure{
        let c = CubicPureRecord::regressed_soave(0.91196, 0.0468e-3, 594.8, 0.4644, None);

        let a = AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1]);

        let m = CPAPureRecord::new(c, a);
        PureRecord::new(0.0, "acetic_acid".to_string(), m)

    }

    pub fn water4c_acetic1a()->Binary{

        // let c = CubicBinaryRecord::TemperatureIndependent { kij: -0.222 };
        let c = CubicBinaryRecord{ kij: -0.222, lij: 0.0 };
        
        let a = AssociationBinaryRecord {  epsilon: None, kappa: None, combining_rule: CombiningRule::ECR };
        
        let b = CPABinaryRecord::full(c, a);

        BinaryRecord::new(b, "water", "acetic_acid")

    }

    pub fn co2()->Pure{

        let c=CubicPureRecord::regressed_soave(0.35079, 0.0272e-3, 304.12, 0.7602, None);

        let a=AssociationPureRecord::solvate(
            [0,1,0]);

        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "co2".to_string(), m)

    }

    pub fn water4c_co2()->Binary{

        // let c = CubicBinaryRecord::TemperatureDependent { aij: -0.15508 , bij: 0.000877 };
        
        let c = CubicBinaryRecord{ kij: -0.15508 , lij: 0.000877 };
        
        let a = AssociationBinaryRecord {epsilon: None, kappa: Some(0.1836), combining_rule: CombiningRule::default() };
        
        let b = CPABinaryRecord::full(c, a);
        
        BinaryRecord::new(b, "water", "co2")

    }

    pub fn octane()->Pure{
        let c=CubicPureRecord::regressed_soave(34.8750e-1, 0.1424e-3, 568.7, 0.99415, None);
        let a=AssociationPureRecord::inert();


        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "octane".to_string(), m)    
    }
    pub fn acoh_octane()->Binary{


        // let c = CubicBinaryRecord::TemperatureIndependent { kij: 0.064 };
        let c = CubicBinaryRecord{ kij: 0.064, lij: 0.0 };

        let a = AssociationBinaryRecord {epsilon: None, kappa: None, combining_rule: CombiningRule::default() };
        let b = CPABinaryRecord::full(c, a);
        
        BinaryRecord::new(b, "acetic_acid", "octane")    

    } 
    pub fn methanol3b()->Pure{

            let c=
            CubicPureRecord::regressed_soave(
                4.5897e-1, 
                0.0334e-3, 
                1.0068,
                513.,
            None);

            let a=AssociationPureRecord::associative(
                160.70e2, 
                34.4e-3, 
                [2,1,0],
            );

        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "methanol".to_string(), m)    

    } 

}
