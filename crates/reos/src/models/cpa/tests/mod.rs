mod pure;
mod mixture;
mod association;
pub mod recipes{

    use crate::{models::{associative::{parameters::{AssociationBinaryRecord, AssociationPureRecord}}, cpa::{Binary, CPA, Pure, parameters::{CPABinaryRecord, CPAOptions, CPAParameters, CPAPureRecord}, rdf::RDFmodelOption}, cubic::{models::{CubicModelOption, SRK}, options::CubicOptions, parameters::{CubicBinaryRecord, CubicPureRecord}}}, parameters::{BinaryRecord, Parameters, PureRecord}};

    use crate::models::associative::combining_rule::CombiningRuleOption;
    
    pub fn scpa(pure_records: Vec<Pure>, binary_records: Vec<Binary>)->Result<CPA, Box<dyn std::error::Error>>{

        let options = CPAOptions::default();
        let parameters = CPAParameters::new(pure_records, binary_records, options)?;
        println!("{}",&parameters);

        Ok(parameters.into())
    }
    pub fn cpa_cs(pure_records: Vec<Pure>, binary_records: Vec<Binary>)->Result<CPA, Box<dyn std::error::Error>>{

        let options = CPAOptions { rdf_model: RDFmodelOption::CS, cubic_options: CubicOptions::classic_soave(CubicModelOption::SRK) };
        let parameters = CPAParameters::new(pure_records, binary_records, options)?;
        println!("{}",&parameters);

        Ok(parameters.into())
    }


    pub fn water4c()->Pure{

        let c = CubicPureRecord::regressed_soave(0.12277, 0.014515e-3, 647.29, 0.67359, None);

        let a = AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0]);

        let m = CPAPureRecord::new(c, a);
            
        PureRecord::new(0.0, "water".to_string(), m)

    }

    pub fn acetic1a()->Pure{
        let c = CubicPureRecord::regressed_soave(0.91196, 0.0468e-3, 594.95, 0.4644, None);

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
        
        let a = AssociationBinaryRecord::CombiningRule(CombiningRuleOption::ECR);
        // let a = AssociationBinaryRecord {  epsilon: None, kappa: None, combining_rule: Some(CombiningRuleOption::ECR) };
        
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
        
        // let a = AssociationBinaryRecord {epsilon: None, kappa: Some(0.1836), combining_rule: None };
        let a = AssociationBinaryRecord::CombiningRule(CombiningRuleOption::MCR1 { kappa: 0.1836 });
        
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

        let b = CPABinaryRecord::only_c(c);
        
        BinaryRecord::new(b, "acetic_acid", "octane")    

    } 
    pub fn methanol3b()->Pure{

            let c=
            CubicPureRecord::regressed_soave(
                4.5897e-1, 
                0.0334e-3, 
                512.64,
                1.0068,
            None);

            let a=AssociationPureRecord::associative(
                16069.999999999998, 
                34.4e-3, 
                [2,1,0],
            );

        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "methanol".to_string(), m)    

    } 
    pub fn ethanol2b()->Pure{

            let c=
            CubicPureRecord::regressed_soave(
                0.85164, 
                0.0491e-3, 
                513.92,
                0.7502,
            None);

            let a=AssociationPureRecord::associative(
                 21500.0, 
                0.0083, 
                [1,1,0],
            );

        let m = CPAPureRecord::new(c, a);

        PureRecord::new(0.0, "methanol".to_string(), m)    

    } 

}
