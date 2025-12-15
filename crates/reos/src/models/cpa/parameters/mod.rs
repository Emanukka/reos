use serde::{Deserialize, Serialize};

use crate::models::associative::parameters::{ CombiningRule, AssociationPureRecord, AssociativeParameters};
use crate::models::cpa::{CPA, SCPAsrkCR1, SCPAsrkECR};
use crate::models::cpa::rdf::{CS, RdfModel, Kontogeorgis};

use crate::models::cubic::{Cubic, CubicModel, SRK};
use crate::models::cubic::parameters::{CubicParameters, CubicPureRecord};
use crate::parameters::{Parameters, PureRecord};
use crate::state::E;


#[derive(Serialize,Deserialize)]
pub struct CPAPureRecord{
    #[serde(flatten)]
    pub cub:CubicPureRecord,
    #[serde(flatten)]
    pub asc:AssociationPureRecord

}
impl PureRecord for CPAPureRecord {
    
}

impl CPAPureRecord{

    pub fn new(cub:CubicPureRecord,asc:AssociationPureRecord)->Self{
        Self { cub, asc }
    }

}

pub struct CPAParameters<C:CubicModel>{
    pub cubic: CubicParameters<C>, 
    pub assoc: AssociativeParameters,
}


impl<C: CubicModel> Parameters<CPAPureRecord> for CPAParameters<C> {

    fn from_records(records:Vec<CPAPureRecord>)->Self {
        
        let mut cubic_records:Vec<CubicPureRecord>=Vec::new();
        let mut asc_records:Vec<AssociationPureRecord>=Vec::new();

        for record in records.into_iter(){
            cubic_records.push(record.cub);
            asc_records.push(record.asc);
        }

        let cubic_params=CubicParameters::<C>::new(C::model(), cubic_records);
        let asc_params=AssociativeParameters::from_records(asc_records);

        CPAParameters{
            cubic:cubic_params,
            assoc:asc_params,
            }
    }
}


#[cfg(test)]
mod tests{
    use serde_json::from_str;

    use crate::{models::{associative::parameters::AssociativeParameters, cpa::parameters::CPAPureRecord}, parameters::Parameters};


    #[test]
    fn test1(){

        let data1 = r#"
        {
            "a0":   0.12277,
            "b":    0.0145e-3, 
            "kappa":0.6736, 
            "tc":   647.14,
            "na":   2,
            "nb":   2
        }
        "#;
        
        let data2 = r#"
        {
            "a0":   0.35079, 
            "b":    0.0272e-3, 
            "kappa":0.7602, 
            "tc":   304.12,
            "nb":   1
        }
        "#;

        let c1:CPAPureRecord = from_str(data1).unwrap();
        let c2:CPAPureRecord = from_str(data2).unwrap();

        // let c1_string = serde_json::to_string_pretty(&c1).unwrap();

        // println!("{}",c1_string)

        // let p = AssociativeParameters::from_records(vec![c1,c2]);

        // let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
        // let c2=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

        // let a1=AssociationPureRecord::associative(
        //     166.55e2, 
        //     0.0692, 
        //     [2,2,0],
        // );
        // let a2=AssociationPureRecord::solvate(
        //     [0,1,0]);

        // let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
        // let mut parameters = CPAParameters::from_records(records);
        // // parameters.cubic.set_kij(0, 1, -0.222);

        // parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
        // parameters.assoc.set_binary_from_owners(0, 1, None, Some(0.1836));

        // let cpa = SCPAsrkCR1::from_parameters(parameters);
        // //Create new State
        // E::from_residual(cpa)
    }
}