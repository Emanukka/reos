mod pure;
mod mixture;

pub mod recipes {
    


    use crate::{models::associative::{Associative, parameters::{AssociationBinaryRecord, AssociationPureRecord, AssociativeParameters}, sites::CombiningRule}, parameters::{BinaryRecord, Parameters, records::PureRecord}};


    type Pure = PureRecord<AssociationPureRecord>;
    type Binary = BinaryRecord<AssociationBinaryRecord>;
    
    pub fn water() -> Pure {

        let m = AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0]);

            
        let pr = PureRecord::new(0.0, "water".to_string(), m);    
        
        pr
    }
    
    pub fn carbon_dioxide() -> Pure {
        PureRecord::new(0.,"carbon dioxide", AssociationPureRecord::solvate(
            [0,1,0]))

    }

    pub fn water_carbon_dioxide() -> Binary {

        BinaryRecord::new(
            AssociationBinaryRecord {epsilon: None, kappa: Some(0.1836), combining_rule: CombiningRule::default()},
            "water",
            "carbon dioxide"
        )

    }


    pub fn acetic_acid() -> Pure {


        let m = AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1]);

        PureRecord::new(0.0, "acetic_acid", m)
    }

    pub fn water_acetic_acid() -> Binary {
        
        let b = AssociationBinaryRecord {  epsilon: None, kappa: None, combining_rule: CombiningRule::ECR };
        

        BinaryRecord::new(b, "water", "acetic_acid")
    }
    
    pub fn methanol()-> Pure{

        let m=AssociationPureRecord::associative(
            160.70e2, 
            34.4e-3, 
            [2,1,0],);


        let pr = PureRecord::new(0.0, "methanol".to_string(), m);    
        // let p = AssociativeParameters::new(vec![pr], vec![], ()).unwrap();
        // let asc = Associative::from_parameters(p);

        // asc   
        pr
    } 
    

}