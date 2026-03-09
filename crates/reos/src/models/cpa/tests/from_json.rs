use crate::models::{associative::{combining_rule::CombiningRuleOption, parameters::AssociationBinaryRecord}, cpa::parameters::CPABinaryRecord, cubic::parameters::{CubicBinaryRecord, Kij}};


#[test]
fn water_co2_binary(){
    

    let s = r#"
    {   

        "kij": {
            "a":-0.15508, 
            "b":0.000877
        },
        "assoc_rule": {
            "mcr1": {"kappa":0.9}
        }

        
    }
    "#;

    let val:CPABinaryRecord = serde_json::from_str(s).unwrap();
    let c = CubicBinaryRecord{ kij:Kij{a:-0.15508, b:0.000877}};
    let a = AssociationBinaryRecord::AssocRule(CombiningRuleOption::MCR1 { kappa: 0.9 });
    let tgt = CPABinaryRecord::new(Some(c), Some(a));
    assert_eq!(val.to_string(), tgt.to_string());
}