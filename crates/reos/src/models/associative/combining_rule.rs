
use serde::{Serialize, Deserialize};

#[derive(Serialize,Clone,Copy,PartialEq,Debug,Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CombiningRuleOption{

    CR1,
    ECR,
    MCR1{kappa:f64},

}


impl std::fmt::Display for CombiningRuleOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self{
            CombiningRuleOption::CR1 => write!(f,"CR1"),
            CombiningRuleOption::MCR1{kappa} => write!(f,"m-CR1(kappa={kappa}"),
            CombiningRuleOption::ECR => write!(f,"ECR"),
        }
    }
}


impl Default for CombiningRuleOption {
    fn default() -> Self {
        CombiningRuleOption::CR1
    }
}

impl Into<CombiningRuleOption> for Option<CombiningRuleOption> {

    fn into(self) -> CombiningRuleOption {

        self.unwrap_or_default()
    }
}