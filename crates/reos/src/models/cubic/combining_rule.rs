use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Serialize};



#[enum_dispatch]
pub trait CombiningRuleModel {

    fn apply(&self, ai_aj:[f64;2], bi_bj:[f64;2], kij:f64) -> [f64;2];

    fn dt(&self,  ai_aj:[f64;2], dai_daj:[f64;2], kij:f64, dkij:f64) -> f64;

    fn to_string(&self) -> String;
}   
#[derive(Clone, Debug, PartialEq)]
pub struct Classic;

impl CombiningRuleModel for Classic {

    fn to_string(&self) -> String {
        "Classic".into()
    }

    fn apply(&self, ai_aj:[f64;2], bi_bj:[f64;2], kij:f64) -> [f64;2] {
        
        let [ai, aj] = ai_aj;
        let [bi, bj] = bi_bj;

        let bij = 0.5 * (bi + bj);
        
        let sqrt = (ai * aj).sqrt();
        let aij = (1. - kij) * sqrt;

        [aij, bij]
    }

    fn dt(&self, a:[f64;2], da:[f64;2], kij:f64, dkij:f64) -> f64 {
        
        let [ai, aj] = a;
        let [dai, daj] = da;

        let sqrt = (ai * aj).sqrt();
        0.5 * (1. - kij) * ( dai * aj + daj * ai) / sqrt - sqrt * dkij

    }
}

#[enum_dispatch(CombiningRuleModel)]
#[derive(Clone, Debug, PartialEq)]
pub enum CombiningRule {
    Classic,
}



#[derive(Clone, Debug, Serialize, Deserialize,PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum CombiningRuleOption {
    Classic,
}

impl Default for CombiningRuleOption {
    fn default() -> Self {
        Self::Classic
    }
}

impl From<CombiningRuleOption> for CombiningRule {
    
    fn from(value: CombiningRuleOption) -> Self {
        
        match value {

            CombiningRuleOption::Classic => Classic.into(),

        }    
    }
}

#[cfg(test)]
mod tests {

    use super::*;


    #[test]
    fn assert_model_parse() {

        let classic:CombiningRule = CombiningRuleOption::Classic.into();

        assert_eq!(classic, CombiningRule::Classic(Classic));
        
    }
}