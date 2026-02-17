use enum_dispatch::enum_dispatch;

use crate::models::cubic::options::OptionsParseError;


#[enum_dispatch]
pub trait CombiningRuleModel {

    fn apply(&self, ai_aj:[f64;2], bi_bj:[f64;2], kij:f64) -> [f64;2];

    fn dt(&self,  ai_aj:[f64;2], dai_daj:[f64;2], kij:f64, dkij:f64) -> f64;

}
#[derive(Clone, Debug)]
pub struct Classic;

impl CombiningRuleModel for Classic {

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
#[derive(Clone, Debug)]
pub enum CombiningRule {
    Classic,
}

impl Default for CombiningRule {
    fn default() -> Self {
        Self::Classic(Classic)
    }
}

impl std::str::FromStr for CombiningRule {
    type Err = OptionsParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {

        match s.to_lowercase().as_str() {
            "classic" => Ok(Self::Classic(Classic)),
            "" => Ok(Self::default()),
            _ => Err(OptionsParseError(format!("{} is not a combining rule implemented",s)))
        }
    }
}