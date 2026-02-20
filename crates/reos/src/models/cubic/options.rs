use serde::{Deserialize, Serialize};


pub use super::models::CubicModelOption;
pub use super::alpha::AlphaOption;
pub use super::combining_rule::CombiningRuleOption;
pub use super::mixing_rule::MixingRuleOption;

// #[derive(Clone, Debug)]
// pub struct CubicParamet{
//     pub model:CubicModels,
//     pub alpha:Alpha,
//     pub combr:CombiningRule,
//     pub mix:MixingRule,
// }

#[derive(Clone, Debug, Serialize, Deserialize,PartialEq)]

pub struct CubicOptions{
    pub cubic_model:CubicModelOption,
    #[serde(default)]
    pub alpha_model:AlphaOption,
    #[serde(default)]
    pub combining_rule:CombiningRuleOption,
    #[serde(default)]
    pub mixing_rule:MixingRuleOption,
    
}

// impl From<CubicOptions> for CubicOptions {

//     fn from(value: CubicOptions) -> Self {
        
//         CubicOptions { model: value.model.into(), alpha: value.alpha.into(), combr: value.combr.into(), mix: value.mix.into() }

//     }
    
// }
impl CubicOptions {

    pub fn new(cubic_model: CubicModelOption, alpha_model: AlphaOption,  combining_rule: CombiningRuleOption, mixing_rule:MixingRuleOption,) -> Self {

        Self { cubic_model, alpha_model, combining_rule, mixing_rule}
    }

    pub fn classic(cubic_model: CubicModelOption, alpha_model: AlphaOption) -> Self{

        Self { cubic_model, alpha_model, combining_rule: CombiningRuleOption::default(), mixing_rule: MixingRuleOption::default() }

    }

    pub fn classic_soave(model:CubicModelOption) -> Self {

        Self::classic(model, AlphaOption::Soave.into())
        // Self { model, alpha: Alpha::soave(), combr: CombiningRule::default(), mix: MixingRule::default() }
    }

}

// pub struct CubicOptionsHelper{
//     pub model:String,
//     pub alpha:String,
//     pub combr:String,
//     pub mix:String

// }


// impl TryFrom<CubicOptionsHelper> for CubicOptions {
//     type Error = OptionsParseError;
//     fn try_from(value: CubicOptionsHelper) -> Result<Self, Self::Error> {
        
//         let model = CubicModels::from_str(&value.model)?;
//         let combr = CombiningRule::from_str(&value.combr)?;
//         let alph
//     }
// }