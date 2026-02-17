use std::str::FromStr;

use thiserror::Error;

use crate::models::cubic::{alpha::Alpha, combining_rule::CombiningRule, mixing_rule::MixingRule, models::CubicModels};

#[derive(Error, Debug)]
#[error("{}",self.0)]
pub struct OptionsParseError(pub String);

#[derive(Clone, Debug)]
pub struct CubicOptions{
    pub model:CubicModels,
    pub alpha:Alpha,
    pub combr:CombiningRule,
    pub mix:MixingRule,
    
}

impl CubicOptions {

    pub fn new(model: CubicModels, alpha: Alpha,  combr: CombiningRule, mix:MixingRule,) -> Self {

        Self { model, alpha, combr, mix}
    }

    pub fn classic(model: CubicModels, alpha: Alpha) -> Self{

        Self { model, alpha, combr: CombiningRule::default(), mix: MixingRule::default() }

    }
    pub fn classic_soave(model:CubicModels) -> Self {

        Self::classic(model, Alpha::soave())
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