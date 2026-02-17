use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use paste::paste;

pub mod twu91;
pub mod soave;
// pub mod twu91;

use soave::Soave;
use twu91::Twu91;

use crate::models::cubic::{models::CubicModels,};

pub trait AlphaModel {

    fn tag() -> Self;

    fn build<A: AsRef<str>>(names:&[A], records:Vec<AlphaRecord>, model: &CubicModels) -> Result<Alpha, AlphaError>;

    fn alpha(&self, t:f64, tc:&[f64]) -> Vec<f64>;

    fn dalpha_dt(&self, t:f64, tc:&[f64]) -> Vec<f64>;
    
    fn to_string(&self) -> String;

}   

#[derive(Serialize, Deserialize, Clone, Debug)]
#[serde(untagged)]
pub enum AlphaRecord {
    Soave{w:f64},
    SoaveRegressed{c1:f64},
    Twu91{l:f64, m:f64, n:f64} 
}

impl std::fmt::Display for AlphaRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        match self {
            
            Self::Soave { w } => write!(f, "w={w}"),
            Self::SoaveRegressed { c1 } => write!(f, "c1={c1}"),
            Self::Twu91 { l, m, n } => write!(f, "L={l}, M={m}, N={n}")

        }
    }

}
#[derive(Error, Debug, PartialEq)]
pub enum AlphaError {
    #[error("Components [{0}] alpha parameters didnt match the choosen Alpha Model")]
    AssociatedRecord(String)
}


macro_rules! impl_alpha {
    ($($var:ident),*) => {
        #[derive(Debug,Clone)]
        pub enum Alpha { $(
            $var($var),
        )*}

        impl Alpha {

            pub fn build<A: AsRef<str>>(&self, names:&[A], records:Vec<AlphaRecord>, model: &CubicModels) -> Result<Alpha, AlphaError>{
                match &self {
                    $(Alpha::$var(_) => $var::build(names, records, model),)*
                }
            }
            pub fn alpha(&self, t:f64, tc:&[f64]) -> Vec<f64> {
                match &self {
                    $(Alpha::$var(variant) => variant.alpha(t, tc),)*
                }
            }

            pub fn dalpha_dt(&self,t:f64, tc:&[f64]) -> Vec<f64> {
                match &self {
                    $(Alpha::$var(variant) => variant.dalpha_dt(t, tc),)*
                }
            }
            
            pub fn to_string(&self) -> String {
                match &self {
                    $(Alpha::$var(variant) => variant.to_string(),)*
                }
            }
            
            paste! {
                    $(
                        pub fn [<$var:lower>]() -> Alpha{
                            Alpha::$var($var::tag())
                        }
                    )*
            }
        }   
    }
}

impl_alpha!{Soave, Twu91}
