use enum_dispatch::enum_dispatch;
use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::models::cubic::{alpha::soave::SoaveFlavor, parameters::CubicOptions};

pub mod models;
pub mod twu91;
pub mod soave;
// pub mod twu91;

use soave::Soave;
use twu91::Twu91;

#[enum_dispatch(Alphaa)]
pub trait AlphaModel {

    fn alpha(&self, t:f64, tc:&[f64]) -> Vec<f64>;

    fn dalpha_dt(&self, t:f64, tc:&[f64]) -> Vec<f64>;
    fn to_string(&self) -> String;

}   
#[enum_dispatch]
#[derive(Clone,Debug)]
pub enum Alphaa {
    Soave,
    Twu91

}

#[derive(Serialize, Deserialize, Clone, Debug)]
#[serde(rename_all = "lowercase")]
#[serde(rename = "alpha")]
#[serde(untagged)]
pub enum AlphaRecord {
    Soave(SoaveFlavor),
    Twu91{l:f64, m:f64, n:f64} 
}

#[derive(Error, Debug)]
pub enum AlphaError {
    #[error("Component '{0}' alpha parameters didnt match the choosen Alpha Model")]
    AssociatedRecord(String)
}

pub enum AlphaOption {
    Soave,
    Twu91
}

impl Alphaa {

    pub fn build<A: AsRef<str>>(names:&[A], records: Vec<AlphaRecord>,options: &CubicOptions) -> Result<Alphaa, Vec<AlphaError>>{

        match options.alpha {

            AlphaOption::Soave => Soave::from_records(names, records, options),
            AlphaOption::Twu91 => Twu91::from_records(names, records, options),

        }
    }
}