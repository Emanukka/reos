use std::f64::consts::SQRT_2;


use crate::models::IDEAL_GAS_CONST as R;


pub const SRK_KAPPA_FACTORS:  [f64; 3] = [0.480000, 1.57400, -0.17600];
pub const PR76_KAPPA_FACTORS: [f64; 3] = [0.374640, 1.54226, -0.26992];
pub const PR78_KAPPA_FACTORS: [f64; 4] = [0.374642, 1.48503, -0.164423, 0.016666];
use enum_dispatch::enum_dispatch;

#[enum_dispatch(CubicModels)]
pub trait CubicModel{

    // fn model() -> Self where Self: Sized;
    
    fn eps(&self)->f64;
    fn sig(&self)->f64;
    fn omega(&self,)->f64;
    fn psi(&self)->f64;
    fn kappa_from_w(&self,w:f64)->f64;

    fn acrit(&self,tc:f64,pc:f64)->f64{
        self.psi()*(R*tc).powf(2.0)/pc
    }
    fn bcrit(&self,tc:f64,pc:f64)->f64{
        self.omega()*R*tc/pc
    }
    fn to_string(&self)-> String;
}

pub struct SRK;
pub struct PR76;
pub struct PR78;

impl CubicModel for SRK{
    
    // fn model()->Self where Self: Sized {
    //     SRK
    // }
    fn omega(&self,)->f64 {
        0.08664
    }
    fn psi(&self)->f64 {
        0.42748
    }
    fn eps(&self)->f64 {
        0.0
    }
    fn sig(&self)->f64 {
        1.0
    }
    fn to_string(&self)-> String {
        "SRK".to_string()
    }

    fn kappa_from_w(&self,w:f64)->f64 {
        let factors = &SRK_KAPPA_FACTORS;
        factors[0] + w*factors[1] + w.powi(2)*factors[2]

    }

}

impl CubicModel for PR76 {
    
    fn omega(&self,)->f64 {
        0.07780
    }
    fn psi(&self)->f64 {
        0.45724
    }
    fn eps(&self)->f64 {
        1.0 - SQRT_2
    }

    fn sig(&self)->f64 {
        1.0 + SQRT_2
    }

    fn to_string(&self)-> String {
        "PR76".to_string()    
    }

    fn kappa_from_w(&self,w:f64)->f64 {
        let factors = &PR76_KAPPA_FACTORS;
        factors[0] + w*factors[1] + w.powi(2)*factors[2]

    }
}

impl CubicModel for PR78 {
    
    fn omega(&self,)->f64 {
        PR76.omega()
    }
    fn psi(&self)->f64 {
        PR76.psi()
    }
    fn eps(&self)->f64 {
        PR76.eps()
    }

    fn sig(&self)->f64 {
        PR76.sig()
    }

    fn to_string(&self)-> String {
        "PR78".to_string()    
    }

    fn kappa_from_w(&self,w:f64)->f64 {

        let factors76 = &PR76_KAPPA_FACTORS;
        let factors78 = &PR78_KAPPA_FACTORS;

        if w < 0.491 {
            factors76[0] + w*factors76[1] + w.powi(2) * factors76[2]

        } else{
            factors78[0] + w*factors78[1] + w.powi(2) * factors78[2] + w.powi(3) * factors78[3]
        }
        
    }
}

// #[derive(Serialize)]
#[enum_dispatch]
pub enum CubicModels{
    SRK,
    PR76,
    PR78,
}
// 

impl Default for CubicModels {

    fn default() -> Self {
        SRK.into()
    }

}

#[derive(Debug)]
pub struct CubicModelParseError;

impl std::fmt::Display for CubicModelParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "CubicModelParseError: invalid cubic model string")
    }
}
impl std::error::Error for CubicModelParseError {}



impl TryFrom<&str> for CubicModels {
    type Error = CubicModelParseError;

    
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value.to_lowercase().as_str(){
            "srk"=> Ok(SRK.into()),
            "pr76"=> Ok(PR76.into()),
            "pr78"=> Ok(PR78.into()),
            _=> Err(CubicModelParseError)
        }
    }
}




