


// Enum necessário p\ interface em python
// (pyo3 não suporta structs com atributos genéricos,
// - python é run-time comp.)
// 

use reos::{Array1, models::{cpa::{CPA, rdf::{CarnahanStarlingRDF, ElliotRDF}}, cubic::Cubic}, residual::Residual, state::eos::EosResult};

#[derive(Clone)]
pub enum ResidualModel{

    sCPA(CPA<ElliotRDF>),
    CPA(CPA<CarnahanStarlingRDF>),
    Cubic(Cubic),

} 



impl Residual for ResidualModel {
    
    fn components(&self)->usize {
        todo!()
    }
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        match self{
            ResidualModel::sCPA(eos)=>{
                eos.pressure(t, rho, x)
            },

            ResidualModel::CPA(eos)=>{
                eos.pressure(t, rho, x)
            },
            ResidualModel::Cubic(eos)=>{
                eos.pressure(t, rho, x)
            }


        }
    }

    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        match self{
            ResidualModel::sCPA(eos)=>{
                eos.residual_chemical_potential(t, rho, x)
            },
            ResidualModel::CPA(eos)=>{
                eos.residual_chemical_potential(t, rho, x)
            },
            ResidualModel::Cubic(eos)=>{
                eos.residual_chemical_potential(t, rho, x)
            },

        }
    }

    fn residual_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        match self{
            ResidualModel::sCPA(eos)=>{
                eos.residual_helmholtz(t, rho, x)
            },
            ResidualModel::CPA(eos)=>{
                eos.residual_helmholtz(t, rho, x)
            },
            ResidualModel::Cubic(eos)=>{
                eos.residual_helmholtz(t, rho, x)
            },

        }        
    }

    fn bmix(&self,x:&Array1<f64>)->f64 {
        match self{
            ResidualModel::sCPA(eos)=>{
                eos.bmix(x)
            },
            ResidualModel::CPA(eos)=>{
                eos.bmix(x)
            },
            ResidualModel::Cubic(eos)=>{
                eos.bmix(x)
            },

        }
    }
}
