


// Enum necessário p\ interface em python
// (pyo3 não suporta structs com atributos genéricos,
// - python é run-time comp.)
// 

use reos::{models::{associative::Associative, cpa::CPA}, residual::Residual, state::eos::EosResult, Array1};

#[derive(Clone)]
pub enum ResidualModel{

    CPA(CPA),

} 



impl Residual for ResidualModel {
    
    fn components(&self)->usize {
        todo!()
    }
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        match self{
            ResidualModel::CPA(eos)=>{
                eos.pressure(t, rho, x)
            },


        }
    }

    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        match self{
            ResidualModel::CPA(eos)=>{
                eos.residual_chemical_potential(t, rho, x)
            },

        }
    }

    fn bmix(&self,x:&Array1<f64>)->f64 {
        match self{
            ResidualModel::CPA(eos)=>{
                eos.bmix(x)
            },

        }
    }
}
