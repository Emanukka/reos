use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use reos::models::cpa::parameters::ASCParameters;
use reos::models::cubic::parameters::CubicParameters;
// use reos::models::{cpa::}
use crate::py_parameters::pyrecords::{PyAssocRecord, PyCubicRecord};
use reos::models::cpa::parameters::AssociationRule;
pub mod pyrecords;

#[derive(Clone)]
#[pyclass(name = "CPAParameters")]
pub struct PyCpaParameters{
    pub asc:ASCParameters,
    pub cub:CubicParameters,
}

#[pymethods]
impl PyCpaParameters {



    // pub fn from_json(json:String,components:Vec<String>)->Self{

    //     let c:Vec<&str>=components.iter().map(|a|a.as_str()).collect();
        
    //     PyCpaParameters{
    //         asc:ASCParameters::from_json(&json, c.clone()),
    //         cub:CubicParameters::from_json(&json, c)
    //     }


    // }

/// Parameters
/// ----------
///     cubic: List[CubicRecord],
///     assoc: List[AssociationRecord]
/// 
/// Returns
/// -------
    /// Return SRK-CPA Parameters from cubic & assoc records.
    #[staticmethod]
    #[pyo3(
        signature = (cubic,assoc),
        text_signature = "(cubic, assoc)"
    )]
    pub fn from_records(cubic:Vec<PyCubicRecord>,assoc:Vec<PyAssocRecord>)->Self{

        Self 
        {asc: ASCParameters::from_records(assoc.iter().map(|u|u.0.clone()).collect()), 
         cub: CubicParameters::from_records(cubic.iter().map(|u|u.0.clone()).collect())  }

    }
    #[pyo3(
    signature = (i,j,kij_a=0.0,kij_b=0.0),
    text_signature = "(i,j,kij_a=0.0,kij_b=0.0)",
    )]

    pub fn set_cubic_binary(&mut self,i:usize,j:usize,kij_a:f64,kij_b:f64){
        self.cub.set_binary(i, j, Some(kij_a), kij_b);
    }
    
    #[pyo3(
    signature = (i,j,rule,eps=None,beta=None),
    text_signature = "(i,j,rule,eps=None,beta=None)",
    )]
    pub fn set_assoc_binary(&mut self,i:usize,j:usize,rule:String,eps:Option<f64>,beta:Option<f64>)->PyResult<()>{

        match AssociationRule::try_from(rule){
            Ok(rule)=>{
                self.asc.set_binary_interaction(i, j, rule, eps, beta);
                Ok(())
            }
        Err(e)=> {
            Err(PyErr::new::<PyValueError, _>(
                e.to_string(),
            ))}
        }
    }

    pub fn as_string(&self)->String{

        let asc =format!("{}",&self.asc);
        let cub =format!("{}",&self.cub);
        cub + &asc

    }
}



#[derive(Clone)]
#[pyclass(name = "CubicParameters")]
pub struct PyCubicParameters(
    pub CubicParameters);

#[pymethods]
impl PyCubicParameters {

/// Parameters
/// ----------
/// 
/// List[CubicRecord],
/// 
/// Returns
/// -------
    #[staticmethod]
    #[pyo3(
        signature = (parameters),
        text_signature = "(parameters)"
    )]
    pub fn from_records(parameters:Vec<PyCubicRecord>)->Self{

        let p =CubicParameters::from_records(parameters.iter().map(|u|u.0.clone()).collect());
        
        Self(
            p
        )

    }

    pub fn as_string(&self)->String{

        let cub =format!("{}",&self.0);
        cub 

    }

    #[pyo3(
    signature = (i,j,kij_a=0.0,kij_b=0.0),
    text_signature = "(i,j,kij_a=0.0,kij_b=0.0)",
    )]

    pub fn set_cubic_binary(&mut self,i:usize,j:usize,kij_a:f64,kij_b:f64){
        self.0.set_binary(i, j, Some(kij_a), kij_b);
    }
}