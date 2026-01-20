
use crate::{impl_pure_record,impl_binary_record,impl_parameters,impl_eos};
use crate::parameters::records::*;
use numpy::{PyArray1, PyArrayMethods, ToPyArray};
use reos::models::cpa::CPA;
use reos::models::cpa::rdf::Kontogeorgis;
use reos::models::cubic::models::CubicModels;

use reos::models::cpa::{ 
    parameters::{CPAParameters,CPABinaryRecord, CPAPureRecord,},
};
use pyo3::{PyErr, pymethods};

use crate::{contribution::PyContribution, eos::PyEquationOfState};
use reos::state::eos::EquationOfState;
use reos::parameters::{Parameters,PureRecord,BinaryRecord};

impl_pure_record!(PyCpaPureRecord,CPAPureRecord,"CPAPureRecord", "../../docs/cpa/pr.md");

impl_binary_record!(PyCPABinaryRecord,CPABinaryRecord,"CPABinaryRecord", "../../docs/cpa/br.md");

impl_parameters!(PyCPAParameters,CPAParameters,PyCpaPureRecord,PyCPABinaryRecord,"CPAParameters", CubicModels ,"../../docs/cpa/parameters.md",);

type SCPA = CPA<Kontogeorgis>;

impl_eos!(SCPA, SCPA, PyCPAParameters, scpa, "../../docs/cpa/eos_scpa.md");

#[pymethods]
impl PyEquationOfState {

    // #[doc = include_str!("../../docs/cpa/unbonded_sites_fraction.md")]
    #[pyo3(text_signature = "(temperature, density, x)")]
    pub fn unbonded_sites_fraction<'py>(&self,temperature:f64, density:f64, x: Bound<'py,PyArray1<f64>>)-> PyResult<Bound<'py, PyArray1<f64>>> {
            
        let py = x.py();    
        let contribution = self.0.residual();
        let x = x.to_owned_array();

        match &contribution{
            PyContribution::SCPA(inner)=>{
                
                Ok(inner.unbonded_sites(temperature, density, &x).to_pyarray(py))

            }
            PyContribution::CPA(inner)=>{
                
                Ok(inner.unbonded_sites(temperature, density, &x).to_pyarray(py))

            }

            _ =>  Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Only eos with Association contribution have this method",)
                )
        }

    } 

    #[pyo3(text_signature = "(temperature, density, x)")]
    pub fn association_constants<'py>(&self,temperature:f64, density:f64, x: Bound<'py,PyArray1<f64>>)-> PyResult<Bound<'py, numpy::PyArray2<f64>>> {
        
        let py = x.py();
        let contribution = self.0.residual();
        let x = x.to_owned_array();

        match &contribution{
            PyContribution::SCPA(inner)=>{
                
                Ok(inner.association_constants(temperature, density, &x).to_pyarray(py))

            }
            PyContribution::CPA(inner)=>{
                
                Ok(inner.association_constants(temperature, density, &x).to_pyarray(py))

            }

            _ =>  Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Only EoS with Association contribution have this method",)
                )
        }

    } 
}