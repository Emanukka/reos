use std::sync::Arc;

use numpy::{PyArray1, PyArrayMethods};
use pyo3::{Bound, PyResult, pyclass, pymethods};
use reos::{phase_equilibrium::tpd::MinTPD, state::density_solver::DensityInitialization};

use crate::{contribution::PyContribution, state::PyState};



use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;

// use crate::py_parameters::PyParameters;


#[pymethods]
impl PyState {


    #[doc = include_str!("../../docs/tpd/min_tpd_init.md")]
    #[pyo3(
    text_signature = "(xphase, xguess, tol=1e-8, it_max=100)"
    )]
    #[pyo3(signature = (xphase,xguess,tol=1e-8,it_max=100))]
    pub fn min_tpd<'py>(&self,xphase:&str,xguess:&Bound<'py, PyArray1<f64>>,tol:Option<f64>,it_max:Option<i32>)->PyResult<PyMinTPD>{

        let state = &self.0;
        let xphase = DensityInitialization::from_str(xphase);
        let xguess=xguess.to_owned_array();

        if xphase.is_err(){
            return Err(PyErr::new::<PyValueError, _>(
                                "`density_initialization` must be 'vapor' or 'liquid'.".to_string(),
                            ))
        }

        match state.min_tpd(xphase.unwrap(), xguess, tol, it_max) {
            Ok(tpd) => Ok(PyMinTPD(tpd)),
            Err(err) => Err(PyErr::new::<PyRuntimeError, _>(err.to_string()))
        }


    }

}
#[pyclass(name = "MinTPD")]
pub struct PyMinTPD(pub MinTPD<PyContribution>);


#[pymethods]
impl PyMinTPD {

    #[getter]
    pub fn dg(&self) -> f64 {
        self.0.dg
    }

    #[getter]
    pub fn state(&self) -> PyState {

        PyState(Arc::clone(&self.0.state))

    }

    #[getter]
    pub fn phase(&self) -> String {
        self.0.phase.to_string()
    }
}