use std::sync::Arc;

use pyo3::{PyResult, pyclass, pymethods};

// use pyo3::exceptions::{PyRuntimeError};
// use pyo3::prelude::*;
// use numpy::{IntoPyArray, PyArray1, PyArrayMethods };
// use reos::{state::{S, State, density_solver::DensityInitialization}};

use reos::scaling::dehlouz::{Scaling, ScalingRecord};

use crate::{contribution::PyContribution, state::PyState};
// 
/// A Thermodynamic State 
#[pyclass(name="Dehlouz")]
pub struct PyDehlouz(pub Scaling<PyContribution>);

#[pymethods]
impl PyDehlouz {

    #[new]
    pub fn new(state: PyState, parameters: ScalingRecord) -> Self{
        
        let inner = Scaling::new(&state.0, parameters);

        Self(inner)

    }

    /// Computes viscosity [Pa.s]
    pub fn viscosity(&self,sc:f64,) -> f64 {

        self.0.viscosity(sc)

    }

    /// Computes Y variable of entropy scaling: ln(ηreal/ηref)
    /// 
    /// Args: 
    /// - sc: float,
    ///       Tv-Residual entropy at (Tc, Pc) of EOS (?)
    pub fn yscaling(&self,sc:f64,) -> f64 {

        self.0.yscaling(sc)

    }
    /// Computes X variable of entropy scaling
    /// 
    /// Args: 
    /// - sc: float,
    ///       Tv-Residual entropy at (Tc, Pc) of EOS (?)
    pub fn xscaling(&self,sc:f64,) -> f64 {

        self.0.xscaling(sc)

    }

    pub fn __repr__(&self) -> PyResult<String> {
        Ok(
            self.0.to_string()
        )
    }

}
// Initializers
// #[pymethods]
// impl PyState {

    // #[doc = include_str!("../../docs/state/tpx.md")]
    // #[staticmethod]
    // #[pyo3(
    // text_signature = "(eos, t, p, x, phase = None)"
    // )]
    // #[pyo3(signature = (eos,temperature,pressure,x, phase = None))]
// }