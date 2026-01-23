//! Thermodynamic State API
//! 
//! Provides Python bindings for the `State` struct from the `reos` crate.
//! Allows creating and manipulating thermodynamic states in Python.

use std::sync::Arc;

use pyo3::exceptions::{PyRuntimeError};
use pyo3::prelude::*;
use numpy::{IntoPyArray, PyArray1, PyArrayMethods };
use reos::{state::{S, State, density_solver::DensityInitialization}};

use crate::eos::{PyEquationOfState, PyProperties};
use crate::contribution::PyContribution;


/// A Thermodynamic State 
#[pyclass(name="State")]
pub struct PyState(pub Arc<S<PyContribution>>);

// Initializers
#[pymethods]
impl PyState {

    #[doc = include_str!("../../docs/state/tpx.md")]
    #[staticmethod]
    #[pyo3(
    text_signature = "(eos, t, p, x, phase = None)"
    )]
    #[pyo3(signature = (eos,temperature,pressure,x, phase = None))]
    pub fn tpx<'py>
    (
        eos: Bound<'py, PyEquationOfState>,
        temperature:f64,
        pressure:f64,
        x:&Bound<'py, PyArray1<f64>>,
        phase:Option<&str>,
    )->PyResult<Self>{
        
        let eos:PyEquationOfState = eos.extract()?; 
        let x = x.to_owned_array();
        let s = phase.unwrap_or("stable");
        let phase = DensityInitialization::from_str(s);

        if let Err(e) = phase {

            return Err(PyErr::new::<PyRuntimeError, _>(e.to_string()))

        } else {
            
            let res= State::new_tpx(Arc::clone(&eos.0), temperature, pressure, x, Some(phase.unwrap()));
            
            match res {
            Ok(state) => Ok(PyState(Arc::new(state))),
            Err(e) =>  Err(PyErr::new::<PyRuntimeError, _>(e.to_string()))
            }
        
        }
    }
    
    #[doc = include_str!("../../docs/state/tdx.md")]
    #[staticmethod]
    #[pyo3(
    text_signature = "(eos, t, d, x)"
    )]
    #[pyo3(signature = (eos,temperature,density,x))]

    pub fn tdx<'py>
    (
        eos: Bound<'py,PyEquationOfState>,
        temperature:f64,
        density:f64,
        x:&Bound<'py, PyArray1<f64>>,
        
    )->PyResult<Self>{

        let eos:PyEquationOfState = eos.extract()?;
        let x = x.to_owned_array();
        let s= State::new_trx(Arc::clone(&eos.0), temperature, density, x);

        Ok(PyState(Arc::new(s)))

    }


}

#[pymethods]
impl PyState {



    /// Computes the logarithmic fugacity coefficient.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    pub fn lnphi<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {

        let state=&self.0;
        state.lnphi().into_pyarray(py)

    }

    /// Computes the max-density of the mixture.
    ///
    /// Returns
    /// -------
    /// float
    pub fn max_density(&self) -> f64 {
        self.0.max_density()
        }

    /// Gets the Pressure.
    ///
    /// Returns
    /// -------
    /// float
    /// 
    #[getter]
    pub fn pressure(&self) -> f64 {
        self.0.p
    }
    /// Gets the Temperature.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    pub fn temperature(&self) -> f64 {
        self.0.t
    }
    /// Gets the volume.
    ///
    /// Returns
    /// -------
    /// float
    #[getter]
    pub fn volume(&self) -> f64 {
        1./self.0.d
    }
    /// Gets the mass density.
    /// 
    /// d = ρ * M, M = <MW, x>
    /// 
    /// Returns
    /// -------
    /// float
    #[getter]
    pub fn mass_density(&self) -> f64 {
        self.0.mass_density()
    }

    /// Gets the density.
    ///
    /// Returns
    /// -------
    /// float
    /// 
    #[getter]
    pub fn density(&self) -> f64 {
        self.0.d
    }

    /// Gets the composition z of the State.
    ///
    /// Returns
    /// -------
    /// np.ndarray[float]
    #[getter]
    pub fn composition<'py>(&self,py:Python<'py>) -> Bound<'py, PyArray1<f64>> {

        self.0.x.clone().into_pyarray(py)
    }

    pub fn get_properties(&self) -> PyResult<PyProperties> {

        let props = self.0.eos.get_properties();
        Ok(PyProperties(props.clone()))
    }

    /// Gets the molar weight of the components.
    ///
    /// Returns
    /// -------
    /// np.ndarray[float]
    #[getter]
    pub fn molar_weight<'py>(&self,py:Python<'py>) -> Bound<'py, PyArray1<f64>> {

        self.0.molar_weight().clone().into_pyarray(py)
    }
    /// Compressibility factor Z.
    ///
    /// Returns
    /// -------
    /// float
    pub fn compressibility(&self) -> f64 {

        self.0.eos.compressibility(self.0.t, self.0.d, &self.0.x)

    }
    
    /// Residual Gibbs free energy in J / mol
    pub fn gibbs(&self) -> f64 {

        self.0.gibbs()

    }

    /// Residual Entropy in J / mol / K
    pub fn entropy(&self) -> f64 {

        self.0.entropy()

    }
    /// Residual Helmholtz free energy in J / mol
    pub fn helmholtz(&self) -> f64 {

        self.0.eos.helmholtz_isov(self.0.t, self.0.d, &self.0.x)

    }



    pub fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "State(t = {:.3} K, p = {:.6} Pa, ρ = {:.6} mol/m³)",
            self.temperature(),
            self.pressure(),
            self.density()
        ))
    }


}
