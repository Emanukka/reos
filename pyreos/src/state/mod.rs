//! Thermodynamic State API
//! 
//! Provides Python bindings for the `State` struct from the `reos` crate.
//! Allows creating and manipulating thermodynamic states in Python.

use std::sync::Arc;

use pyo3::exceptions::{PyRuntimeError};
use pyo3::prelude::*;
use numpy::{IntoPyArray, PyArray1, PyArrayMethods,  };
use reos::{state::{S, State, density_solver::DensityInitialization}};

use crate::eos::{PyEquationOfState, PyProperties};
use crate::contribution::PyContribution;


/// A Thermodynamic State 
#[derive(Clone)]
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

    #[staticmethod]
    #[pyo3(signature = (eos,temperature,pressure, phase = None))]
    pub fn pure_tp<'py>(
        eos: Bound<'py,PyEquationOfState>,
        temperature:f64,
        pressure:f64,
        phase:Option<&str>
    )->PyResult<Self>{


        
        let s = phase.unwrap_or("stable");
        let phase = DensityInitialization::from_str(s);

        if let Err(e) = phase {

            return Err(PyErr::new::<PyRuntimeError, _>(e.to_string()))

        } else {
            
            let eos:PyEquationOfState = eos.extract()?;
            let res= State::new_tp(Arc::clone(&eos.0), temperature, pressure, Some(phase.unwrap()));            
            
            match res {
                Ok(state) => Ok(PyState(Arc::new(state))),
                Err(e) =>  Err(PyErr::new::<PyRuntimeError, _>(e.to_string()))
            }
        
        }
    }



}

#[pymethods]
impl PyState {



    /// Natural logarithm of the fugacity coefficient, defined as:
    /// 
    /// `
    /// ln(ϕᵢ) = ∂F/∂nᵢ - ln(Z)
    /// `
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

    /// Pressure `[Pa]`, defined as:
    /// 
    /// `
    /// P(T,ρ,x) = ZρRT
    /// `
    ///
    /// Returns
    /// -------
    /// float
    /// 
    #[getter]
    pub fn pressure(&self) -> f64 {
        self.0.p
    }
    /// Temperature.
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
    /// Compressibility factor, defined as:
    /// 
    /// `
    /// Z(T,ρ,x) = Zⁱᵍ + Zʳ = 1 - V∂F/∂V
    /// `
    ///
    /// Returns
    /// -------
    /// float
    pub fn compressibility(&self) -> f64 {

        self.0.eos.compressibility(self.0.t, self.0.d, &self.0.x)

    }
    
    /// Residual TP Gibbs energy in `[J / mol]`, defined as:
    /// 
    /// `
    /// Gʳ(T,P,x) = Aʳ + RTZʳ - RTln(Z)
    /// `
    pub fn gibbs(&self) -> f64 {

        self.0.gibbs()

    }

    /// Residual TP Entropy `[J / mol / K]`, defined as:
    /// 
    /// `
    /// Sʳ(T,P,x) =  Sʳ(T,ρ,x) + Rln(Z)
    /// `
    pub fn tp_entropy(&self) -> f64 {

        self.0.tp_entropy()

    }
    /// Residual TV Entropy `[J / mol / K]`, defined as:
    /// 
    /// `
    /// Sʳ(T,ρ,x) = R(- F - T∂F/∂T)
    /// `
    pub fn tv_entropy(&self) -> f64 {

        self.0.tv_entropy()

    }
    /// Residual TV Helmholtz free energy `[J / mol]`, defined as:
    /// 
    /// `
    /// Aʳ(T,ρ,x) = RT ⋅ F
    /// ` 
    pub fn helmholtz(&self) -> f64 {

        self.0.eos.helmholtz(self.0.t, self.0.d, &self.0.x)

    }



    pub fn __repr__(&self) -> PyResult<String> {
        Ok(self.0.to_string())
    }


}
