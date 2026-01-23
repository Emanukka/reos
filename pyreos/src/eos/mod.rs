//! Equation of State API for Python.
//!
//! Defines the `PyEquationOfState` struct and its methods to be used individually
//! or to create a `PyState` object.
//! 
//! # The `impl_eos!` macro 
//! 
//! This macro is used to implement constructors for different the variants of `PyContribution`.
//! 
//! 
//! Example:
//! ```rust
//! impl_eos!(SRK, SRKContribution, PySrkParameters, "srk");
//! ```

use std::sync::Arc;

use numpy::{PyArray, PyArray1, PyArrayMethods, ToPyArray};
use pyo3::prelude::*;
use reos::{parameters::Properties, state::eos::EquationOfState};
// use std::sync::Arc;


use crate::contribution::PyContribution;

/// A Equation of State.
/// 
/// Can be used individually or to create a `State` object.
#[pyclass(name="EquationOfState")]
#[derive(Clone)]
pub struct PyEquationOfState( pub Arc<EquationOfState<PyContribution>> );

/// Properties of the Equation of State.
///
/// Components names and molar weights. 
#[pyclass(name="Properties")]
#[derive(Clone)]
pub struct PyProperties(pub Properties);

#[pymethods]
impl PyProperties {
    
        
    #[getter]
    pub fn molar_weight<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        // todo: modify to not copy the array data
        let molar_weight = &self.0.molar_weight;
        let mw = molar_weight.to_pyarray(py);
        mw
        
    }

    #[getter]
    pub fn components(&self) -> Vec<String> {
        self.0.names.clone()
    }

    pub fn __repr__(&self) -> PyResult<String> {
        Ok(format!(
            "Properties(components: {:?}, molar_weight: {:?})",
            self.0.names, self.0.molar_weight.as_slice()
        ))
    }
}
#[pymethods]
impl PyEquationOfState {
    
    /// Calculate the ideal gas pressure:
    /// 
    /// `P_ig = d R T`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    pub fn ideal_gas_pressure(&self,t: f64,d: f64)->f64{
        self.0.ideal_gas_pressure(t, d)
    }
    
    /// Calculate pressure:
    /// 
    /// `P = P_res + P_ig`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// float
    ///     Pressure [Pa]
    pub fn pressure<'py>(&self,t: f64,d: f64,x: &Bound<'py, PyArray1<f64>>)->f64{

        self.0.pressure(t, d, &x.to_owned_array())
    }

    /// Calculate residual molar helmholtz energy 
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// float
    ///     Residual molar Helmholtz energy [J/mol]
    pub fn helmholtz<'py>(&self,t: f64,d: f64,x: &Bound<'py, PyArray1<f64>>)->f64 {

        self.0.helmholtz_isov(t, d, &x.to_owned_array())
    }

    /// Calculate residual molar entropy
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// float
    ///     Residual molar Entropy [J/(mol K)]
    pub fn entropy<'py>(&self,t: f64,d: f64,x: &Bound<'py, PyArray1<f64>>)->f64 {

        self.0.entropy(t, d, &x.to_owned_array())
    }

    /// Calculate compressibility factor:
    /// `Z = P / P_ig`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    pub fn compressibility<'py>(&self,t:f64,d:f64,x:&Bound<'py, PyArray1<f64>>)->f64{
        self.0.compressibility(t, d, &x.to_owned_array())

    }
    /// Calculate the logarithm of fugacity coefficient of each component in the mixture:
    /// 
    /// `lnphi[i] = mu_res[i]/RT - ln_Z`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// numpy.ndarray[float]
    ///     Logarithm of fugacity coefficients
    pub fn lnphi<'py>(&self,t:f64,d:f64,x:&Bound<'py, PyArray1<f64>>)->Bound<'py, PyArray1<f64>>{

        self.0.lnphi(t, d, &x.to_owned_array()).to_pyarray(x.py())
        
    }

    /// Calculate the residual molar chemical potential of each component in the mixture:
    /// 
    /// `mu_res[i] = RT lnphi[i]`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// numpy.ndarray[float]
    ///     Residual molar chemical potential [J/(mol)]
    pub fn chem_pot<'py>(&self,t:f64,d:f64,x:&Bound<'py, PyArray1<f64>>)->Bound<'py, PyArray1<f64>>{

        self.0.chem_pot(t, d, &x.to_owned_array()).to_pyarray(x.py())
        
    }

    /// Calculate the residual molar gibbs energy of the mixture:
    /// 
    /// `Gres = Ares + RT ( Zres - lnZ )`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// numpy.ndarray[float]
    ///     Residual molar gibbs energy [J/(mol)]
    pub fn gibbs<'py>(&self,t:f64, d:f64, x:&Bound<'py, PyArray1<f64>>)->f64{

        self.0.gibbs(t, d, &x.to_owned_array())
        
    }

    // pub fn properties(&self) -> Vec<String> {
        // self.0.properties()
    // }
}

#[macro_export]
macro_rules! impl_eos {

    (
        $variant:ident , $type:ident, $param:ident, $pyname: ident, $docdir: expr
    ) => {

        #[pymethods]
        impl PyEquationOfState {
            #[doc = include_str!($docdir)]
            #[staticmethod]
            #[pyo3(signature = (parameters))]
            pub fn $pyname(parameters: $param) -> Self {
                
                let p = parameters.0;
                let model = $type::from_parameters(p);
                let wrapper = PyContribution::$variant(model);
                let e = EquationOfState::from_residual(wrapper);
                PyEquationOfState(std::sync::Arc::new(e))

            }
        }
   
    };
}




