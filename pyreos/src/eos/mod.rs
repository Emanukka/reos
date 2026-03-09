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

use numpy::{PyArray1, PyArrayMethods, ToPyArray};
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
    /// `Pⁱᵍ(T,ρ) = ρ R T`
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
    
    /// Pressure `[Pa]`, defined as:
    /// 
    /// `P(T,ρ,x) = ZρRT`
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

    /// Residual Helmholtz free energy `[J / mol]`, defined as:
    /// 
    /// `
    /// Aʳ(T,ρ,x) = RT ⋅ F
    /// ` 
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

        self.0.helmholtz(t, d, &x.to_owned_array())
    }

    /// Residual TV Entropy `[J / mol / K]`, defined as:
    /// 
    /// `
    /// Sʳ(T,ρ,x) = R(- F - T∂F/∂T)
    /// `
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
    pub fn tv_entropy<'py>(&self,t: f64, d: f64,x: &Bound<'py, PyArray1<f64>>)->f64 {

        self.0.tv_entropy(t, d, &x.to_owned_array())
    }

    /// Compressibility factor, defined as:
    /// 
    /// `
    /// Z(T,ρ,x) = Zⁱᵍ + Zʳ = 1 - V∂F/∂V
    /// `
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
    /// Natural logarithm of the fugacity coefficient, defined as:
    /// 
    /// `
    /// ln(ϕᵢ) = ∂F/∂nᵢ - ln(Z)
    /// `
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
    /// Residual TP Chemical potential in `[J / mol]`, defined as:
    /// 
    /// `
    /// μᵢʳ(T,P,x) = RTln(ϕᵢ)
    /// `
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

    /// Residual Gibbs energy in `[J / mol]`, defined as:
    /// 
    /// `
    /// Gʳ(T,P,x) = Aʳ + RTZʳ - RTln(Z)
    /// `
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
}

#[macro_export]
macro_rules! impl_eos {

    (
        $name:ident, $docdir: expr
    ) => {

        paste::paste!{
            #[pyo3::pymethods]
            impl crate::eos::PyEquationOfState {
                #[doc = include_str!($docdir)]
                #[staticmethod]
                #[pyo3(signature = (parameters))]
                pub fn [<$name: lower>](parameters: [<Py $name Parameters>]) -> Self {
                    
                    let p = parameters.0;
                    let model = $name::from(p);
                    let wrapper = crate::contribution::PyContribution::$name(model);
                    let e = reos::state::eos::EquationOfState::from_residual(wrapper);
                    crate::eos::PyEquationOfState(std::sync::Arc::new(e))

                }
            }
        }

   
    };
}




