use numpy::{IntoPyArray, PyArray1};
use pyo3::exceptions::{PyTypeError};
use pyo3::prelude::*;
use reos::models::cpa::associative::Associative;

use reos::models::cpa::CPA;
use reos::models::cubic::{Cubic, CubicModel};
use reos::state::eos::EquationOfState;
use std::sync::Arc;
use crate::py_eos::cpa::py_association::PyAssociation;
use crate::py_eos::py_residual::ResidualModel;
use crate::py_eos::PyEquationOfState;
use crate::py_parameters::PyCpaParameters;
use crate::py_state::PyState;

pub mod py_association;

#[pymethods]
impl PyEquationOfState {


    /// A CPA Equation Of State.
    /// 
    /// Parameters
    /// ----------
    /// 
    /// parameters : CPAParameters object
    /// 
    /// Returns
    /// -------
    /// CPA Equation Of State
    #[staticmethod]
    pub fn cpa(parameters: Bound<'_,PyCpaParameters>)->PyResult<Self>{

        let p: PyCpaParameters=parameters.extract()?;

        let cpa=CPA::new(
            Cubic::new(p.cub.clone(), CubicModel::SRK),
            Associative::new(p.asc.clone())
        );

        let r=ResidualModel::CPA(cpa);

        let arc=Arc::new(EquationOfState::from_residual(r));

        Ok(
            PyEquationOfState(arc)
        )
    }

    /// Returns the Association contribution object of EOS.
    /// 
    /// The object returned can calculate associative properties
    /// from the properties converged in the 'State initialization'
    /// e.g: calculate the fraction of non-bonded sites X
    /// 
    /// Warning! If the EoS doens't have an associative contribution,
    /// an error will be thrown
    fn get_association(&self)->PyResult<PyAssociation>{

        let residual=&self.0.residual;

        if let ResidualModel::CPA(cpa) = residual{
            Ok(
                PyAssociation(cpa.assoc.clone())
            )
        }
        else {

            Err(PyErr::new::<PyTypeError, _>("Error! EOS doens't contai Associative Contribution."))
        }
    }
}
#[pymethods]

impl PyState {
    
    /// Get the fraction of non-bonded sites at the current state.
    /// 
    /// Each element at position k represents the unbonded fraction of site k.
    /// Each site has a owner and a type.
    /// To know that, get the association contribution from the EquationOfState
    /// used to create the current State and use the function 'get_fmap'
    /// 
    /// This function will return an array of size equal to the total number of
    /// distinct sites.
    ///  
    /// Warning! If the EoS doens't have an associative contribution,
    /// an error will be thrown.
    /// 
    /// Returns
    /// -------
    /// Fraction of non-bonded sites X of site

    fn non_bonded_sites<'py> (&self,py:Python<'py>)->PyResult<Bound<'py, PyArray1<f64>>>{

        let residual=&self.0.eos.residual;


        if let ResidualModel::CPA(cpa) = residual{
            let t=self.0.t;
            let rho=self.0.t;
            let x=&self.0.x;

            Ok
            (
                cpa.assoc.X_tan(t, rho, x).unwrap().into_pyarray(py)
                
            )
        }else {
            Err(PyErr::new::<PyTypeError, _>("Error! State doens't contain Associative Contribution."))
        }

    }
}