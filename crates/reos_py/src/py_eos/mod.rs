use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;
use reos::models::associative::Associative;
use reos::models::cpa::CPA;
use reos::models::cubic::{Cubic, CubicModel};
use reos::parameters::{JsonStruct, Parameters};
use reos::state::eos::EquationOfState;
use std::sync::Arc;
use crate::py_eos::py_association::PyAssociation;
// // use reos::residual::{Residual};
use crate::py_parameters::PyCpaParameters;


pub mod py_residual;
pub mod py_association;


use py_residual::ResidualModel;
#[derive(Clone)]

#[pyclass(name="EquationOfState")]
pub struct PyEquationOfState( pub Arc<EquationOfState<ResidualModel>> );

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




