use pyo3::types::PyAnyMethods;
use pyo3::{pymethods, Bound, PyResult};
use reos::models::cubic::{Cubic, CubicModel};
use reos::state::eos::EquationOfState;
use std::sync::Arc;
use crate::py_eos::py_residual::ResidualModel;
use crate::py_eos::PyEquationOfState;
// // use reos::residual::{Residual};
use crate::py_parameters::{ PyCubicParameters};


#[pymethods]
impl PyEquationOfState {


    /// A SRK Equation Of State.
    /// 
    /// Parameters
    /// ----------
    /// 
    /// parameters : CubicParameters object
    /// 
    /// Returns
    /// -------
    /// SRK Equation Of State
    #[staticmethod]
    pub fn srk(parameters: Bound<'_,PyCubicParameters>)->PyResult<Self>{

        let p: PyCubicParameters=parameters.extract()?;

        let eos=Cubic::new(p.0,CubicModel::SRK);

        let r=ResidualModel::SRK(eos);

        let arc=Arc::new(EquationOfState::from_residual(r));

        Ok(
            PyEquationOfState(arc)
        )
    }


}