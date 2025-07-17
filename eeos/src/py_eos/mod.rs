use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use reeos::models::associative::Associative;
use reeos::models::cpa::CPA;
use reeos::models::cubic::{Cubic, CubicModel};
use reeos::parameters::{JsonStruct, Parameters};
use reeos::state::eos::EquationOfState;
use std::sync::Arc;
use reeos::residual::{Residual,ResidualModel};

use crate::py_parameters::PyCpaParameters;


#[derive(Clone)]

#[pyclass(name="EquationOfState")]
pub struct PyEquationOfState( pub Arc<EquationOfState<ResidualModel>> );

#[pymethods]
impl PyEquationOfState {

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


}



