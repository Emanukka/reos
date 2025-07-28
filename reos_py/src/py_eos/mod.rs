use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use reos::models::associative::Associative;
use reos::models::cpa::CPA;
use reos::models::cubic::{Cubic, CubicModel};
use reos::parameters::{JsonStruct, Parameters};
use reos::state::eos::EquationOfState;
use std::sync::Arc;
use reos::residual::{Residual,ResidualModel};
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



