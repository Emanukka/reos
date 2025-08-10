use pyo3::prelude::*;
use reos::state::eos::EquationOfState;
use std::sync::Arc;

pub mod cpa;
pub mod py_residual;
pub mod cubic;

use py_residual::ResidualModel;
#[derive(Clone)]

#[pyclass(name="EquationOfState")]
pub struct PyEquationOfState( pub Arc<EquationOfState<ResidualModel>> );








