use pyo3::prelude::*;

use crate::{py_eos::PyEquationOfState, py_eq::PyPhaseEquilibrium, py_parameters::{pyrecords::{PyAssocRecord, PyCubicRecord}, PyCpaParameters, PyCubicParameters}, py_state::PyState};
mod py_eos;
mod py_parameters;
mod py_state;
mod py_eq;

#[pymodule]
fn reos(m: &Bound<'_, PyModule>) -> PyResult<()> {
    
    m.add_class::<PyState>()?;
    m.add_class::<PyEquationOfState>()?;


    m.add_class::<PyCpaParameters>()?;
    m.add_class::<PyCubicParameters>()?;
    m.add_class::<PyCubicRecord>()?;
    m.add_class::<PyAssocRecord>()?;


    m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
