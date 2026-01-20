// ! Consts module exposing physical constants to Python

use pyo3::{pyclass, pymethods};
pub use reos::models::IDEAL_GAS_CONST as R;

#[pyclass]
pub struct Consts;

#[pymethods]
impl Consts {
    /// Ideal gas constant [J/(mol K)]
    #[staticmethod]
    pub fn ideal_gas_const()->f64{
        R
    }
}
