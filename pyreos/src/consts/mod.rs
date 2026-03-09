// ! Consts module exposing physical constants to Python

use pyo3::{pyclass, pymethods};
pub use reos::models::R_GAS;

#[pyclass]
pub struct Consts;

#[pymethods]
impl Consts {
    /// Ideal gas constant [J/(mol K)]
    #[staticmethod]
    pub fn ideal_gas_const()->f64{
        R_GAS
    }
}
