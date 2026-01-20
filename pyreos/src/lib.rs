use pyo3::prelude::*;


pub mod eos;
pub mod parameters;
pub mod state;
pub mod phase_equilibrium;
pub mod models;
pub mod contribution;
pub mod consts;

/// Macro to create and register a PyO3 submodule with multiple classes.
/// see: <https://github.com/PyO3/pyo3/issues/1517#issuecomment-3402169068>
/// 
/// Example:
/// ```rust
/// add_pymodule!(py, m, "my_module", [ClassA, ClassB])?;
/// ```
#[macro_export]
macro_rules! add_pymodule {
    ($py:expr, $parent:expr, $name:expr, [$($cls:ty),* $(,)?]) => {{
        let sub_module = pyo3::types::PyModule::new($py, $name)?;
        // add a lista de classes ao submodulo
        $(
            sub_module.add_class::<$cls>()?;
        )*

        //add o submodulo ao modulo mae
        $parent.add_submodule(&sub_module)?;
        $py.import("sys")?
            .getattr("modules")?
            .set_item(format!("reos.{}", $name), &sub_module)?;
        // add o submodulo ao sys.modules

        Ok::<_, pyo3::PyErr>(())
    }};
}


#[pymodule]
fn reos(m: &Bound<'_, PyModule>) -> PyResult<()> {
    
    // m.add_class::<PyState>()?;
    // m.add_class::<PyEquationOfState>()?;

    // m.add_submodule::<models::cubic::parameters::PyCubicPureRecord>()?;
    // m.add_class::<PyCpaParameters>()?;
    // m.add_class::<PyCubicParameters>()?;
    // m.add_class::<PyCubicPureRecord>()?;
    // m.add_class::<PyAssocRecord>()?;

    // Consts module
    add_pymodule!(m.py(), m, "consts", 
    [consts::Consts])?;

    //State module
    add_pymodule!(m.py(), m, "state", 
    [state::PyState])?;
    //Cubic module
    add_pymodule!(m.py(), m, "cubic", 
    [
    models::cubic::PyCubicPureRecord,
    models::cubic::PyCubicBinaryRecord,
    models::cubic::PyCubicParameters
    ])?;

    //CPA module
    add_pymodule!(m.py(), m, "cpa", 
    [models::cpa::PyCpaPureRecord,
    models::cpa::PyCPABinaryRecord,
    models::cpa::PyCPAParameters,
    
    ])?;

    //EoS module
    add_pymodule!(m.py(), m, "eos", 
    [eos::PyEquationOfState])?;


    // m.add_class::<PyPhaseEquilibrium>()?;
    Ok(())
}
