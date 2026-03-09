use std::default;

// use crate::{impl_pure_record,impl_binary_record,impl_parameters,impl_eos}
use  crate::impl_eos;
use crate::impl_py_binary_record;
use crate::impl_py_parameters;
use crate::impl_py_pure_record;
use reos::models::cubic::Cubic;
use reos::models::cubic::{parameters::{CubicParameters,CubicBinaryRecord, CubicPureRecord}};

// use crate::{contribution::PyContribution, eos::PyEquationOfState};
use reos::state::eos::EquationOfState;
// use reos::parameters::{Parameters,PureRecord,BinaryRecord};






impl_py_pure_record!("CubicPureRecord", Cubic, "../../docs/cubic/pr.md");

impl_py_binary_record!("CubicBinaryRecord",Cubic, "../../docs/cubic/br.md");

impl_py_parameters!("CubicParameters", Cubic, "../../docs/cubic/parameters.md");

impl_eos!(Cubic, "../../docs/cubic/eos.md");


// use super::parameters::*;

// #[derive(Clone)]
// #[pyclass(name = "CubicPureRecord")]
// pub struct PyCubicPureRecord (PureRecord<CubicPureRecord>);

// #[derive(Clone)]
// #[pyclass(name = "CubicBinaryRecord")]
// pub struct PyCubicBinaryRecord (BinaryRecord<CubicBinaryRecord>);

// #[pymethods]
// impl PyCubicBinaryRecord {


//     #[staticmethod]
//     #[pyo3(signature = (**data))]
//     /// Cubic Binary Record Initializer
//     /// 
//     /// Parameters
//     /// ---------------
//     /// id1: str,
//     ///     name of component 1
//     /// id2: str,
//     ///     name of component 2
//     /// kij: dict[str, float], optional
//     ///     can have 2 fields: `a`, `b`, such that `kij = a + bT` 
//     /// 
//     pub fn new(data: Option<&Bound<'_, PyDict>>) -> PyResult<Self> {

//         match data {
//             Some(x) => Ok(Self(pythonize::depythonize::<BinaryRecord<CubicBinaryRecord>>(x)?)),
//             None => Err(pyo3::PyErr::new::<PyTypeError, _>("empty data"))
//         }
//     }

//     fn __repr__(&self) -> String {
//         self.0.to_string()
//     }

// }


// #[derive(Clone)]
// #[pyclass(name = "CubicParameters")]
// pub struct PyCubicParameters(CubicParameters);
// #[pymethods]
// impl PyCubicParameters {


//     #[staticmethod]
//     #[pyo3(text_signature = "(pure_records, binary_records = [], **opt)")]
//     #[pyo3(signature = (pure_records, binary_records = vec![], **opt))]
//     /// CubicParameters Initializer
//     /// 
//     /// Parameters
//     /// ---------------
//     /// args:
//     ///     pure_records: list[CubicPureRecord],
//     ///     binary_records: list[CubicBinaryRecord],
//     ///     opt: kwargs,
//     ///     available fields:
//     ///         - `cubic_model`: `srk`, `pr78`
//     ///         - `alpha_model`: `soave`, `twu91`, (default `soave`)
//     fn from_records<'py>(py:Python, pure_records: Vec<PyCubicPureRecord>, binary_records: Vec<PyCubicBinaryRecord>, opt:Option<&Bound<'_, PyDict>>) -> PyResult<Self> {

//         let pure_records = pure_records.into_iter()
//             .map(|r|r.0)
//             .collect();
        
//         let binary_records = binary_records.into_iter()
//             .map(|r| r.0)
//             .collect();

//         let default: Bound<'_, PyDict> = PyDict::new(py);
//         let opt: &Bound<'_, PyDict> = opt.unwrap_or(&default);

//         let opt = pythonize::depythonize(opt)?;
//         let res = CubicParameters::new(pure_records, binary_records, opt);
        
//         match res {
//             Ok(x) => Ok(Self(x)),
//             Err(e) => Err(pyo3::PyErr::new::<PyException, _>(e.to_string()))
//         }

//     }

//     fn __repr__(&self) -> String {
//         self.0.to_string()
//     }

// }
// impl Parameters


// impl_pure_record!(PyCubicPureRecord,CubicPureRecord,"CubicPureRecord", "../../docs/cubic/pr.md");

// impl_binary_record!(PyCubicBinaryRecord,CubicBinaryRecord,"CubicBinaryRecord", "../../docs/cubic/br.md");

// type PR78Parameters = CubicParameters<PR78_MARKER>;
// type PR76Parameters = CubicParameters<PR76>;

// impl_parameters!(PyCubicParameters, CubicParameters, PyCubicPureRecord, PyCubicBinaryRecord, "CubicParameters", CubicOptions, "../../docs/cubic/parameters.md",);
// 