
pub use pyo3::exceptions::PyTypeError;
pub use pyo3::types::PyDict;
pub use pyo3::{Bound, PyResult, pyclass, pyclass_init, pymethods};

/// Macro to generate concrete Python wrapper-like structs for a `PureRecord<T>` generic, 
/// given the corresponding concrete pure record model  `T`.
#[macro_export]
macro_rules! impl_pure_record {
    ($name: ident, $type: ident, $pyname: expr, $docdir: expr) => {
        #[pyclass(name = $pyname)]
        pub struct $name (
            pub PureRecord<$type>
        );
        
        impl Into<PureRecord<$type>> for $name {
            fn into(self) -> PureRecord<$type>{
                self.0
            }
        }

        impl Clone for $name {
            fn clone(&self) -> Self {
                Self(self.0.clone())
            }
        }

   
        #[pymethods]
        impl $name {

            #[doc = include_str!("../../docs/serializer/to_json_vec.md")]
            #[staticmethod]
            #[pyo3(text_signature = "(file_name, records, build = false)")]
            #[pyo3(signature = (file_name, records, build = false))]
            pub fn to_json_vec(file_name: &str, records:Vec<Self>, build: bool) -> PyResult<String> {

                let records = records.into_iter().map(|r| r.0).collect::<Vec<PureRecord<$type>>>();

                match reos::parameters::writer::to_json_vec::<PureRecord<$type>>(file_name, records, build) {
                    Ok(json_str) => Ok(json_str),
                    Err(e) => Err(pyo3::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Serialization error: {}", e))),
                }
            }
            fn __repr__(&self) -> String {
                self.0.to_string()
            }

            #[doc = include_str!($docdir)]
            #[staticmethod]
            #[pyo3(signature = (name, molar_weight, **parameters))]
            pub fn new(name:&str, molar_weight:f64, parameters:Option<&Bound<'_, PyDict>>)->PyResult<Self>{
                
                if let Some(x) = parameters {

                    let model:$type = pythonize::depythonize(x)?;
                    
                    let pr = PureRecord::new(molar_weight, name.to_string(), model);

                    return  Ok(Self(pr))
                    
                } else {
                Err(pyo3::PyErr::new::<PyTypeError, _>("parameters empty!"))
                }

            }
        }
                
    };
}
/// Macro to generate concrete Python wrapper-like structs for a `BinaryRecord<T>` generic, 
/// given the corresponding concrete binary record model  `T` .
#[macro_export]
macro_rules! impl_binary_record {
    ($name: ident, $type: ident, $pyname: expr,$docdir: expr) => {

        #[pyclass(name = $pyname)]
        pub struct $name (
            pub BinaryRecord<$type>
        );
        
        impl Into<BinaryRecord<$type>> for $name {
            fn into(self) -> BinaryRecord<$type>{
                self.0
            }
        }
        impl Clone for $name {
            fn clone(&self) -> Self {
                Self(self.0.clone())
            }
        }

        #[pymethods]
        impl $name {
            
            #[doc = include_str!("../../docs/serializer/to_json_vec.md")]
            #[staticmethod]
            #[pyo3(text_signature = "(file_name, records, build = false)")]
            #[pyo3(signature = (file_name, records, build = false))]            
            pub fn to_json_vec(file_name: &str, records:Vec<Self>, build: bool) -> PyResult<String> {

                let records = records.into_iter().map(|r| r.0).collect::<Vec<BinaryRecord<$type>>>();
                match reos::parameters::writer::to_json_vec::<BinaryRecord<$type>>(file_name, records, build) {
                    Ok(json_str) => Ok(json_str),
                    Err(e) => Err(pyo3::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!("Serialization error: {}", e))),
                }
            }
            fn __repr__(&self) -> String {
                self.0.to_string()
            }

            #[doc = include_str!($docdir)]
            #[staticmethod]
            #[pyo3(signature = (id1, id2, **parameters))]
            pub fn new(id1:&str,id2:&str, parameters:Option<&Bound<'_, PyDict>>)->PyResult<Self>{
                
                if let Some(x) = parameters {

                    let model:$type = pythonize::depythonize(x)?;
                    
                    let pr = BinaryRecord::new(model, id1.into(), id2.into());

                    return  Ok(Self(pr))
                    
                } else {
                Err(pyo3::PyErr::new::<PyTypeError, _>("parameters empty!"))
                }

            }
        }
                
    };
}
