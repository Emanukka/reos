//! Parameters API for Equation Of State models.
//! 
//! Provides macros to generate Python classes 
//! for Parameters, PureRecord and BinaryRecord, since
//! this objects depend on generic model types from `reos` crate.
//! 
//! Current macros:
//! - `impl_py_parameters!`
//! - `impl_py_pure_record!`
//! - `impl_py_binary_record!`
//! 

#[macro_export]
macro_rules! impl_py_pure_record {
    ($pyname:literal, $name: ident, $docdir: expr) => {

        use reos::parameters::PureRecord;

        paste::paste!{
            #[derive(Clone)]
            #[pyo3::pyclass(name = $pyname)]
            pub struct [<Py $name PureRecord>](pub PureRecord<[<$name PureRecord>]>);

            #[pyo3::pymethods]
            impl [<Py $name PureRecord>] {

                    #[staticmethod]
                    #[pyo3(signature = (**data))]
                    #[doc = include_str!($docdir)]
                    pub fn new(data: Option<&pyo3::Bound::<'_,pyo3::types::PyDict>>) -> pyo3::PyResult<Self> {

                        match data {
                            Some(x) => Ok(Self(pythonize::depythonize::<PureRecord<[<$name PureRecord>]>>(x)?)),
                            None => Err(pyo3::PyErr::new::<pyo3::exceptions::PyTypeError, _>("empty data"))
                        }
                    }

                    fn __repr__(&self) -> String {
                        self.0.to_string()
                    }
            }
        }
    }
}
#[macro_export]
macro_rules! impl_py_binary_record {
    ($pyname:literal, $name: ident, $docdir: expr) => {

        use reos::parameters::BinaryRecord;

        paste::paste!{
            #[derive(Clone)]
            #[pyo3::pyclass(name = $pyname)]
            pub struct [<Py $name BinaryRecord>](pub BinaryRecord<[<$name BinaryRecord>]>);
            
            #[pyo3::pymethods]
            impl [<Py $name BinaryRecord>] {
                
                    #[staticmethod]
                    #[pyo3(signature = (**data))]
                    #[doc = include_str!($docdir)]

                    pub fn new(data: Option<&pyo3::Bound::<'_,pyo3::types::PyDict>>) -> pyo3::PyResult<Self> {

                        match data {
                            Some(x) => Ok(Self(pythonize::depythonize::<BinaryRecord<[<$name BinaryRecord>]>>(x)?)),
                            None => Err(pyo3::PyErr::new::<pyo3::exceptions::PyTypeError, _>("empty data"))
                        }
                    }

                    fn __repr__(&self) -> String {
                        self.0.to_string()
                    }
            }
        }
    }
}

#[macro_export]
macro_rules! impl_py_parameters {
    ($pyname:literal, $name: ident, $docdir: expr) => {

        use reos::parameters::Parameters;

        paste::paste!{
            #[derive(Clone)]
            #[pyo3::pyclass(name = $pyname)]
            pub struct [<Py $name Parameters> ](pub [<$name Parameters> ]);
            
            #[pyo3::pymethods]
            impl [<Py $name Parameters>] {
                
                #[staticmethod]
                #[pyo3(text_signature = "(pure_records, binary_records = [], **opt)")]
                #[pyo3(signature = (pure_records, binary_records = vec![], **opt))]
                #[doc = include_str!($docdir)]
                fn from_records(py:pyo3::Python, pure_records: Vec<[<Py $name PureRecord>]>, binary_records: Vec<[<Py $name BinaryRecord>]>, opt:Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<Self> {

                    let pure_records = pure_records.into_iter()
                        .map(|r|r.0)
                        .collect();
                    
                    let binary_records = binary_records.into_iter()
                        .map(|r| r.0)
                        .collect();

                    let default: pyo3::Bound<'_, pyo3::types::PyDict> = pyo3::types::PyDict::new(py);
                    let opt: &pyo3::Bound<'_, pyo3::types::PyDict> = opt.unwrap_or(&default);

                    let opt = pythonize::depythonize(opt)?;
                    let res = [<$name Parameters>]::new(pure_records, binary_records, opt);
                    
                    match res {
                        Ok(x) => Ok(Self(x)),
                        Err(e) => Err(pyo3::PyErr::new::<pyo3::exceptions::PyException, _>(e.to_string()))
                    }

                }
                /// Initialize Parameters from JSON files.
                /// 
                /// Parameters
                /// ----------
                /// names : List[str]
                ///     List of component names.
                /// ppath : str
                ///     Path to the pure component JSON file.
                /// bpath : Optional[str]
                ///     Path to the binary interaction JSON file.
                /// opt : kwargs, optional
                ///     Options string.
                ///     Defines the model options values to be used.
                /// 
                /// Returns
                /// -------
                ///     Parameters object initialized from the JSON files.
                ///
                #[staticmethod]
                #[pyo3(text_signature = "(names, ppath, bpath = None)")]
                #[pyo3(signature = (names, ppath, bpath = None, **opt) )]
                fn from_json(py: pyo3::Python, names: Vec<String>, ppath: String, bpath: Option<String>, opt:Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<Self> {

                    let default: pyo3::Bound<'_, pyo3::types::PyDict> = pyo3::types::PyDict::new(py);
                    let opt: &pyo3::Bound<'_, pyo3::types::PyDict> = opt.unwrap_or(&default);
                    let opt = pythonize::depythonize(opt)?;

                    let res = [<$name Parameters>]::from_json(&names, &ppath, bpath.as_ref(), opt);  

                    match res {
                        Ok(x) => Ok(Self(x)),
                        Err(e) => Err(pyo3::PyErr::new::<pyo3::exceptions::PyException, _>(e.to_string()))
                    }

                }


                /// Initialize Parameters from multiple JSON files.
                /// 
                /// Parameters
                /// ----------
                /// - sets : List[List[str]]
                ///     List of component names.
                ///     Each inner list defines a set of components.
                /// 
                /// - ppaths : List[str]
                ///     Paths to the pure component JSON files.
                ///     Its size must be the same of sets.
                /// 
                /// - bpaths : Optional[List[str]]
                ///     Paths to the binary interaction JSON files.
                ///     Cases:
                ///         - None : No binary files provided.
                ///         - [str1] : Finds binary records at `str1` for all sets.
                ///         - [str1, str2, ...] : For each set, finds binary records at the corresponding binary path.
                ///                               If its size is different from sets, an error is raised.
                ///     
                /// - opt : kwargs
                ///     Options string.
                ///     Defines the model options values to be used.
                /// 
                /// Returns
                /// -------
                ///     Parameters object initialized from the JSON files.
                ///
                #[staticmethod]
                #[pyo3(text_signature = "(sets = [names1, names2, ...], ppaths = [ppath1, ppath2, ...], bpaths = None, opt = None)")]
                #[pyo3(signature = (sets, ppaths, bpaths = None, opt = None) )]
                fn from_multiple_jsons(py: pyo3::Python, sets: Vec<Vec<String>>, ppaths: Vec<String>, bpaths: Option<Vec<String>>, opt:Option<&pyo3::Bound<'_, pyo3::types::PyDict>>) -> pyo3::PyResult<Self> {

                    let default: pyo3::Bound<'_, pyo3::types::PyDict> = pyo3::types::PyDict::new(py);
                    let opt: &pyo3::Bound<'_, pyo3::types::PyDict> = opt.unwrap_or(&default);
                    let opt = pythonize::depythonize(opt)?;

                    let res = [<$name Parameters>]::from_multiple_jsons(&sets, &ppaths, bpaths.as_deref(), opt);  

                    match res {
                        Ok(x) => Ok(Self(x)),
                        Err(e) => Err(pyo3::PyErr::new::<pyo3::exceptions::PyException, _>(format!("{:?}", e)))
                    }

                }

                fn __repr__(&self) -> String {
                    self.0.to_string()
                }
            }
        }
    }
}