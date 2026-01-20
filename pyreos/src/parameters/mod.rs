//! Parameters API for Equation Of State models.
//! 
//! Provides macros to generate Python classes 
//! for Parameters, PureRecord and BinaryRecord, since
//! this objects depend on generic model types from `reos` crate.
//! 
//! Current macros:
//! - `impl_parameters!`
//! - `impl_pure_record!`
//! - `impl_binary_record!`
//! 
pub mod records;

/// Macro to generate concrete Python wrapper-like structs for a generic Parameters type, 
/// given the corresponding PureRecord and BinaryRecord.
#[macro_export]
macro_rules! impl_parameters {
    ($name: ident, $ptype: ident, $pr_type: ident, $br_type: ident, $pyname: expr, $options_type: ident, $docdir: expr, ) => {

        #[pyclass(name = $pyname)]
        #[doc = stringify!($pyname)]
        pub struct $name (
            pub $ptype
        );
        
        impl Into<$ptype> for $name {
            fn into(self) -> $ptype{
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
            
            #[doc = include_str!($docdir)]
            #[staticmethod]
            #[pyo3(text_signature = "(pure_records, binary_records = [], opt = None)")]
            #[pyo3(signature = (pure_records, binary_records = vec![], opt = None) )]
            pub fn from_records(pure_records: Vec<$pr_type>, binary_records:Vec<$br_type>, opt: Option<&str>)-> PyResult<Self>{

                let pure_records = pure_records
                .into_iter()
                .map(|r|r.0)
                .collect();

                let binary_records = binary_records
                .into_iter()
                .map(|r| r.0)
                .collect();
                
                if let Some(opt) = opt {

                    match $options_type::try_from(opt) {
                        Ok(_) => {
                            let options:$options_type = $options_type::try_from(opt).unwrap();
                            let new = $ptype::new(pure_records, binary_records, options);

                            return Ok(Self(new));
                        },
                        Err(e) => {
                            return Err(crate::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                                e.to_string(),
                            ));
                        }
                    }
   
                } else {
                    let new = $ptype::new(pure_records, binary_records, $options_type::default());

                    return Ok(Self(new));
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
            /// opt : Optional[str]
            ///     Options string.
            ///     Defines the model type to be used.
            /// 
            /// Returns
            /// -------
            ///     Parameters object initialized from the JSON files.
            ///
            #[staticmethod]
            #[pyo3(text_signature = "(names, ppath, bpath = None, opt = None)")]
            #[pyo3(signature = (names, ppath, bpath = None, opt = None) )]
            fn from_json(names: Vec<String>, ppath: String, bpath: Option<String>, opt: Option<&str>) -> PyResult<Self> {

                
                if let Some(opt) = opt {

                    match $options_type::try_from(opt) {
                        Ok(_) => {
                            let options:$options_type = $options_type::try_from(opt).unwrap();
                            let result = $ptype::from_json(&names, &ppath, bpath.as_ref(), options);

                            if let Err(e) = result {
                                return Err(crate::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                                    e.to_string(),
                                ));
                            }

                            return Ok(Self(result.unwrap()));
                        },
                        Err(e) => {
                            return Err(crate::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                                e.to_string(),
                            ));
                        }
                    }
   
                } else {

                    let result = $ptype::from_json(&names, &ppath, bpath.as_ref(), $options_type::default());

                    if let Err(e) = result {
                        return Err(crate::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                            e.to_string(),
                        ));
                    }

                    return Ok(Self(result.unwrap()));
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
            /// - opt : Optional[str]
            ///     Options string.
            ///     Defines the model type to be used.
            /// 
            /// Returns
            /// -------
            ///     Parameters object initialized from the JSON files.
            ///
            #[staticmethod]
            #[pyo3(text_signature = "(sets = [names1, names2, ...], ppaths = [ppath1, ppath2, ...], bpaths = None, opt = None)")]
            #[pyo3(signature = (sets, ppaths, bpaths = None, opt = None) )]
            fn from_multiple_json(sets: Vec<Vec<String>>, ppaths: Vec<String>, bpaths: Option<Vec<String>>, opt: Option<&str>) -> PyResult<Self> {

                if let Some(opt) = opt {

                    match $options_type::try_from(opt) {
                        Ok(_) => {

                            let options:$options_type = $options_type::try_from(opt).unwrap();
                            let result = $ptype::from_multiple_jsons(&sets, &ppaths, bpaths.as_deref(), options);

                            if let Err(e) = result {
                                return Err(crate::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                                    format!("{:?}", e)
                                ));
                            }

                            Ok(Self(result.unwrap()))
                        },
                        Err(e) => {

                            Err(crate::PyErr::new::<pyo3::exceptions::PyValueError, _>(
                                e.to_string(),
                            ))
                        }
                    }
   
                } else {
                    // let sets: Vec<Vec<&str>> = sets.iter().map(|v| v.iter().map(|s| s.as_str()).collect()).collect();
                    let result = $ptype::from_multiple_jsons(&sets, &ppaths, bpaths.as_deref(), $options_type::default());
                    
                    if let Err(e) = result {
                        return Err(crate::PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                            format!("{:?}", e)
                        ));
                    }
                    Ok(Self(result.unwrap()))
                }


            }
            fn __repr__(&self) -> String {
                self.0.to_string()
            }
            
        }
    }
}



// impl_parameters!(PyPr76Parameters,PR76Parameters,PyCubicPureRecord,PyCubicBinaryRecord,"PR76Parameters");
// impl_parameters!(PyPr78Parameters,PR78Parameters,PyCubicPureRecord,PyCubicBinaryRecord,"PR78Parameters");


// #[pyclass(name = "SrkParameters")]
// // pub struct Foo (pub CubicParameters<SRK>);


// #[pymethods]
// impl Foo {

//     #[staticmethod]
//     fn from_records(records: Vec<PyCubicPureRecord>, binary:Vec<PyCubicBinaryRecord>){


//         let pure_records = records
//         .into_iter()
//         .map(|r|r.0)
//         .collect();

//         let binary_records= binary
//         .into_iter()
//         .map(|r| r.0)
//         .collect();

//         let new = CubicParameters::<SRK>::new(pure_records, binary_records);

//         serde_json::to_string(&new).unwrap();
//     }
// }
// #[pyclass(name = "CPAParameters")]
// pub struct PySrkCpaParameters(pub CubicParameters<SRK>, pub AssociativeParameters);


// impl Parameters<PyCPARecord,_> for PySrkCpaParameters {

//     fn raw(records: Vec<PyCPARecord>, binary: Vec<BinaryParameter<B>>, properties: Option<reos::parameters::Properties>)->Self {
        
//         unimplemented!()
//     }

//     fn get_properties(&self)->&reos::parameters::Properties {
        
//         unimplemented!()
//     }
// }

// #[pymethods]
// impl PyCpaParameters {


// /// Parameters
// /// ----------
// ///     cubic: List[CubicRecord],
// ///     assoc: List[AssociationRecord]
// /// 
// /// Returns
// /// -------
//     /// Return CPA Parameters from cubic & assoc records.
//     #[staticmethod]
//     #[pyo3(
//         signature = (cubic,assoc),
//         text_signature = "(cubic, assoc)"
//     )]
//     pub fn from_records(cubic:Vec<PyCubicRecord>,assoc:Vec<PyAssocRecord>)->Self{

//         Self 
//         {asc: ASCParameters::from_records(assoc.iter().map(|u|u.0.clone()).collect()), 
//          cub: CubicParameters::from_records(cubic.iter().map(|u|u.0.clone()).collect())  }

//     }
//     #[pyo3(
//     signature = (i,j,kij_a=0.0,kij_b=0.0),
//     text_signature = "(i,j,kij_a=0.0,kij_b=0.0)",
//     )]

//     pub fn set_cubic_binary(&mut self,i:usize,j:usize,kij_a:f64,kij_b:f64){
//         self.cub.set_binary(i, j, Some(kij_a), kij_b);
//     }
    
//     #[pyo3(
//     signature = (i,j,rule,eps=None,beta=None),
//     text_signature = "(i,j,rule,eps=None,beta=None)",
//     )]
//     pub fn set_assoc_binary(&mut self,i:usize,j:usize,rule:String,eps:Option<f64>,beta:Option<f64>)->PyResult<()>{

//         match AssociationRule::try_from(rule){
//             Ok(rule)=>{
//                 self.asc.set_binary_interaction(i, j, rule, eps, beta);
//                 Ok(())
//             }
//         Err(e)=> {
//             Err(PyErr::new::<PyValueError, _>(
//                 e.to_string(),
//             ))}
//         }
//     }

//     #[pyo3(
//     signature = (sitej,sitel,eps,beta),
//     text_signature = "(sitej,sitel,eps,beta)",
//     )]
//     pub fn change_sites_p(&mut self,sitej:usize,sitel:usize,eps:f64,beta:f64)->PyResult<()>{

//         self.asc.change_sites_parameters(sitej, sitel, eps, beta);
        
//         return Ok(());

//         // Err(e)=> {
//         //     Err(PyErr::new::<PyValueError, _>(
//         //         e.to_string(),
//         //     ))}
//         // }
//     }
//     pub fn as_string(&self)->String{

//         let asc =format!("{}",&self.asc);
//         let cub =format!("{}",&self.cub);
//         cub + &asc

//     }
// }



// #[derive(Clone)]
// #[pyclass(name = "CubicParameters")]
// pub struct PyCubicParameters(
//     pub CubicParameters);

// #[pymethods]
// impl PyCubicParameters {

// /// Parameters
// /// ----------
// /// 
// /// List[CubicRecord],
// /// 
// /// Returns
// /// -------
//     #[staticmethod]
//     #[pyo3(
//         signature = (parameters),
//         text_signature = "(parameters)"
//     )]
//     pub fn from_records(parameters:Vec<PyCubicRecord>)->Self{

//         let p =CubicParameters::from_records(parameters.iter().map(|u|u.0.clone()).collect());
        
//         Self(
//             p
//         )

//     }

//     pub fn as_string(&self)->String{

//         let cub =format!("{}",&self.0);
//         cub 

//     }

//     #[pyo3(
//     signature = (i,j,kij_a=0.0,kij_b=0.0),
//     text_signature = "(i,j,kij_a=0.0,kij_b=0.0)",
//     )]

//     pub fn set_cubic_binary(&mut self,i:usize,j:usize,kij_a:f64,kij_b:f64){
//         self.0.set_binary(i, j, Some(kij_a), kij_b);
//     }
// }
// #[cfg(test)]
// pub mod tests{
//     use std::sync::Arc;

    
//     use reos::models::cpa::parameters::water_co2;
//     #[test]
//     fn show_sites(){
//         let eos=water_co2();
//         let sites=&eos.residual.assoc.parameters.f;

//         // println!("Sites = {}" ,sites);
//         // let eos=water_octane_acetic_acid();
//         // let sites=&eos.residual.assoc.parameters.f;

//         // println!("Sites = {}" ,sites);
//     }
// }