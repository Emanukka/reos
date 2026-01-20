//! Directory for python models API
//!
//! Each submodule is seen at python as submodule of `reos`.
//! Therefore, each submodule uses the macros defined in this crate 
//! to create the corresponding PureRecord, BinaryRecord, Parameters 
//! and EquationOfState python objects.
//! 
pub mod cubic;
pub mod cpa;