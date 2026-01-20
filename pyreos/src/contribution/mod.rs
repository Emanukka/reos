//! Contribution API for Python.
//!
//! Defines the `PyContribution` enum, which wraps concrete residual contribution models from `reos`.
//! The choice for enum instead of diferents structs (such as made for Parameters and Records) is
//! that would be necessary to create a equal number of PyState wrappers for each contribution model.
//! 
//! # The `impl_contribution!` macro 
//! 
//! This macro is used to implement the `Residual` trait for all variants of `PyContribution`
//! enabling the `EquationOfState<PyContribution>` be a concrete type, such that every call 
//! to a function belonging the `Residual` trait is properly dispatched through a pattern matching in each
//! `PyContribution` variant.
//! 
//! 
//! Example:
//! ```rust
//! impl_eos!(SRK, SRKContribution, PySrkParameters, "srk");
//! ```
//! 

use reos::residual::Residual;
use reos::models::cpa::rdf::{CS,Kontogeorgis};
use reos::models::cubic::{Cubic};
use reos::models::cpa::CPA;


use reos::Array1;

/// Wrapper enum for the current contribution models.
/// This object won't be exposed to python.
// #[derive(Clone)]
pub enum PyContribution{
    
    Cubic(Cubic),
    SCPA(CPA<Kontogeorgis>),
    CPA(CPA<CS>),
    // CPA(CPA<SRK,Kontogeorgis>),

} 

    

macro_rules! impl_contribution {

    (
        $enum:ident, { $( $variant:ident ),+ $(,)? }
    ) => {
        impl Residual for $enum {
            
            fn molar_weight(&self) -> &Array1<f64> {
                match self {
                    $( $enum::$variant(model) => model.molar_weight(), )+
                }
            }
            fn r_pressure(&self,t:f64,d:f64,x:&Array1<f64>) -> f64 {
                match self {
                    $( $enum::$variant(model) => model.r_pressure(t,d,x), )+
                }
            }

            fn r_helmholtz(&self,t:f64,d:f64,x:&Array1<f64>) -> f64 {
                match self {
                    $( $enum::$variant(model) => model.r_helmholtz(t,d,x), )+
                }
            }
            
            fn r_chemical_potential(&self,t:f64,d:f64,x:&Array1<f64>) -> Array1<f64> {
                match self {
                    $( $enum::$variant(model) => model.r_chemical_potential(t,d,x), )+
                }
            }

            fn r_entropy(&self,t:f64,d:f64,x:&Array1<f64>) -> f64 {
                match self {
                    $( $enum::$variant(model) => model.r_entropy(t,d,x), )+
                }
            }

            fn max_density(&self, x:&Array1<f64>) -> f64 {
                match self {
                    $( $enum::$variant(model) => model.max_density(x), )+
                }
            }

            fn components(&self) -> usize {
                match self {
                    $( $enum::$variant(model) => model.components(), )+
                }
            }
            
            fn get_properties(&self) -> &reos::parameters::Properties {
                match self {
                    $( $enum::$variant(model) => model.get_properties(), )+
                }
            }
        }
    };
}


impl_contribution!( 
    PyContribution, 
    {SCPA,CPA,Cubic}
);

