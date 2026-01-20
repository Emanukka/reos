pub mod tpd;

// use numpy::{PyArray1, PyArrayMethods, ToPyArray};
// use pyo3::{exceptions::PyValueError, pyclass, pymethods, types::PyAnyMethods, Bound, PyResult, Python};
// use reos::{phase_equilibrium::{Antoine, PhaseEquilibrium}, state::{density_solver::DensityInitialization, eos::EosError}};
// use pyo3::PyErr;
// use crate::py_eos::{py_residual::ResidualModel, PyEquationOfState};

// #[pyclass(name = "Antoine")]
// pub struct PyAntoine(Antoine);
// #[pyclass(name = "PhaseEquilibrium")]
// pub struct PyPhaseEquilibrium(pub PhaseEquilibrium<ResidualModel>);

// #[pymethods]
// impl PyPhaseEquilibrium {

//     #[new]
//     pub fn new<'py> (eos:Bound<'py,PyEquationOfState>)->PyResult<Self> {

//         let eos:PyEquationOfState=eos.extract()?;
//         Ok(
//         PyPhaseEquilibrium(
//         PhaseEquilibrium::new(&eos.0, None)))
//     }


//     #[pyo3(
//     signature = (t,x,tol_p=1e-7,tol_y=1e-6),
//     text_signature = "(t,x,tol_p=1e-7,tol_y=1e-6)",
//     )]
//     pub fn bbpy<'py>(&self,t:f64,x:&Bound<'py, PyArray1<f64>>,tol_p:f64,tol_y:f64)->PyResult<(f64,Vec<f64>)>{
//         let x = x.to_owned_array();

//         match self.0.bbpy(t, x, Some(tol_p), Some(tol_y)){

//             Ok(r)=>{Ok((r.0,r.1.to_vec()))},
//             Err(e)=> {
//                 Err(PyErr::new::<PyValueError, _>(
//                     e.to_string(),
//                 ))}
//         }
//     }
//     /// Parameters
//     /// ----------
//     ///     t,
//     ///     p,
//     ///     z: MotherPhase Composition,
//     ///     x: DaughterPhase ('liquid' or 'vapor'),
//     ///     xguess: DaughterPhase Guess Composition
//     /// Returns
//     /// -------
//     /// ΔG formation of the incipient phase from mother phase at (T,P,z) condition.
//     /// and incipient phase composition x.
//     #[pyo3(
//     signature = (t,p,z,incipient_phase,xguess,tol,it_max),
//     text_signature = "(t,p,z,incipient_phase,xguess,tol=1e-8,it_max=100)",
//     )]
//     pub fn tpd<'py>(
//         &self,
//         t:f64,
//         p:f64,
//         z:&Bound<'py, PyArray1<f64>>,
//         incipient_phase:String,
//         xguess:&Bound<'py, PyArray1<f64>>,
//         tol:f64,
//         it_max:i32)->PyResult<(f64,Bound<'py, PyArray1<f64>>)>{

//             let xphase=DensityInitialization::from_str(&incipient_phase);
//             if xphase.is_err(){
//                 return Err(PyErr::new::<PyValueError, _>(
//                                     "`density_initialization` must be 'vapor' or 'liquid'.".to_string(),
//                                 ))
//             }

//             let py=z.py();
//             let xguess = xguess.to_owned_array();
//             let z = z.to_owned_array();


//             match self.0.tpd(t, p, z, xphase.unwrap(), xguess, Some(tol), Some(it_max)){
//                 Ok(r)=>{Ok(
//                     (r.0,PyArray1::from_array(py, &r.1))
//                 )},
//                 Err(e)=> {
//                     Err(PyErr::new::<PyValueError, _>(
//                         e.to_string(),
//                     ))}
//                 }

        
//         }
//     // fn set_antoine<'py>(&mut self,antoine:Bound<'py,PyAntoine>){

//     //     self.0.set_antoine(antoine);
//     // }
// }