use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// use crate::py_parameters::PyParameters;

use numpy::{IntoPyArray, PyArray1, PyArrayMethods };
use reeos::{residual::ResidualModel, state::{density_solver::DensityInitialization, State, S}};

use crate::py_eos::PyEquationOfState;




/// A Thermodynamic State 
/// 
/// Units in SI.
/// 
/// Parameters
/// ----------
/// eos : EquationOfState
///     The equation of state to use.
/// 
/// temperature : Kelvin
/// 
/// pressure : Pascal 
/// 
/// x : numpy.ndarray[float]
///     Molar fraction of each component.
/// 
/// Returns
/// -------
/// State : state at (T,P,x) or (T,ρ,x)
#[pyclass(name="State")]
pub struct PyState(pub S<ResidualModel>);

#[pymethods]
impl PyState {

    /// A Thermodynamic State at (T,P,x).
    /// 
    /// Units in SI.
    /// 
    /// Parameters
    /// ----------
    /// eos : EquationOfState
    ///     The equation of state to use.
    /// 
    /// temperature
    /// 
    /// pressure 
    /// 
    /// x : numpy.ndarray[float]
    ///     Molar fraction of each component.
    /// 
    /// Returns
    /// -------
    /// State : state at (T,P,x) 
    #[staticmethod]
    #[pyo3(
    text_signature = "(eos,temperature,pressure,x,density_initialization)"
    )]
    #[pyo3(signature = (eos,temperature,pressure,x,density_initialization))]

    pub fn tpx<'py>
    (
        eos: Bound<'py,PyEquationOfState>,
        temperature:f64,
        pressure:f64,
        x:&Bound<'py, PyArray1<f64>>,
        density_initialization: &str,
    )->PyResult<Self>{
        let eos:PyEquationOfState=eos.extract()?;
        let x = x.to_owned_array();
        let phase = DensityInitialization::from_str(density_initialization);

        if phase.is_err(){
            return Err(PyErr::new::<PyValueError, _>(
                                "`density_initialization` must be 'vapor' or 'liquid'.".to_string(),
                            ))
        }

        let res= State::new_tpx(&eos.0, temperature, pressure, x, phase.unwrap());
        match res{
        Ok(state)=> Ok(PyState(state)),
        Err(e)=> {
            Err(PyErr::new::<PyValueError, _>(
                e.to_string(),
            ))
            }
        }
    }
    /// A Thermodynamic State at (T,ρ,x)
    /// 
    /// Units in SI.
    /// 
    /// Parameters
    /// ----------
    /// eos : EquationOfState
    ///     The equation of state to use.
    /// 
    /// temperature
    /// 
    /// density 
    /// 
    /// x : numpy.ndarray[float]
    ///     Molar fraction of each component.
    /// 
    /// Returns
    /// -------
    /// State : state at (T,ρ,x)
    #[staticmethod]
    #[pyo3(
    text_signature = "(eos,temperature,density,x)"
    )]
    #[pyo3(signature = (eos,temperature,density,x))]

    pub fn trx<'py>
    (
        eos: Bound<'py,PyEquationOfState>,
        temperature:f64,
        density:f64,
        x:&Bound<'py, PyArray1<f64>>,
    )->PyResult<Self>{
        let eos:PyEquationOfState=eos.extract()?;
        let x = x.to_owned_array();

        let res= State::new_trx(&eos.0, temperature, density, x);
        match res{
        Ok(state)=> Ok(PyState(state)),
        Err(e)=> {
            Err(PyErr::new::<PyValueError, _>(
                e.to_string(),
            ))
            }
        }
    }

    /// Return logarithmic fugacity coefficient.
    ///
    /// Returns
    /// -------
    /// numpy.ndarray
    pub fn ln_phi<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {

        let state=&self.0;
        state.lnphi().unwrap().into_pyarray(py)

    }

    /// bmix.
    ///
    /// Returns
    /// -------
    /// float
    pub fn bmix(&self) -> f64 {
        
        self.0.bmix()

        }

    /// Return Pressure.
    ///
    /// Returns
    /// -------
    /// float
    pub fn pressure(&self) -> f64 {
        self.0.pressure().unwrap()
    }
    /// Return volume.
    ///
    /// Returns
    /// -------
    /// float
    pub fn volume(&self) -> f64 {
        1./self.0.rho
    }



    // pub fn tpd<'py>(&self,xphase:&str,xguess:&Bound<'py, PyArray1<f64>>)->PyResult<f64>{

    //     let xguess = xguess.to_owned_array();

    //     let xphase = 
    //     if let Ok(val)=DensityInitialization::from_str(xphase){
    //         val
    //     }else{
    //         return Err(PyErr::new::<PyValueError, _>(
    //             "xphase must be 'liquid' or 'vapor' "))
    //     };
        
    //     let state=&self.0;
    //     match state{
    //         StateCPA::SRK(eos)=>{
    //             match eos.tpd(xphase, xguess) {
    //             Ok(dg)=>Ok(dg),

    //             Err(e)=>Err(PyErr::new::<PyValueError, _>(
    //                 e.to_string()))
    //             } 

    //         }
    //         StateCPA::PR(eos)=>{
    //             match eos.tpd(xphase, xguess) {
    //             Ok(dg)=>Ok(dg),
                
    //             Err(e)=>Err(PyErr::new::<PyValueError, _>(
    //                 e.to_string()))
    //             }
    //         }
 
    //     }
    // }

    // pub fn water_saturation(&self)->PyResult<f64>
    // {
    //     let state=&self.0;
    //     match state{
    //         StateCPA::SRK(eos)=>{
    //         match eos.water_saturation(){
    //             Ok(sat)=>Ok(sat),

    //             Err(e)=>Err(PyErr::new::<PyValueError, _>(
    //                 e.to_string()))
    //         }

    //         }
    //         StateCPA::PR(eos)=>{
    //             match eos.water_saturation(){
    //                 Ok(sat)=>Ok(sat),

    //                 Err(e)=>Err(PyErr::new::<PyValueError, _>(
    //                     e.to_string()))
    //             }
    //         }
    //     }
    // }

}
