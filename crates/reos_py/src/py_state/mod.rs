use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// use crate::py_parameters::PyParameters;

use numpy::{IntoPyArray, PyArray1, PyArrayMethods };
use reos::{state::{density_solver::DensityInitialization, State, S}};

use crate::py_eos::{PyEquationOfState,py_residual::ResidualModel};



/// A Thermodynamic State 
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
    /// State at (T,P,x) 
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
    /// State at (T,ρ,x)
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

    /// Return the Hard-Shere volume 'b' of the mixture.
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
        self.0.p
    }
    /// Return Temperature.
    ///
    /// Returns
    /// -------
    /// float
    pub fn temperature(&self) -> f64 {
        self.0.t
    }
    /// Return volume.
    ///
    /// Returns
    /// -------
    /// float
    pub fn volume(&self) -> f64 {
        1./self.0.rho
    }
    /// Return density.
    ///
    /// Returns
    /// -------
    /// float
    pub fn density(&self) -> f64 {
        self.0.rho
    }

    /// Return the composition z of the State.
    ///
    /// Returns
    /// -------
    /// np.ndarray[float]
    pub fn composition<'py>(&self,py:Python<'py>) -> Bound<'py, PyArray1<f64>> {

        self.0.x.clone().into_pyarray(py)
    }

    /// Return the compressibility factor Z.
    ///
    /// Returns
    /// -------
    /// float
    pub fn compressibility(&self) -> f64 {

        self.0.eos.compressibility(self.0.t, self.0.rho, &self.0.x).unwrap()
    }


    /// Return the mininum of TPD and the incipient phase state.
    ///
    /// Parameters
    /// -------
    /// 
    /// xphase: 'liquid','vapor'
    ///     incipiente phase 
    /// 
    /// xguess: numpy.ndarray[float] 
    ///     initial guess for the xphase's composition 
    ///     
    /// 
    /// Returns
    /// -------
    /// ΔG formation of the incipient phase from mother phase State at (T,P,z) condition.
    /// and incipient phase State at (T,P,x).
    #[pyo3(
    text_signature = "(xphase,xguess,tol=1e-8,it_max=100)"
    )]
    #[pyo3(signature = (xphase,xguess,tol=1e-8,it_max=100))]

    pub fn min_tpd<'py>(&self,xphase:&str,xguess:&Bound<'py, PyArray1<f64>>,tol:Option<f64>,it_max:Option<i32>)->PyResult<(f64,PyState)>{

        let state=&self.0;
        let xphase = DensityInitialization::from_str(xphase);
        let xguess=xguess.to_owned_array();
        if xphase.is_err(){
            return Err(PyErr::new::<PyValueError, _>(
                                "`density_initialization` must be 'vapor' or 'liquid'.".to_string(),
                            ))
        }

        let (dg,new_phase)=state.min_tpd(xphase.unwrap(), xguess, tol, it_max).unwrap();

        Ok(
        (dg,PyState(new_phase))
        )
    }


}
