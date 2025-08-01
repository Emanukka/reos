use numpy::{IntoPyArray, PyArray, PyArray1, PyArray2, PyArrayMethods};
use pyo3::{exceptions::PyTypeError, pyclass, pymethods, Bound, Py, PyErr, PyResult, Python};
use reos::{models::associative::Associative};

use crate::py_eos::py_residual::ResidualModel;

#[derive(Clone)]

#[pyclass(name="Association")]
pub struct  PyAssociation(pub Associative);

#[pymethods]

impl PyAssociation {
    
    fn non_bonded_sites<'py>(&self,t:f64,rho:f64,x: &Bound<'py, PyArray1<f64>>)->Bound<'py, PyArray1<f64>>{

        let assoc=&self.0;
        let token=x.py();

        let x=x.to_owned_array();
        let x=assoc.X_tan(t, rho, &x);
        x.unwrap().into_pyarray(token)

    }
    
    fn get_tmat<'py>(&self,py: Python<'py>)->Bound<'py, PyArray2<f64>>{

        let assoc=&self.0;
        let tmat=assoc.parameters.tmat.clone();
        tmat.into_pyarray(py)

    }
    // fn get_fmap<'py>(&self,py: Python<'py>)->Bound<'py, PyArray1<f64>>{

    //     let assoc=&self.0;
    //     // assoc.parameters.f.clone().into_pyarray(py)

    // }

    fn get_nassoc<'py>(&self,py: Python<'py>)->Bound<'py, PyArray1<usize>>{

        let assoc=&self.0;
        assoc.parameters.nassoc.clone().into_pyarray(py)


    }
    
}


