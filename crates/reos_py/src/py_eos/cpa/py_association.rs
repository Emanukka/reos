use numpy::{ IntoPyArray, PyArray1, PyArray2};
use pyo3::{ pyclass, pymethods, types::PyTuple, Bound, Python};
use reos::models::cpa::associative::Associative;



#[derive(Clone)]

#[pyclass(name="Association")]
pub struct  PyAssociation(pub Associative);

#[pymethods]

impl PyAssociation {
    

    // /// Get the fraction of non-bonded sites X.
    // /// 
    // /// Parameters
    // /// ----------
    // /// 
    // /// t: f64
    // ///    Temperature
    // /// 
    // /// rho: f64
    // ///    Molar Density
    // /// 
    // /// x : numpy.ndarray[float]
    // ///     Molar fraction of each component.
    // /// 
    // /// Returns
    // /// -------
    // /// Fraction of non-bonded sites X as np.ndarray
    // #[pyo3(
    // text_signature = "(t,rho,x)"
    // )]
    // #[pyo3(signature = (t,rho,x))]
    // fn non_bonded_sites<'py>(&self,t:f64,rho:f64,x: &Bound<'py, PyArray1<f64>>)->Bound<'py, PyArray1<f64>>{

    //     let assoc=&self.0;
    //     let token=x.py();

    //     let x=x.to_owned_array();
    //     let x=assoc.X_tan(t, rho, &x);
    //     x.unwrap().into_pyarray(token)

    // }
    
    /// Get the transformation matrix T.
    /// 
    /// Row's lenght is the number of associative componentes 'n',
    /// Column's lenght is the number of distinctive sites 'NS'
    /// 
    /// Each element indicates if component i (row) has site k (column)
    /// 
    /// Returns
    /// -------
    /// T matrix
    fn get_tmat<'py>(&self,py: Python<'py>)->Bound<'py, PyArray2<f64>>{

        let assoc=&self.0;
        let tmat=assoc.parameters.tmat.clone();
        tmat.into_pyarray(py)

    }
    /// Get the map F.
    /// 
    /// Map's size is the number of distinctive sites (NS) in mixture.
    /// 
    /// At the k position in F[k], it's returned a tuple (j,i):
    /// 
    /// - j is the type of the site (A(0),B(1) or C (2)),
    /// - i is the owner of the site.
    /// 
    /// The map is useful to determine in the X which is the site's owner at k position.
    /// 
    /// Returns
    /// -------
    /// F array
    fn get_fmap<'py>(&self,py: Python<'py>)->Bound<'py,PyTuple>{

        let assoc=&self.0;
        let sites:Vec<(usize,usize)>=assoc.
        parameters.
        f.clone().
        iter().map(|s|s.clone().into())
        .collect(); 

        let tup=PyTuple::new(py,sites).unwrap();
        tup
    }

    /// Get array (lenght 'n') with indices of all associative components in the mixture.
    /// 
    /// Returns
    /// -------
    /// nassoc array
    fn get_nassoc<'py>(&self,py: Python<'py>)->Bound<'py, PyArray1<usize>>{

        let assoc=&self.0;
        assoc.parameters.nassoc.clone().into_pyarray(py)


    }
    
}


