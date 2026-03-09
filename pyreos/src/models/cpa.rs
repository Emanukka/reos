use crate::{contribution::PyContribution, impl_eos, impl_py_binary_record, impl_py_parameters, impl_py_pure_record};

use numpy::{PyArray1, PyArrayMethods, ToPyArray};
use pyo3::{Bound, PyAny, PyErr, PyResult, Python, types::{IntoPyDict, PyDict, PyFloat, PyInt, PyList, PyString}};
use reos::models::cpa::{ 
    CPA, parameters::{CPABinaryRecord, CPAParameters, CPAPureRecord}
};

// use reos::parameters::{Parameters,PureRecord,BinaryRecord};

impl_py_pure_record!("CPAPureRecord", CPA, "../../docs/cpa/pr.md");

impl_py_binary_record!("CPABinaryRecord", CPA, "../../docs/cpa/br.md");

impl_py_parameters!("CPAParameters", CPA ,"../../docs/cpa/parameters.md");


impl_eos!(CPA, "../../docs/cpa/eos.md");


#[pyo3::pymethods]
impl PyCPAParameters {

    /// Returns a list of sites
    /// 
    /// The content of each site `j` is inside a dictionary, whose the fields are:
    /// 
    /// - `type`: type of site j
    /// - `idx`: index of  site j
    /// - `owner`: owner of site j
    /// - `mul`: multiplicity of site j
    /// - `epsilon`: association binding energy of site j
    /// - `kappa`: association binding volume of site j
    /// 
    /// Returns
    /// -------
    /// List[Dict]
    pub fn get_sites<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyList>> {

        let sites = &self.0.assoc.sites;
        let mut v1: Vec<Bound<'_, PyDict>> = vec![];
        
        for j in 0..sites.len() {

            let sj = &sites[j];
            let py_idx = PyInt::new(py, sj.idx);
            let py_owner = PyInt::new(py, sj.owner);
            let py_mul = PyFloat::new(py, sj.mul);
            let py_eps = PyFloat::new(py, sj.epsilon);
            let py_kappa = PyFloat::new(py, sj.kappa);
            let py_type = PyString::new(py,sj.typ.to_string().as_str() );

            let arr: Vec<(&str, &Bound<'_, PyAny>)> = vec![
                (("type", py_type.as_any())),
                (("idx"), py_idx.as_any()),
                (("owner", py_owner.as_any())),
                (("mul", py_mul.as_any())),
                (("epsilon", py_eps.as_any())),
                (("kappa", py_kappa.as_any()))
            ];
            let py_dct = arr.into_py_dict(py)?;
            v1.push(py_dct);
        }

        PyList::new(py, v1)
        

    }
}
#[pyo3::pymethods]
impl crate::eos::PyEquationOfState {

    /// Returns association variables
    /// 
    /// The variables are returned inside a dictionary, where the fields are:
    /// 
    /// - `X`:  `Xⱼ` unbonded fraction of site j
    /// - `Delta`: `Δⱼₗ` matrix of association strength between sites j and l 
    /// - `K`: matrix of association constants of sites j and l, defined as `Kⱼₗ = ρmⱼmₗΔⱼₗ`
    /// - `m`: `mⱼ` sites mole fraction of site j, such that `mⱼ = xᵢ₍ⱼ₎`
    /// 
    /// Parameters
    /// ----------
    /// t : float
    ///     Temperature [K]
    /// d : float
    ///     Density [mol/m3]
    /// x : numpy.ndarray[float]
    ///     Mole fractions
    /// Returns
    /// -------
    /// Dict[str, np.ndarray]

    #[pyo3(text_signature = "(temperature, density, x)")]
    pub fn get_assoc_calcs<'py>(&self, temperature:f64, density:f64, x: &Bound<'py,PyArray1<f64>>)-> PyResult<Bound<'py, PyDict>> {

        let py = x.py();
        let contribution = self.0.residual();
        let x = x.to_owned_array();



        match &contribution{
            PyContribution::CPA(inner)=>{
                
                // allocate memory on python heap
                let k = inner.assoc_consts(temperature, density, &x);
                let u = inner.assoc.assoc.unbonded_sites_fraction(&x, &k).to_pyarray(py);
                let delta= inner.delta_assoc(density, &x, &k).to_pyarray(py);
                let k = k.to_pyarray(py);

                let m = inner.assoc.assoc.sites_mole_frac(&x).to_pyarray(py);

                let arr: Vec<(&str, &Bound<'_, PyAny>)> = vec![
                    (("Delta", delta.as_any())),
                    (("X"), u.as_any()),
                    (("K", k.as_any())),
                    (("m", m.as_any())),
                ];

                arr.into_py_dict(py)

            }

            _ =>  Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Only eos with `Association` contribution have this method",)
                )
        }

    } 
}
