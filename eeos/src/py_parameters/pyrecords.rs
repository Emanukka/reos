use pyo3::{pyclass, pyclass_init, pymethods};
use reeos::parameters::{association::AssociationPureRecord, cubic::CubicPureRecord};


#[derive(Clone)]
#[pyclass(name = "CubicRecord")]
pub struct PyCubicRecord(pub CubicPureRecord);

#[derive(Clone)]
#[pyclass(name = "AssociationRecord")]
pub struct PyAssocRecord(pub AssociationPureRecord);


#[pymethods]
impl PyCubicRecord {
    
    #[new]
    pub fn new(a0:f64,b:f64,c1:f64,tc:f64)->Self{
        PyCubicRecord(
        CubicPureRecord::new(a0, b, c1, tc)
        )
    }

}
#[pymethods]
impl PyAssocRecord {
    
    #[staticmethod]
    pub fn associative(eps:f64,beta:f64,b:f64,na:usize,nb:usize,nc:usize)->Self{
        PyAssocRecord(
        AssociationPureRecord::associative(eps, beta, [na,nb,nc], b)
        )
    }
    #[staticmethod]
    pub fn solvate(b:f64,na:usize,nb:usize,nc:usize)->Self{
        PyAssocRecord(
        AssociationPureRecord::solvate( [na,nb,nc], b)
        )
    }
    #[staticmethod]
    pub fn inert(b:f64)->Self{
        PyAssocRecord(
        AssociationPureRecord::inert(b)
        )
    }
}

