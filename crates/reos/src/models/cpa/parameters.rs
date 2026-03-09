use std::collections::HashMap;
use std::str::FromStr;
use serde::{Deserialize, Serialize};

use crate::models::associative::parameters::*;

use crate::models::cpa::rdf::{RDF, RDFcpa, RDFmodel, RDFmodelOption};
use crate::models::cubic::models::CubicModelOption;
use crate::models::cubic::options::CubicOptions;
use crate::models::cubic::parameters::{CubicBinaryRecord, CubicParameters, CubicPureRecord};
use crate::parameters::Parameters;


#[derive(Serialize,Deserialize, Debug, Clone)]
pub struct CPAPureRecord{
    #[serde(flatten)]
    pub c:CubicPureRecord,
    #[serde(flatten)]
    pub a:AssociationPureRecord

}

impl CPAPureRecord{

    pub fn new(c:CubicPureRecord,a:AssociationPureRecord)->Self{
        Self { c, a }
    }

}


#[derive(Clone)]
pub struct CPAParameters{
    pub cubic: CubicParameters, 
    pub assoc: AssociativeParameters,
    pub rdf: RDFcpa,
}

#[derive(Serialize,Deserialize,Clone,Debug)]
pub struct CPABinaryRecord{
    #[serde(flatten)]
    pub association: Option<AssociationBinaryRecord>,
    #[serde(flatten)]
    pub cubic: Option<CubicBinaryRecord>

}
#[derive(Serialize,Deserialize)]
pub struct CPAOptions{

    pub rdf_model: RDFmodelOption,
    #[serde(flatten)]
    pub cubic_options: CubicOptions
}


impl CPABinaryRecord{

    pub fn new(c:Option<CubicBinaryRecord>,a:Option<AssociationBinaryRecord>)->Self{
        Self { cubic:c, association:a }
    }

    pub fn full(c:CubicBinaryRecord, a:AssociationBinaryRecord)->Self{
        
        Self::new(Some(c), Some(a))
    }

    pub fn only_c(c:CubicBinaryRecord)->Self{
        
        Self::new(Some(c), None)
    }

    pub fn only_a(a:AssociationBinaryRecord)->Self{
        
        Self::new(None, Some(a))
    }

}





impl Parameters for CPAParameters {

    type Pure = CPAPureRecord;
    type Binary = CPABinaryRecord;
    type Options = CPAOptions;

    fn from_raw(pure:Vec<Self::Pure>, binary: crate::parameters::BinaryMap<Self::Binary>, properties: Option<crate::parameters::Properties>, opt: Self::Options) -> Result<Self, Box<dyn std::error::Error>> {
        
        let n = pure.len();

        let mut c_pure= Vec::with_capacity(n);
        let mut a_pure = Vec::with_capacity(n);

        let mut c_binary = HashMap::new();
        let mut a_binary = HashMap::new();
        
        for (key, b) in binary{

            if let Some(c) = b.cubic { 

                c_binary.insert(key,c);

            }

            if let Some(a) = b.association {
                a_binary.insert(key,a);

            }
                
        }
        for record in pure{
            c_pure.push(record.c);
            a_pure.push(record.a);
        }

        let cubic = CubicParameters::from_raw(c_pure, c_binary, properties, opt.cubic_options)?;
        let assoc = AssociativeParameters::from_raw(a_pure, a_binary, None, ())?;
        let rdf_model: RDF = opt.rdf_model.into();


        let rdf = RDFcpa { 
            b: cubic.bij.diag().to_owned(),
            // bij: cubic.bij.clone(),
            model: rdf_model };
            
        Ok(CPAParameters{
            cubic,
            assoc,
            rdf
        })

    }
}

impl std::fmt::Display for CPAPureRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "CPAPureRecord(cubic={}, assoc={})", self.c, self.a)
        
    }
}

impl Default for CPAOptions {

    fn default() -> Self {
        
        let rdf_model = RDFmodelOption::default();
        let cubic_options = CubicOptions::classic_soave(CubicModelOption::SRK);
        Self { rdf_model, cubic_options }
    }
    
}

impl std::fmt::Display for CPABinaryRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let mut data = String::from_str("CPABinaryRecord(\n").unwrap();


        if let Some(c) = &self.cubic {
            data.push_str(concat!("cubic="));
            data.push_str(c.to_string().as_str());
            // let _ = write!(f, "CPABinaryRecord(cubic={:#?}, assoc={:#?})", c);
        } else {
            data.push_str("cubic=None");
        }
        data.push_str("\n");

        if let Some(a) = &self.association {
            data.push_str(concat!(" assoc="));
            data.push_str(a.to_string().as_str());            
        }
        
        else {
            data.push_str(" assoc=None");
        }

        data.push_str("\n)");

        write!(f, "{}", data)

    }
}

impl std::fmt::Display for CPAParameters {
    
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(f, "CPAParameters(\n  rdf_model='{}',\n  {},\n  {}\n)", self.rdf.model.to_string(),self.cubic, self.assoc,  )

    }

}