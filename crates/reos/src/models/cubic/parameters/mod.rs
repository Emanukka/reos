use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize, ser::SerializeStruct};
use crate::{models::cubic::{CubicModel}, parameters::{Parameters, PureRecord}};


#[derive(Debug)]
pub struct CubicParameters<T:CubicModel>{
    pub ncomp:usize,
    pub a0:   Array2<f64>,
    pub b:    Array2<f64>,
    pub tc:   Array2<f64>,
    pub kappa:Array2<f64>,
    pub aij:  Array2<f64>, //kij = Aij + Bij*T
    pub bij:  Array2<f64>,
    // pub w: Option<Array<f64>>
    pub model: T

}
#[derive(Clone,Debug,Serialize, Deserialize)]

pub struct CubicPureRecord{
    pub a0:f64,
    pub b:f64,
    pub c1:f64,
    pub tc:f64,
}

impl CubicPureRecord {
    
    pub fn new(a0:f64,b:f64,c1:f64,tc:f64)->Self{
        Self{
            a0,b,c1,tc  
        }
    }
}

impl PureRecord for CubicPureRecord {
    
}

impl<T:CubicModel> CubicParameters<T>{
    pub fn new(model:T,records:Vec<CubicPureRecord>)->Self {
        let ncomp =records.len();

        let mut a0:Array2<f64>    = Array2::<f64>::zeros((ncomp,1));
        let mut b:Array2<f64>     = Array2::<f64>::zeros((ncomp,1));
        let mut tc:Array2<f64>    = Array2::<f64>::zeros((ncomp,1));
        let mut kappa:Array2<f64> = Array2::<f64>::zeros((ncomp,1));
        let aij:Array2<f64>       = Array2::<f64>::zeros((ncomp,ncomp));
        let bij:Array2<f64>       = Array2::<f64>::zeros((ncomp,ncomp));


        for (i,record) in records.iter().enumerate(){
            a0[(i,0)]=record.a0;
            b[(i,0)] =record.b;
            kappa[(i,0)]=record.c1;
            tc[(i,0)]=record.tc;
        }    
        Self{
            ncomp,
            a0,
            b:b.clone(),
            kappa,
            tc,
            aij,
            bij,
            model
        } 
    }

    pub fn set_kij(&mut self,i:usize,j:usize,kij:f64){

        self.aij[(i,j)] = kij;
        self.aij[(j,i)] = kij;
    }
    pub fn set_kij_temperature_dependent(&mut self,i:usize,j:usize,aij:f64,bij:f64){

        self.aij[(i,j)] = aij;
        self.aij[(j,i)] = aij; 
 
        self.bij[(i,j)] = bij;
        self.bij[(j,i)] = bij;
    }

    pub fn from_tpw(
        model:T,
        tc:Array1<f64>,
        pc:Array1<f64>,
        w:Array1<f64>,)->Self{
            let ncomp =tc.len();
            assert!(pc.len() == ncomp,"pc length mismatch");
            assert!(w.len()  == ncomp,"w length mismatch");

            let a0:Array2<f64>    = model.acrit(&tc, &pc).to_shape((ncomp,1)).unwrap().to_owned();
            let b:Array2<f64>     = model.bcrit(&tc, &pc).to_shape((ncomp,1)).unwrap().to_owned();
            let tc:Array2<f64>    = tc.to_shape((ncomp,1)).unwrap().to_owned();
            let aij:Array2<f64>   = Array2::<f64>::zeros((ncomp,ncomp));
            let bij:Array2<f64>   = Array2::<f64>::zeros((ncomp,ncomp));
            let kappa:Array2<f64> = model.kappa_from_w(&w).to_shape((ncomp,1)).unwrap().to_owned();

            Self{
                ncomp,
                a0,
                b:b.clone(),
                kappa,
                tc,
                aij,
                bij,
                model
            }   
        }
}

impl<C:CubicModel> Parameters<CubicPureRecord> for CubicParameters<C> {
    fn from_records(records:Vec<CubicPureRecord>)->Self {
        
        Self::new(C::model(), records)
    }
}

impl<T:CubicModel> Serialize for CubicParameters<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: serde::Serializer {
        
        let mut state = serializer.serialize_struct("CubicParameters", 5)?;
        state.serialize_field("cubic_model", &self.model.which())?;

        state.serialize_field("a0", &self.a0
        .flatten()
        .iter()
        .map(|v| format!("{:.6e}", v))
        .collect::<Vec<_>>())?;
        
        state.serialize_field("b", &self.b
        .flatten()
        .iter()
        .map(|v| format!("{:.6e}", v))
        .collect::<Vec<_>>())?;

        state.serialize_field("kappa", &self.kappa
        .flatten()
        .iter()
        .map(|v| format!("{:.6e}", v))
        .collect::<Vec<_>>())?;
        
        state.serialize_field("tc", &self.tc
        .flatten()
        .iter()
        .map(|v| format!("{:.6e}", v))
        .collect::<Vec<_>>())?;

        state.end()        
    }
}


#[cfg(test)]
mod tests{
    use ndarray::array;

    use crate::models::cubic::{PR76, PR78};


    #[test]
    fn test1(){
        let tc = array![ 507.317, 860.372];
        let pc = array![ 3290000., 1246000.];
        let w = array![ 0.275 , 1.022  ];
        let params = super::CubicParameters::from_tpw(PR78, tc, pc, w);

        let json = serde_json::to_string_pretty(&params).unwrap();

        
        println!("{}", json);

    }

    #[test]
    fn test2(){
        let tc = array![ 507.317, 860.372];
        let pc = array![ 3290000., 1246000.];
        let w = array![ 0.275 , 1.022  ];
        let params = super::CubicParameters::from_tpw(PR76, tc, pc, w);

        let json = serde_json::to_string_pretty(&params).unwrap();
        // serde_json::to_writer_pretty(&mut std::fs::File::create("pr76_params.json").unwrap(), &params).unwrap();
        println!("{}", json);

    }

}