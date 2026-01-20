use std::{sync::Arc};

use ndarray::{Array1};

use crate::{residual::Residual, state::E};
pub mod tpd;
pub mod vle;

//Contém EOS
// Antoine Parameters
pub struct PhaseEquilibrium<R>{
    eos:Arc<E<R>>,
    correl:Option<Antoine>
}


impl<R:Residual> PhaseEquilibrium<R> {
    
    pub fn new(eos:&Arc<E<R>>,correl:Option<Antoine>)->Self{
        Self { eos:eos.clone(), correl }
    }
}
// se o record foi montado em outra ordem (ou apenas uma subs da mistura tenha parametro
// de antoine, ent )

// melhor solução é ser retornado um map de pressoes (i,f64)
// de modo que, caso i nao esteja no map, será retornado pressao 1bar



// pub struct Antoine(HashMap<usize,AntoineRecord>);
pub struct Antoine(Array1<AntoineRecord>);

impl Antoine {

    fn default(shape:usize)->Self{
        Antoine(Array1::default(shape))
    }

    fn from_records(records:Vec<AntoineRecord>)->Self{
        
        // let mut map:HashMap<usize,AntoineRecord>=HashMap::new();
        let mut arr:Array1<AntoineRecord>=Array1::default(records.len());
        for (i,r) in records.iter().enumerate(){
            arr[i]=r.clone();
            // map.insert(i, r.clone());
        }

        Self(arr)
    }


    fn psat(&self,t:f64)->Array1<f64>{

        let mut pres=Array1::zeros(self.0.len());
        for (j,i) in self.0.iter().enumerate(){
            let (a,b,c)=(i.a,i.b,i.c);
            let log= a-b/(t+c);
            match i.base{
                LogBase::Log10=>{ pres[j]=10.0_f64.powf(log)}
                LogBase::LogE=> { pres[j]=log.exp()         }
            }
        }
        //bar -> pa
        pres*1e5
    }
}


#[derive(Clone)]
pub struct AntoineRecord{
    a:f64,
    b:f64,
    c:f64,
    base:LogBase
}

#[derive(Clone)]
pub enum LogBase{
    Log10,
    LogE
}
impl AntoineRecord {
    
    fn record(a:f64,b:f64,c:f64,log:LogBase)->Self{
        Self{a,b,c,base:log}
    }
}

impl Default for AntoineRecord{
    
    fn default() -> Self {
        Self { a: 10.0, b: 0.0, c: 0.0, base: LogBase::Log10 }
    }
}


