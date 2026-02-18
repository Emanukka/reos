use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use crate::{models::cubic::{alpha::{Alpha, AlphaRecord}, combining_rule::CombiningRule, mixing_rule::{MixingRule, MixingRuleModel}, models::{ CubicModel, CubicModels, PR76, PR78}, options::CubicOptions},
parameters::{ BinaryMap, Parameters, records::{BinaryParameter, Properties}}};


#[derive(Clone,Debug,Serialize, Deserialize)]
pub struct CubicBinaryRecord{
    pub kij:f64,
    #[serde(default)]
    pub lij:f64
}

impl std::fmt::Display for CubicBinaryRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let mut ret = vec!["("];        
        let mut kl = vec![];

        if self.kij != 0. {
            kl.push(format!("kij={}", self.kij));
        }
        
        if self.lij != 0. {
            kl.push(format!("lij={}", self.lij));
        }
        
        let mut s = kl.join(", ");
        ret.push(&mut s);
        ret.push(")");
        

        write!(f, "{}", ret.join(""))
    }
}
impl Default for CubicBinaryRecord {
    fn default() -> Self {
        Self { kij: 0., lij: 0. }
    }
} 

#[derive(Clone,Debug,Serialize, Deserialize)]
pub struct CubicPureRecord{

    pub tc:f64,
    pub c:Option<f64>, //Volume translation 
    #[serde(flatten)]
    pub parameters: PureParameters,
    #[serde(flatten)]
    pub alpha:AlphaRecord,

}

#[derive(Clone,Debug,Serialize, Deserialize)]
#[serde(untagged)]
pub enum PureParameters{

    Classic{
        pc:f64,
    },

    Regressed{
        a:f64,
        b:f64,
    }
}

impl std::fmt::Display for PureParameters {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        match self{
            Self::Classic { pc } => write!(f, "{pc} Pa",),
            Self::Regressed { a, b } => write!(f, "a={a}, b={b}")
        }
        
    }
}
impl CubicPureRecord {

    pub fn classic(tc:f64, pc:f64, alpha:AlphaRecord, c:Option<f64>,) -> Self {
        let parameters = PureParameters::Classic { pc };
        Self { tc, parameters, c, alpha }
    }
    pub fn regressed(a:f64, b:f64, tc:f64, alpha:AlphaRecord, c:Option<f64>,) -> Self {
        let parameters = PureParameters::Regressed { a, b } ;
        Self { tc, parameters, c, alpha }
    }

    pub fn classic_soave(tc:f64, pc:f64, w:f64, c:Option<f64>) -> Self {

        Self::classic(tc, pc, AlphaRecord::Soave { w }, c )
    }

    pub fn regressed_soave(a:f64, b:f64, tc:f64, c1:f64, c:Option<f64>) -> Self {
        Self::regressed(a, b, tc, AlphaRecord::SoaveRegressed { c1 }, c)
    }
    
}


impl std::fmt::Display for CubicPureRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "CubicPureRecord(tc={} K, {}, {})", self.tc, self.parameters, self.alpha)

    }
}


#[derive(Clone,Debug)]
pub struct CubicParameters {
    pub a: Vec<f64>,
    pub b: Vec<f64>,
    pub c: Vec<f64>,
    pub tc:Vec<f64>,
    pub options:CubicOptions,
    pub binary:Array2<CubicBinaryRecord>,
    pub properties: Properties,
}



impl Parameters for CubicParameters {
    type Pure = CubicPureRecord;
    type Binary = CubicBinaryRecord;
    type Options = CubicOptions;

    fn from_raw(
        pure:Vec<Self::Pure>, 
        binary: crate::parameters::BinaryMap<Self::Binary>, 
        properties: Option<Properties>, 
        opt: Self::Options) -> Result<Self, Box<dyn std::error::Error>> {

        let n = pure.len();
        
        let model = opt.model;
        let alpha = opt.alpha;
        let combr = opt.combr;
        let mix = opt.mix;
        
        let properties = properties.unwrap_or_default();

        let mut a_ = Vec::with_capacity(n);
        let mut b_ = Vec::with_capacity(n);
        let mut c = Vec::with_capacity(n);
        let mut tc= Vec::with_capacity(n);
        let mut alpha_records = Vec::with_capacity(n);
        let mut binary_ = Array2::default((n,n));

        for r in pure.into_iter(){
            
            match r.parameters {

                PureParameters::Classic { pc } => {
                    a_.push(model.acrit(r.tc, pc));
                    b_.push(model.bcrit(r.tc, pc));
                }
                PureParameters::Regressed { a, b } => {
                    a_.push(a);
                    b_.push(b);
                }
            }

            c.push(r.c.unwrap_or_default());
            tc.push(r.tc);
            alpha_records.push(r.alpha);

        }

        binary.into_iter().for_each(|((i,j),b)| {
            binary_[(i,j)] = b.clone();
            binary_[(j,i)] = b;
        });

        let alpha = alpha.build(&properties.names, alpha_records, &model)?;

        let options = CubicOptions::new(model, alpha, combr, mix);
        Ok(
        Self{
            a:a_,
            b:b_,
            c,
            tc,
            options,
            binary:binary_,
            properties,
        })


    }
}


impl std::fmt::Display for CubicParameters {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        if self.a.len() == 1 {
            write!(f, "CubicParameters(\n\tmodel={:?},\n\tac={:?},\n\tbc={:?},\n\ttc={:?},\n\t{})",
                self.options.model.to_string(),
                self.a,
                self.b,
                self.tc,
                self.options.alpha.to_string(),
            )

        } else {

            let sbin = self.binary
            .rows()
            .into_iter()
            .map(|row| row.to_string())
            .collect::<Vec<_>>()
            .join("\n\t       ");

            write!(f, "CubicParameters(\n\tmodel={:?},\n\tmixing={:?},\n\tac={:?},\n\tbc={:?},\n\ttc={:?},\n\t{},\n\tbinary=[{}]\n)",
                self.options.model.to_string(),
                self.options.mix.to_string(),
                self.a,
                self.b,
                self.tc,
                self.options.alpha.to_string(),
                sbin
            )

        }

    }
}

