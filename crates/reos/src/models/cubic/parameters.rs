use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};

use super::alpha::{Alpha, AlphaRecord};
use super::models::{CubicModels,CubicModel};
use super::mixing_rule::{MixingRule,MixingRuleModel};
use super::combining_rule::{CombiningRule,CombiningRuleModel};

use super::options::*;



#[derive(Clone, Serialize, Deserialize, Debug, PartialEq)]
#[serde(rename_all = "lowercase")]
#[serde(deny_unknown_fields)]
/// kij(T) = a + b * T
pub struct Kij{ 
    #[serde(default)]
    pub a:f64,
    #[serde(default)]
    pub b:f64   
}


#[derive(Clone,Debug,Serialize, Deserialize)]
pub struct CubicBinaryRecord{

    pub kij: Kij,
    // pub kij:f64,
    // #[serde(default)]
    // pub lij:f64
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
        #[serde(alias="a0")]
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
    pub aij: Array2<f64>,
    pub bij: Array2<f64>,
    pub c: Vec<f64>,
    pub tc:Vec<f64>,
    pub kij:Array2<Kij>,
    pub model:CubicModels,
    pub combr:CombiningRule,
    pub mix: MixingRule,
    pub alpha:Alpha,
    pub properties: crate::parameters::Properties,

}



impl crate::parameters::Parameters for CubicParameters {
    type Pure = CubicPureRecord;
    type Binary = CubicBinaryRecord;
    type Options = CubicOptions;

    fn from_raw(
        pure:Vec<Self::Pure>, 
        binary: crate::parameters::BinaryMap<Self::Binary>, 
        properties: Option<crate::parameters::Properties>, 
        opt: Self::Options) -> Result<Self, Box<dyn std::error::Error>> {

        let n = pure.len();
        
        let model: CubicModels = opt.cubic_model.into();
        let alpha: Alpha = opt.alpha_model.into();
        let combr: CombiningRule = opt.combining_rule.into();
        let mix = opt.mixing_rule.into();
        
        let properties = properties.unwrap_or_default();

        let [mut a_, mut b_] = [Vec::with_capacity(n), Vec::with_capacity(n)];

        let mut c = Vec::with_capacity(n);
        let mut tc= Vec::with_capacity(n);
        let mut alpha_records = Vec::with_capacity(n);
        // let mut binary_ = Array2::default((n,n));
        let mut kij = Array2::default((n,n));

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

            kij[(i,j)] = b.kij.clone();
            kij[(j,i)] = b.kij;

        });

        let alpha = alpha.build(&properties.names, alpha_records, &model)?;
        let [aij, bij] = combr.apply(a_, b_);
        
        // let options = CubicOptions::new(model, alpha, combr, mix);
        Ok(
        Self{
            aij,
            bij,
            // a:a_,
            // b:b_,

            c,
            tc,
            model,
            mix,
            combr,
            alpha,
            // binary:binary_,
            kij,
            properties,
        })


    }
}


impl std::fmt::Display for Kij {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(f, "({},{})", self.a, self.b)
        
    }
}
impl std::fmt::Display for CubicBinaryRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "kij{}", self.kij)

    }
}

impl Default for Kij {
    fn default() -> Self {
        Self { a: 0., b: 0. }
    }
} 
impl Default for CubicBinaryRecord {
    fn default() -> Self {
        Self { kij: Kij::default() }
    }
} 

impl std::fmt::Display for CubicParameters {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        
        if self.tc.len() == 1 {
            write!(f, "CubicParameters(\n\tmodel={:?},\n\ta={:?},\n\tb={:?},\n\ttc={:?},\n\t{})",
                self.model.to_string(),
                self.aij[(0,0)],
                self.bij[(0,0)],
                self.tc[0],
                self.alpha.to_string(),
            )

        } else {

            let mut ab = vec![];

            for mat in [&self.aij, &self.bij] {
                let s = mat
                .rows()
                .into_iter()
                .map(|row| row.to_string())
                .collect::<Vec<_>>()
                .join(",\n\t    ");
                ab.push(s);
            }

            let mut kij = vec![];

            for mat in [&self.kij] {
                let s = mat
                .rows()
                .into_iter()
                .map(|row| row.to_string())
                .collect::<Vec<_>>()
                .join(",\n\t    ");
                kij.push(s);

            }
            // let sbin = self.kij
            // .rows()
            // .into_iter()
            // .map(|row| row.to_string())
            // .collect::<Vec<_>>()
            // .join("\n\t       ");
            // let saij = self.aij
            // .rows()
            // .into_iter()
            // .map(|row| row.to_string())
            // .collect::<Vec<_>>()
            // .join("\n\t       ");
            // let sbij = self.bij
            // .rows()
            // .into_iter()
            // .map(|row| row.to_string())
            // .collect::<Vec<_>>()
            // .join("\n\t       ");

            write!(f, "CubicParameters(\n\tmodel={:?}, mixing_rule={:?}, combining_rule={:?},\n\ttc={:?},\n\t{},\n\taij=[{}],\n\tbij=[{}],\n\tkij=[{}]\n)",
                self.model.to_string(),
                self.mix.to_string(),
                self.combr.to_string(),
                self.tc,
                self.alpha.to_string(),
                ab[0],
                ab[1],
                kij[0]
            )

        }

    }
}


#[cfg(test)]
mod tests{
    
    use serde::{Serialize, Deserialize};
    #[derive(Serialize, Deserialize, Debug, PartialEq)]
    #[serde(rename = "kij")]
    struct Kij{
        a:f64,
        #[serde(default)]
        b:f64   
    }


    #[test]
    fn test1(){

        let kij: Kij = Kij { a: 1., b: 2. };

        let s = r#"
        {
            "a":1.0,
            "b":2.0
        }

        "#;
        
        assert_eq!(serde_json::from_str::<Kij>(s).unwrap(), kij);
    }   
    #[test]
    fn test2(){

        let kij: Kij = Kij { a: 1., b: 0.};

        let s = r#"
        {
            "a":1.0
        }

        "#;
        
        assert_eq!(serde_json::from_str::<Kij>(s).unwrap(), kij);
    }   

}