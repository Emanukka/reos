use serde::{Deserialize, Serialize};

use super::sites::{CombiningRule, NS, Site, SiteInteraction, SiteType};
use crate::parameters::{Parameters, Properties, records::BinaryParameter};



#[derive(Clone)]
pub struct AssociativeParameters{
    pub sites:        Vec<Site>,
    pub interactions: Vec<SiteInteraction>,
    pub na:           usize,
    pub nb:           usize,
    pub nc:           usize,
}

type Pure = AssociationPureRecord;
type Binary = AssociationBinaryRecord;

impl Parameters<Pure, Binary, ()> for AssociativeParameters {

    fn from_raw(pure:Vec<Pure>, binary: Vec<BinaryParameter<Binary>>, _: Option<Properties>, _: ()) -> Self {
        
        let n = pure.len();

        let mut n_a = 0;
        let mut n_b = 0; 
        let mut n_c = 0; 

        let mut s = 0;
        let mut sites:Vec<Site> = vec![]; 

        for i in 0..n{

            let record = &pure[i];
            
            let na = record.na as f64;
            let nb = record.nb as f64;
            let nc = record.nc as f64;
            

            if na != 0.0 {
                sites.push(Site::new(SiteType::A, i, s,na, record.epsilon,record.kappa));
                n_a += 1;
                s += 1;
            }
            if nb != 0.0 {
                sites.push(Site::new(SiteType::B, i, s,nb, record.epsilon,record.kappa));
                n_b += 1;
                s += 1;
            }
            if nc != 0.0{
                sites.push(Site::new(SiteType::C, i, s,nc,record.epsilon,record.kappa));
                n_c += 1;
                s += 1;
            }

        }

        let interactions = SiteInteraction::interactions_from_sites(&sites, binary); 

        AssociativeParameters{
            na: n_a,
            nb: n_b,
            nc: n_c,
            interactions,
            sites,
        }

    }
}

impl std::fmt::Display for AssociativeParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        let ssites = self.sites.iter()
        .map(|s| s.to_string())
        .collect::<Vec<String>>()
        .join(",\n\t   ");
        
        let sinter = self.interactions.iter()
        .map(|s| s.to_string())
        .collect::<Vec<String>>()
        .join(",\n\t   ");

        write!(f, 
            "AssociativeParameters(\n\tna={}, nb={}, nc={},\n\tsites=[\n\t   {}],\n\tinteractions=[\n\t   {}])",
            self.na,
            self.nb,
            self.nc,
            ssites,
            sinter,
        )


    }
}

pub enum AssociationType{
    Inert,
    // Solvate,
    Associative
}

impl Default for AssociationType {
    
    fn default() -> Self {
        Self::Inert
    }
}

#[derive(Serialize, Deserialize,Debug,Clone,PartialEq)]
pub struct AssociationPureRecord{

    #[serde(default)]
    pub epsilon:f64,
    #[serde(default)]
    pub kappa:f64,
    #[serde(default)]
    pub na:usize,
    #[serde(default)]
    pub nb:usize,
    #[serde(default)]
    pub nc:usize,
    // #[serde(default)]
    // pub typ:AssociationType,

}
impl std::fmt::Display for AssociationPureRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "AssociationPureRecord(epsilon={}, kappa={}, na={}, nb={}, nc={})", self.epsilon, self.kappa, self.na, self.nb, self.nc)
        
    }
}
#[derive(Clone,Serialize,Deserialize,Debug)]
pub struct AssociationBinaryRecord{
    
    #[serde(default)]
    pub epsilon:Option<f64>,
    #[serde(default)]
    pub kappa:Option<f64>,
    #[serde(default= "CombiningRule::default")]
    #[serde(rename = "rule")]
    pub combining_rule:CombiningRule
}


impl AssociationBinaryRecord {
    
    pub fn new(epsilon:Option<f64>,kappa:Option<f64>,combining_rule:CombiningRule)->Self{

        Self { epsilon, kappa, combining_rule }
    }

}


impl std::fmt::Display for AssociationBinaryRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "AssociationBinaryRecord(epsilon={:?}, kappa={:?}, combining_rule={:?})", 
        self.epsilon, 
        self.kappa, 
        self.combining_rule)
        
    }
}

impl AssociationPureRecord {
    
    pub fn new(
        epsilon:f64,
        kappa:f64,
        na:usize,
        nb:usize,
        nc:usize,
    )->Self{
        Self{
        epsilon,
        kappa,
        na,
        nb,
        nc,
        // typ,
        }
    }
    pub fn inert() -> Self {
        AssociationPureRecord{
        epsilon:0.0,
        kappa:0.0,
        na:0,
        nb:0,
        nc:0,
        }
    }

    pub fn solvate(sites:[usize;NS])-> Self{
        AssociationPureRecord{
        epsilon:0.0,
        kappa:0.0,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        }
    }

    pub fn associative(epsilon:f64,kappa:f64,sites:[usize;NS])->Self{
        AssociationPureRecord{
        epsilon,
        kappa,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        }
    }
    
    pub fn get_type(&self)->AssociationType{
        
        if (self.na == 0) && (self.nb == 0) && (self.nc == 0){
            AssociationType::Inert
        } else {
           AssociationType::Associative
        }

    }

}





#[cfg(test)]
mod tests{

    use serde_json::from_str;
    use crate::parameters::records::{BinaryRecord, PureRecord};

    use super::*;


    #[test]
    fn test_induced_association_json(){

        let data1 = r#"
        {   
            "name": "water",
            "epsilon": 166.55e2, 
            "kappa": 0.0692,
            "na":   2,
            "nb":   2
        }
        "#;
        
        let data2 = r#"
        {   
            "name": "co2",
            "nb":   1
        }
        "#;

        let data3 = r#"
        {
            "kappa": 0.1836,
            "id1": "water",
            "id2": "co2"

        }
        "#;
        let pr1: PureRecord<AssociationPureRecord> = from_str(data1).unwrap();
        let pr2: PureRecord<AssociationPureRecord> = from_str(data2).unwrap();
        let br: BinaryRecord<AssociationBinaryRecord> = from_str(data3).unwrap();

        // let c1_string = serde_json::to_string_pretty(&c1).unwrap();

        // let pr1 = PureRecord::new(0.0, "water".into(), c1);
        // let pr2 = PureRecord::new(0.0, "co2".into(), c2);
        // let induced = AssociationBinaryRecord::new(None,Some(0.1836), None);
        // let br = BinaryRecord::new(induced, "water".into(), "co2".into());

        let p = AssociativeParameters::new(vec![pr1,pr2], vec![br], ());

        let string = p.to_string();

        println!("{}",string);
        // let string = p.to_string();
        let inter = p.interactions;
        let na = p.na;
        let nb = p.nb;

        let self_water = &inter[0];
        let water_co2 = &inter[1];

        assert_eq!(na * nb, 2);

        assert_eq!(self_water.epsilon, 166.55e2);
        assert_eq!(self_water.kappa, 0.0692);
        assert_eq!(water_co2.epsilon, 0.5 * 166.55e2);
        assert_eq!(water_co2.kappa, 0.1836);

        // println!("{}",string)

    }
}