use serde::{Deserialize, Serialize};

use super::sites::{NS, Site, SiteType};
use super::sites_interaction::SitesInteraction;
use super::combining_rule::*;

use crate::{ parameters::{BinaryMap, Parameters, Properties}};



#[derive(Clone)]
pub struct AssociativeParameters{
    pub sites:        Vec<Site>,
    pub interactions: Vec<SitesInteraction>,
    pub na:           usize,
    pub nb:           usize,
    pub nc:           usize,
}

impl Parameters for AssociativeParameters {
    type Pure = AssociationPureRecord;
    type Binary = AssociationBinaryRecord;
    type Options = ();

    fn from_raw(pure:Vec<Self::Pure>, binary: BinaryMap<Self::Binary>, _: Option<Properties>, _: ()) -> Result<Self, Box<dyn std::error::Error>> {
        
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

        let interactions = SitesInteraction::new_interactions(&sites, binary); 

        Ok( 
            AssociativeParameters{na: n_a, nb: n_b, nc: n_c, interactions, sites}
        )

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

}

#[derive(Clone,Serialize,Deserialize,Debug)]
pub enum AssociationBinaryRecord{
    
    Set {epsilon:f64, kappa:f64},
    CombiningRule (CombiningRuleOption)    

}



impl AssociationPureRecord {
    
    pub fn new(epsilon:f64, kappa:f64, na:usize, nb:usize, nc:usize) -> Self{
        Self{ epsilon, kappa, na, nb, nc }
    }

    pub fn inert() -> Self {
        Self::new(0., 0., 0, 0, 0)
    }

    pub fn solvate(sites:[usize;NS])-> Self{
        Self::new(0., 0.,sites[0],sites[1],sites[2],)
    }

    pub fn associative(epsilon:f64,kappa:f64,sites:[usize;NS])->Self{
        Self::new(epsilon, kappa,sites[0],sites[1],sites[2],)
    }

}



impl Default for AssociationBinaryRecord {

    fn default() -> Self {
        Self::CombiningRule(CombiningRuleOption::CR1)
    }
}

impl std::fmt::Display for AssociationPureRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        write!(f, "AssociationPureRecord(epsilon={}, kappa={}, na={}, nb={}, nc={})", self.epsilon, self.kappa, self.na, self.nb, self.nc)
        
    }
}

impl std::fmt::Display for AssociationBinaryRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        match &self{

            Self::Set { epsilon, kappa } => {
                write!(f, "AssociationBinaryRecord(epsilon={}, kappa={})", epsilon, kappa)
            }

            Self::CombiningRule(combining_rule) => {
                write!(f, "AssociationBinaryRecord({})", combining_rule)
                
            }
        }
    }
}
