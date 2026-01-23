use std::{collections::HashMap};

use serde::{Deserialize, Serialize};
use crate::{parameters::records::BinaryParameter};
use super::parameters::{AssociationBinaryRecord};
use super::R;
use super::strength;

const W:[[f64;3];3]=[[0.0,1.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0]];

pub const A: usize = 0;
pub const B: usize = 1;
pub const C: usize = 2;
pub const NS:  usize = 3;
pub const SITES:[usize;3]=[A,B,C];

#[derive(Debug, Clone, PartialEq,Serialize)]
pub enum SiteType{
    A,
    B,
    C,
}


#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct Site{
    /// Type
    pub typ:SiteType,
    /// Owner
    pub owner:usize,
    /// Index
    pub idx:usize,
    /// Multiplicity
    pub mul:f64,

    pub epsilon:f64,
    pub kappa:f64,
}

impl Site {
    
    pub fn owner(&self)->usize{
        self.owner
    }

    pub fn typ_idx(&self)->usize{
        match self.typ{
            SiteType::A => A,
            SiteType::B => B,
            SiteType::C => C,
        }
    }

    pub fn epsilon(&self)->f64{
        self.epsilon
    }
    pub fn kappa(&self)->f64{
        self.kappa
    }

    pub fn has_owner(&self,i:usize)->bool{
        i == self.owner
    }
    
    pub fn is_self_associative(&self)->bool{

        if (self.epsilon() != 0.0) && (self.kappa() != 0.0) { true }
        else { false }
    
    }
    pub fn is_solvate(&self)->bool{

        if !self.is_self_associative(){ true }
        else { false }
    
    }

    pub fn interacts_with(&self,other:&Self)->bool{

        
        if (W[self.typ_idx()][other.typ_idx()] == 1.0) 
        && (self.is_self_associative()) 
        && (other.is_self_associative()) { true }
        else { false }
    }

    pub fn solvated_by(&self,other:&Self)->bool{

        if (W[self.typ_idx()][other.typ_idx()] == 1.0) 
        && (self.is_solvate() && other.is_self_associative())
        { true }
        else { false }
    }

    // pub fn interacts_with(&self,other:&Self)->bool{

        
    //     if (W[self.t()][other.t()] == 1.0) 
    //     && (self.is_associative()) 
    //     && (other.is_associative()) { true }
    //     else { false }
    // }

    pub fn new(typ:SiteType, owner:usize, idx:usize, mul:f64, epsilon:f64, kappa:f64)->Self{
        Self{
            typ,
            owner,
            idx,
            mul,
            epsilon,
            kappa
        }
    }
}


pub trait AssociationStrength: Default{

    fn dimensionless_delta_jl(t:f64,epsilon:f64,kappa:f64)->f64{
        (epsilon/R/t).exp_m1() * kappa
    }

    fn elliot_dimensionless_delta_jl(
        t:f64,
        epsilon_j:f64,
        kappa_j:f64,
        epsilon_l:f64,
        kappa_l:f64,
        )->f64{
        let dj = Self::dimensionless_delta_jl(t, epsilon_j, kappa_j);
        let dl = Self::dimensionless_delta_jl(t, epsilon_l, kappa_l);
        
        (dj*dl).sqrt()
    }

    fn cr1_factor_ik(f_ii:f64,f_kk:f64)->f64{
        0.5 * (f_ii + f_kk)
    }
    fn elliot_factor_ik(f_ii:f64,f_kk:f64)->f64{
        (f_ii * f_kk).sqrt()
    }

    
    fn association_strength_jl(
        &self,
        t:f64,
        f_ii:f64,
        f_kk:f64,
        interaction:&SiteInteraction)->f64;


}
#[derive(Serialize,Clone,Debug)]
pub struct SiteInteraction{
    pub site_j:usize, //Site -> usize 
    pub site_l:usize,
    pub epsilon:f64,
    pub kappa: f64,
    pub combining_rule:CombiningRule
}
#[derive(Serialize,Clone,Copy,PartialEq,Debug,Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CombiningRule{

    CR1,
    ECR
}

impl std::fmt::Display for CombiningRule {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self{
            CombiningRule::CR1 => write!(f,"cr1"),
            CombiningRule::ECR => write!(f,"ecr"),
        }
    }
}
impl Default for CombiningRule {
    fn default() -> Self {
        CombiningRule::CR1
    }
}

impl Into<CombiningRule> for Option<CombiningRule> {

    fn into(self) -> CombiningRule {

        self.unwrap_or_default()
    }
}

impl SiteInteraction {

    pub fn interactions_from_sites(sites:&Vec<Site>,binary:Vec<BinaryParameter<AssociationBinaryRecord>>)->Vec<Self>{

        let binary:HashMap<(usize,usize),AssociationBinaryRecord> = binary
        .into_iter()
        .map(|b| ((b.id1,b.id2),b.model_record)).collect();

        let mut interactions = Vec::<SiteInteraction>::new();
        let s = sites.len();
        
        for j in 0..s{
            for l in j..s{
                
                let site_j = &sites[j];
                let site_l = &sites[l];
                let i = site_j.owner;
                let k = site_l.owner;

                let opt = binary.get(&(i,k));
                
                // only cross-association
                if site_j.interacts_with(site_l){

                    interactions.push(
                        opt.map_or_else(
                        || Self::default(site_j.clone(), site_l.clone()), 
                        |b|{
                            Self::from_sites(site_j.clone(), site_l.clone(), b.epsilon, b.kappa, b.combining_rule)
                            }
                        )
                    )
                    
                // only induced-association
                } else if site_j.solvated_by(site_l) || site_l.solvated_by(site_j) {
                    
                    match opt {
                        
                        Some(b) => {
                            interactions.push(
                            Self::from_sites(site_j.clone(), site_l.clone(), b.epsilon, b.kappa, b.combining_rule)
                            )
                        }
                        None => continue
                    }
                }

            }
        }

        interactions

    }

    pub fn association_strength_jl(
        &self,
        t:f64,
        f_ii:f64,
        f_kk:f64,
        sites:&[Site])->f64{
        
        match &self.combining_rule {
            
            CombiningRule::CR1 => {

                strength::cr1_association_strength_jl(t, f_ii, f_kk, self.epsilon, self.kappa)
            }

            CombiningRule::ECR => {

                let j = self.site_j;
                let l = self.site_l;

                let sj = &sites[j];
                let sl = &sites[l];

                let epsilon_j = sj.epsilon;
                let epsilon_l = sl.epsilon;
                
                let kappa_j = sj.kappa;
                let kappa_l = sl.kappa;

                strength::ecr_association_strength_jl(t, f_ii, f_kk, epsilon_j, epsilon_l, kappa_j, kappa_l)
            }
        }

    }
}

impl SiteInteraction{

    fn arithmetic(epsilon_j:f64,epsilon_l:f64)->f64{
        0.5 * (epsilon_j + epsilon_l)
    }

    fn geometric(kappa_j:f64,kappa_l:f64)->f64{
        (kappa_j * kappa_l).sqrt()
    }

    pub fn from_sites(
        site_j:Site,
        site_l:Site,
        epsilon:Option<f64>,
        kappa:Option<f64>,
        combining_rule:CombiningRule)->Self{


        let epsilon = epsilon.unwrap_or(Self::arithmetic(site_j.epsilon,site_l.epsilon));
        let kappa = kappa.unwrap_or(Self::geometric(site_j.kappa,site_l.kappa));
        let combining_rule = combining_rule;

        Self { site_j: site_j.idx, site_l: site_l.idx, epsilon, kappa, combining_rule }
    }

    fn default(site_j:Site,site_l:Site)->Self{
        Self::from_sites(site_j, site_l, None, None, CombiningRule::default())
    }

    pub fn change_combining_rule(&mut self,combining_rule:CombiningRule){

        self.combining_rule = combining_rule;
        
    }
    pub fn change_cross_parameters(&mut self,epsilon:Option<f64>,kappa:Option<f64>){

        self.epsilon = epsilon.unwrap_or(self.epsilon);

        self.kappa = kappa.unwrap_or(self.kappa);

    } 
    pub fn belongs_to(&self,i:usize,k:usize,map:&[usize])->bool{


        let j = self.site_j;
        let l = self.site_l;

        let owner_j = map[j];
        let owner_l = map[l];

        let owners = (owner_j,owner_l);

        if owners == (i,k){
            true
        } else if owners == (k,i){
            true
        } else {
            false  
        }
    }
    
       
}

impl From<(Site,Site)> for SiteInteraction {
    fn from(value: (Site,Site)) -> Self {

        Self::from_sites(value.0, value.1,None,None,CombiningRule::default())
    }
}

impl std::fmt::Display for Site {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Site(type={},owner={},idx={},mul={},eps={},kappa={})",self.typ,self.owner,self.idx,self.mul,self.epsilon,self.kappa)
    }
}
impl std::fmt::Display for SiteInteraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SiteInteraction(j={},l={},epsilon={},kappa={}, rule='{}')",self.site_j,self.site_l,self.epsilon,self.kappa,self.combining_rule)
    }
}


impl std::fmt::Display for SiteType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self{
            SiteType::A => write!(f,"A"),
            SiteType::B => write!(f,"B"),
            SiteType::C => write!(f,"C"),
        }
    }
}




// impl From<Site> for (usize,usize) {

//     fn from(value: Site) -> Self {
        
//         (value.t(),value.c())
//     }
// }



#[cfg(test)]
pub mod tests{
    

    
    use serde_json::to_string_pretty;


    use super::*;



    #[test]
    fn test_binary_solvation(){
        
        let ew = 166.55e2;
        let bw = 0.0692;
        let wco2 = 0.1836;
        let site1 = Site::new(SiteType::A, 0, 0, 2.,ew, bw);
        let site2 = Site::new(SiteType::B, 0, 1, 2.,ew, bw);
        let site3 = Site::new(SiteType::B, 1, 2, 1.,0.0, 0.0);
        let sites = vec![site1,site2,site3];
        let b = AssociationBinaryRecord::new(None, Some(wco2), CombiningRule::default());
        let interactions = SiteInteraction::interactions_from_sites(&sites,vec![BinaryParameter::new(b, 0, 1)]);
        
        let i1 = &interactions[0];
        let i2 = &interactions[1];

        assert_eq!(i1.epsilon,ew);
        assert_eq!(i2.epsilon,0.5 * ew);
        assert_eq!(i2.kappa,wco2);

        for inter in &interactions{

            let json = to_string_pretty(inter).unwrap();
            println!("{}",json)
        }


    }
    
    #[test]
    fn test_ternary_with_solvation(){
        let ew = 166.55e2;
        let bw = 0.0692;
        let eacoh = 403.23e2;
        let bacoh = 4.5e-3;
        let wco2 = 0.1836;

        let site1 = Site::new(SiteType::A, 0, 0,2., ew, bw);
        let site2 = Site::new(SiteType::B, 0, 1,2., ew, bw);
        let site3 = Site::new(SiteType::C, 1, 2,1., eacoh, bacoh);
        let site4 = Site::new(SiteType::B, 2, 3,1., 0.0, 0.0);
        let b2 = AssociationBinaryRecord::new(None, Some(wco2), CombiningRule::default());
        let b1 = AssociationBinaryRecord::new(None, None, CombiningRule::ECR);
        
        let sites = vec![site1,site2,site3,site4];
        let b = vec![BinaryParameter::new(b1, 0, 1) , BinaryParameter::new(b2, 0, 2)];
        let interactions = SiteInteraction::interactions_from_sites(&sites,b);

        let i1 = &interactions[0];
        let i2 = &interactions[1];
        let i3 = &interactions[2];
        let i4 = &interactions[3];

        for inter in &interactions{

            let json = to_string_pretty(inter).unwrap();
            println!("{},",json);

        }
        let n = interactions.len();
        assert_eq!(n,5);

        assert_eq!(i1.epsilon,ew);
        assert_eq!(i1.kappa,bw);
        assert_eq!(i2.epsilon,0.5 * (ew + eacoh));
        assert_eq!(i2.kappa,(bw*bacoh).sqrt());
        assert_eq!(i4.epsilon,0.5 * (ew + eacoh));
        assert_eq!(i4.kappa,(bw * bacoh).sqrt());
        assert_eq!(i3.epsilon,0.5 * ew);
        assert_eq!(i3.kappa,wco2);
        assert_eq!(i2.combining_rule,CombiningRule::ECR);
        assert_eq!(i4.combining_rule,CombiningRule::ECR);


    }

}