use std::{collections::HashMap, marker::PhantomData, path::Iter};

use serde::Serialize;

use crate::{models::{IDEAL_GAS_CONST as R, associative::parameters::{AssociationBinaryRecord, W}}, state::E};


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

#[derive(Debug, Clone, PartialEq,Serialize)]
pub struct Site{
    pub t:SiteType,
    ///Owner
    pub c:usize,
    ///Index
    pub idx:usize,

    pub epsilon:f64,
    pub kappa:f64,
}

impl Site {
    
    ///Get component of site alpha
    pub fn c(&self)->usize{
        self.c
    }
    ///Get type index site alpha
    pub fn t(&self)->usize{
        match self.t{
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
        i==self.c()
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

        
        if (W[self.t()][other.t()] == 1.0) 
        && (self.is_self_associative()) 
        && (other.is_self_associative()) { true }
        else { false }
    }

    pub fn solvated_by(&self,other:&Self)->bool{

        if (W[self.t()][other.t()] == 1.0) 
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

    pub fn new(t:SiteType,c:usize,j:usize,epsilon:f64,kappa:f64)->Self{
        Self{
            t,
            c,
            idx:j,
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
#[derive(Serialize,Clone)]
pub struct SiteInteraction{
    pub site_j:Site,
    pub site_l:Site,
    epsilon:f64,
    kappa: f64,
    pub combining_rule:CombiningRule
}
#[derive(Serialize,Clone,Copy)]
pub enum CombiningRule{
    CR1,
    ECR
}

impl Default for CombiningRule {
    fn default() -> Self {
        CombiningRule::CR1
    }
}

impl AssociationStrength for CombiningRule {
    

    fn association_strength_jl(
        &self,
        t:f64,
        f_ii:f64,
        f_kk:f64,
        interaction:&SiteInteraction)->f64 {
        match self {

            CombiningRule::CR1 => {

                let epsilon = interaction.epsilon;
                let kappa= interaction.kappa;
                // (epsilon/R/t).exp_m1() * kappa
                let d_jl = Self::dimensionless_delta_jl(t, epsilon, kappa);
                let f_ik = Self::cr1_factor_ik(f_ii,f_kk) ;
                
                f_ik * d_jl

            }

            CombiningRule::ECR => {
                let epsilon_j = interaction.site_j.epsilon();
                let kappa_j = interaction.site_j.kappa();
                let epsilon_l = interaction.site_l.epsilon();
                let kappa_l = interaction.site_l.kappa();

                let d_jl = Self::elliot_dimensionless_delta_jl(t, epsilon_j, kappa_j, epsilon_l, kappa_l);
                let f_ik= Self::elliot_factor_ik(f_ii, f_kk);

                f_ik * d_jl

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
        combining_rule:Option<CombiningRule>)->Self{


        let epsilon = epsilon.unwrap_or(Self::arithmetic(site_j.epsilon,site_l.epsilon));
        let kappa = kappa.unwrap_or(Self::geometric(site_j.kappa,site_l.kappa));
        let combining_rule = combining_rule.unwrap_or_default();

        Self { site_j, site_l, epsilon, kappa, combining_rule }
    }


    pub fn interactions_from_sites(sites:&Vec<Site>,binary:HashMap<(usize,usize),AssociationBinaryRecord>)->Vec<Self>{

        let mut interactions = Vec::<SiteInteraction>::new();
        let s = sites.len();
        for j in 0..s{
            for l in j..s{
                
                let site_j = &sites[j];
                let site_l = &sites[l];
                
                if site_j.interacts_with(site_l){
                    
                    let comp1 = site_j.c();
                    let comp2 = site_l.c();
                    let k1 = (comp1,comp2);
                    let k2 = (comp2,comp1);
                    let site_j = site_j.clone();
                    let site_l = site_l.clone();

                    if binary.contains_key(&k1){
                        let bin = binary.get(&k1).unwrap();
                        let interaction = Self::from_sites(site_j, site_l,bin.epsilon,bin.kappa,bin.combining_rule);
                        interactions.push(interaction)

                    } else if binary.contains_key(&k2)  {
                        let bin = binary.get(&k2).unwrap();
                        let interaction = Self::from_sites(site_j, site_l,bin.epsilon,bin.kappa,bin.combining_rule);
                        interactions.push(interaction)
                    }else {
                        let interaction = Self::from_sites(site_j, site_l,None,None,None);
                        interactions.push(interaction)
                    }
                }
            }
        }
        interactions
    }

    pub fn solvations_from_ik(
        i:usize,
        k:usize,
        sites:&Vec<Site>,
        epsilon:Option<f64>,
        kappa:f64)->Vec<Self>{
        
        let s = sites.len();
        let mut interactions = Vec::new();
        for j in 0..s{
            for l in j..s{
                
                let site_j = sites[j].clone();
                let site_l = sites[l].clone();
                let mut interaction = SiteInteraction::from((site_j,site_l));

                if interaction.belongs_to(i, k) && interaction.is_solvation(){

                    interaction.change_cross_parameters(epsilon, Some(kappa));
                    interactions.push(interaction);

                } else {
                  continue;  
                }

            }
        }
        interactions


    }

    pub fn change_combining_rule(&mut self,combining_rule:CombiningRule){

        self.combining_rule = combining_rule;
        
    }
    pub fn change_cross_parameters(&mut self,epsilon:Option<f64>,kappa:Option<f64>){

        self.epsilon = epsilon.unwrap_or(self.epsilon);

        self.kappa = kappa.unwrap_or(self.kappa);

    } 


    pub fn association_strength_jl(
        &self,
        t:f64,
        f_ii:f64,
        f_kk:f64,)->f64{
        
        self.combining_rule.association_strength_jl(t, f_ii, f_kk, self)
        
    }

    pub fn belongs_to(&self,i:usize,k:usize)->bool{

        let owner_j = self.site_j.c();
        let owner_l = self.site_l.c();
        let owners = (owner_j,owner_l);

        if owners == (i,k){
            true
        } else if owners == (k,i){
            true
        } else {
            false  
        }
    }
    
    pub fn is_solvation(&self)->bool{

        let sj = &self.site_j;
        let sl = &self.site_l;

        if sj.solvated_by(sl) || sl.solvated_by(sj){
            true
        } else {
          false  
        }

    }
       
}

impl From<(Site,Site)> for SiteInteraction {
    fn from(value: (Site,Site)) -> Self {

        Self::from_sites(value.0, value.1,None,None,None)
    }
}
// impl std::fmt::Display for SiteInteraction {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
//         write!(f,"")
//         // match self{
//             // SiteType::A => write!(f,"A"),
//             // SiteType::B => write!(f,"B"),
//             // SiteType::C => write!(f,"C"),
//         // }
//     }
// }

impl std::fmt::Display for SiteType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self{
            SiteType::A => write!(f,"A"),
            SiteType::B => write!(f,"B"),
            SiteType::C => write!(f,"C"),
        }
    }
}

impl std::fmt::Display for Site {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Site(type={},owner={},idx={})",self.t,self.c,self.idx)
    }
}

impl From<Site> for (usize,usize) {

    fn from(value: Site) -> Self {
        
        (value.t(),value.c())
    }
}



#[cfg(test)]
pub mod tests{
    use std::hash::Hash;

    use ndarray::array;
    use serde_json::to_string_pretty;

    use super::*;


    #[test]
    fn test1(){
        let ew = 166.55e2;
        let bw = 0.0692;
        let site1 = Site::new(SiteType::A, 0, 0, ew, bw);
        let site2 = Site::new(SiteType::B, 0, 1, ew, bw);
        let site3 = Site::new(SiteType::B, 1, 2, 0.0, 0.0);
        let sites = vec![site1,site2,site3];
        let mut interactions = SiteInteraction::interactions_from_sites(&sites,HashMap::new());
        
        // println!("no-solv:");
        // for inter in &interactions{

        //     let json = to_string_pretty(inter).unwrap();
        //     println!("{}",json)
        // }

        // dbg!(interactions[0].is_solvation());
        //add co2 solvated site
        let kwco2 = 0.1836;
        let mut solvations = SiteInteraction::solvations_from_ik(0, 1, &sites, None, kwco2);

        interactions.append(&mut solvations);
        
        // println!("solv:");
        for inter in &interactions{

            let json = to_string_pretty(inter).unwrap();
            println!("{}",json)
        }
        let n = interactions.len();
        assert_eq!(n,2);


    }
    
    #[test]
    fn test2(){
        let ew = 166.55e2;
        let bw = 0.0692;
        let eacoh = 403.23e2;
        let bacoh = 4.5e-3;

        let site1 = Site::new(SiteType::A, 0, 0, ew, bw);
        let site2 = Site::new(SiteType::B, 0, 1, ew, bw);
        let site3 = Site::new(SiteType::C, 1, 2, eacoh, bacoh);
        let site4 = Site::new(SiteType::B, 2, 3, 0.0, 0.0);
        
        let sites = vec![site1,site2,site3,site4];
        let mut interactions = SiteInteraction::interactions_from_sites(&sites,HashMap::new());
        
        println!("no-solv:");
        for inter in &interactions{

            let json = to_string_pretty(inter).unwrap();
            println!("{}",json)
        }

        // dbg!(interactions[0].is_solvation());
        //add co2 solvated site
        let kwco2 = 0.1836;
        let mut solvations = SiteInteraction::solvations_from_ik(0, 2, &sites, None, kwco2);

        interactions.append(&mut solvations);
        
        println!("solv:");
        for inter in &interactions{

            let json = to_string_pretty(inter).unwrap();
            println!("{},",json);

        }
        let n = interactions.len();
        // println!("N = {}",n)
        assert_eq!(n,5);

    }

}