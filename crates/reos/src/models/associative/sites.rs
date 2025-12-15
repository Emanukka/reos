use std::{marker::PhantomData, path::Iter};

use serde::Serialize;

use crate::{models::{IDEAL_GAS_CONST as R, associative::parameters::W}, state::E};


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
    pub j:usize,

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
    
    /// false in the case when is solvate (since the site exists)
    pub fn is_associative(&self)->bool{

        if (self.epsilon() != 0.0) && (self.kappa() != 0.0) { true }
        else { false }
    
    }
    pub fn interacts_with(&self,other:&Self)->bool{

        
        if (W[self.t()][other.t()] == 1.0) 
        && (self.is_associative()) 
        && (other.is_associative()) { true }
        else { false }

    }
    pub fn new(t:SiteType,c:usize,j:usize,epsilon:f64,kappa:f64)->Self{
        Self{
            t,
            c,
            j,
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

    fn cr1_volumetric_factor_jl(fv_j:f64,fv_l:f64)->f64{
        0.5 * (fv_j + fv_l)
    }
    fn elliot_volumetric_factor_jl(fv_j:f64,fv_l:f64)->f64{
        (fv_j * fv_l).sqrt()
    }


    fn cpa_volumetric_factor_jl(g:f64,bjl:f64)->f64{
        g * bjl
    }

    fn saft_volumetric_factor_jl(gjl:f64,djl:f64)->f64{
        gjl * djl.powi(3)
    }
    
    fn association_strength_jl(
        &self,
        t:f64,
        fv_j:f64,
        fv_l:f64,
        interaction:&SiteInteraction)->f64;


}
#[derive(Serialize)]
pub struct SiteInteraction{
    site_j:Site,
    site_l:Site,
    epsilon:f64,
    kappa: f64,
    pub combining_rule:CombiningRule
}
#[derive(Serialize)]
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
        fv_j:f64,
        fv_l:f64,
        interaction:&SiteInteraction)->f64 {
        match self {

            CombiningRule::CR1 => {

                let epsilon = interaction.epsilon;
                let kappa= interaction.kappa;
                // (epsilon/R/t).exp_m1() * kappa
                let d_jl = Self::dimensionless_delta_jl(t, epsilon, kappa);
                let fv_jl = Self::cr1_volumetric_factor_jl(fv_j,fv_l) ;
                
                fv_jl * d_jl

            }

            CombiningRule::ECR => {
                let epsilon_j = interaction.site_j.epsilon();
                let kappa_j = interaction.site_j.kappa();
                let epsilon_l = interaction.site_l.epsilon();
                let kappa_l = interaction.site_l.kappa();

                let d_jl = Self::elliot_dimensionless_delta_jl(t, epsilon_j, kappa_j, epsilon_l, kappa_l);
                let fv_jl = Self::elliot_volumetric_factor_jl(fv_j, fv_l);

                fv_jl * d_jl

            }
        }
    }

}

impl SiteInteraction{

    pub fn get_eps(&self)->f64{
        self.epsilon
    }

    pub fn get_kappa(&self)->f64{
        self.kappa
    }

    pub fn new(site_j:Site,site_l:Site)->Self{

        let epsilon = 0.5 * (site_j.epsilon + site_l.epsilon);
        let kappa = (site_j.kappa * site_l.kappa).sqrt();
        let combining_rule = CombiningRule::default();

        Self { site_j, site_l, epsilon, kappa, combining_rule }
    }

    pub fn from_sites(sites:Vec<Site>)->Vec<Self>{

        let mut interactions = Vec::<SiteInteraction>::new();

        for site_j in &sites{
            for site_l in &sites{
                
                if site_j.interacts_with(site_l){
                    let interaction = SiteInteraction::new(site_j.clone(), site_l.clone());
                    interactions.push(interaction);
                }
            }
        }
        interactions
    }
    
    pub fn add_interaction_from_components_ik(
        &mut self,
        sites:&Vec<Site>,
        i:usize,
        k:usize,
        epsilon:f64,
        kappa:f64){
        
        

    }
    pub fn association_strength_jl(
        &self,
        t:f64,
        fv_j:f64,
        fv_l:f64,)->f64{
        
        self.combining_rule.association_strength_jl(t, fv_j, fv_l, self)
        
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
    
    pub fn change_cross_parameters(&mut self,epsilon:Option<f64>,kappa:Option<f64>){

        self.epsilon = epsilon.unwrap_or(self.get_eps());

        self.kappa = kappa.unwrap_or(self.get_kappa());

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
        write!(f, "Site(type={},owner={},idx={})",self.t,self.c,self.j)
    }
}

impl From<Site> for (usize,usize) {

    fn from(value: Site) -> Self {
        
        (value.t(),value.c())
    }
}



#[cfg(test)]
pub mod tests{
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
        
        let interactions = SiteInteraction::from_sites(vec![site1,site2,site3]);
        
        
        for inter in interactions{
            let json = to_string_pretty(&inter).unwrap();
            println!("{}",json)
        }

        //add co2 solvated site
        // let json = interactions.iter().map(|x| to_string_pretty(x).unwrap());

        
        // let data1 = r#"
        // {   
        //     "epsilon": 166.55e2, 
        //     "beta": 0.0692,
        //     "na":   2,
        //     "nb":   2
        // }
        // "#;
        
        // let data2 = r#"
        // {
        //     "nb":   1
        // }
        // "#;

        // let c1: AssociationPureRecord = from_str(data1).unwrap();
        // let c2: AssociationPureRecord = from_str(data2).unwrap();

        // // let c1_string = serde_json::to_string_pretty(&c1).unwrap();

        // // println!("{}",c1_string)

        // let p = AssociativeParameters::from_records(vec![c1,c2]);

        // // let string = to_string(&p).unwrap();
        // let string = p.to_string();


        // println!("{}",string)

        // let s1=Site::A(0);
        // let s2=Site::B(0);

        // println!("Site1 = {}",s1);
        // println!("Site2 = {}",s2);


        // let sites=array![s1,s2];

        // println!("Sites = {}",sites)
        
    }
}