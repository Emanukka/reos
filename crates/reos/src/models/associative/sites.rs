use std::marker::PhantomData;

use serde::Serialize;

use crate::models::{IDEAL_GAS_CONST as R, associative::parameters::W};


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

    fn volumetric_factor_jl(fv_j:f64,fv_l:f64)->f64{
        0.5 * (fv_j + fv_l)
    }
    fn elliot_volumetric_factor_jl(fv_j:f64,fv_l:f64)->f64{
        (fv_j * fv_l).sqrt()
    }

    fn compute_volumetric_factor_jl(&self,fv_j:f64,fv_l:f64)->f64;

    fn compute_dimensionless_delta_jl(&self,t:f64,interaction:&SiteInteraction)->f64;

    
    fn association_strength_jl(
        &self,
        t:f64,
        fv_j:f64,
        fv_l:f64,
        interaction:&SiteInteraction)->f64{
        
        let d_jl = self.compute_dimensionless_delta_jl(t, interaction);
        let fv_jl = self.compute_volumetric_factor_jl(fv_j, fv_l);

        fv_jl * d_jl

        }


}
pub struct SiteInteraction{
    pub site_j:Site,
    pub site_l:Site,
    pub epsilon:f64,
    pub kappa: f64,
    pub combining_rule:CombiningRule
}

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
    
    fn compute_dimensionless_delta_jl(&self,t:f64,interaction:&SiteInteraction)->f64 {

        match self {
            CombiningRule::CR1 => {
                let epsilon = interaction.epsilon;
                let kappa= interaction.kappa;
                // (epsilon/R/t).exp_m1() * kappa
                Self::dimensionless_delta_jl(t, epsilon, kappa)

            }

            CombiningRule::ECR => {
                let epsilon_j = interaction.site_j.epsilon();
                let kappa_j = interaction.site_j.kappa();
                let epsilon_l = interaction.site_l.epsilon();
                let kappa_l = interaction.site_l.kappa();

                Self::elliot_dimensionless_delta_jl(t, epsilon_j, kappa_j, epsilon_l, kappa_l)
            }
        }
    }

    fn compute_volumetric_factor_jl(&self,fv_j:f64,fv_l:f64)->f64 {
        
        match self {
            
            CombiningRule::CR1 => {
                Self::volumetric_factor_jl(fv_j, fv_l)
            }

            CombiningRule::ECR => {
                Self::elliot_volumetric_factor_jl(fv_j, fv_l)
            }
        }
    }
}

impl SiteInteraction{

    fn new(site_j:Site,site_l:Site)->Self{

        let epsilon = 0.5*(site_j.epsilon + site_l.kappa);
        let kappa = (site_j.kappa * site_l.kappa).sqrt();
        let combining_rule = CombiningRule::default();

        Self { site_j, site_l, epsilon, kappa, combining_rule }
    }

    fn from_sites(sites:Vec<Site>)->Vec<Self>{

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

    fn dimensionless_delta_jl(&self,t:f64)->f64{
        self.combining_rule.compute_dimensionless_delta_jl(t, self)
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

    use super::Site;


    #[test]
    fn print_sites(){

        // let s1=Site::A(0);
        // let s2=Site::B(0);

        // println!("Site1 = {}",s1);
        // println!("Site2 = {}",s2);


        // let sites=array![s1,s2];

        // println!("Sites = {}",sites)
        
    }
}