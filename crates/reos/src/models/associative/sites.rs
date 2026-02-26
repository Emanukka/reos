use std::{collections::HashMap};

use serde::{Deserialize, Serialize};
use super::parameters::{AssociationBinaryRecord};
use super::strength;

// allowed site interactions
const W:[[f64;3];3]=[
    [0.0,1.0,1.0],
    [1.0,0.0,1.0],
    [1.0,1.0,1.0]];

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


    pub fn has_owner(&self,i:usize)->bool{
        i == self.owner
    }
    
    pub fn is_self_associative(&self)->bool{

        if (self.epsilon != 0.0) && (self.kappa != 0.0) { true }
        else { false }
    
    }
    pub fn is_solvate(&self)->bool{

        if !self.is_self_associative(){ true }
        else { false }
    
    }

    pub fn cross_associate_with(&self,other:&Self)->bool{

        
        if (W[self.typ_idx()][other.typ_idx()] == 1.0) 
        && (self.is_self_associative()) 
        && (other.is_self_associative()) { true }
        else { false }
    }

    pub fn can_be_solvated_by(&self,other:&Self)->bool{

        if (W[self.typ_idx()][other.typ_idx()] == 1.0) 
        && (self.is_solvate() && other.is_self_associative())
        { true }
        else { false }
    }

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


#[derive(Serialize, Clone, Debug, PartialEq)]
pub struct SiteInteraction{
    pub site_j:usize, 
    pub site_l:usize,
    pub interaction_params:InteractionParams
}


#[derive(Serialize,Clone,Copy,PartialEq,Debug,Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum CombiningRuleOption{

    CR1,
    MCR1{kappa:f64},
    ECR
}

#[derive(Serialize,Clone,Copy,PartialEq,Debug,Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InteractionParams{

    EpsilonKappa{epsilon:f64, kappa:f64},
    ECR
}


impl std::fmt::Display for CombiningRuleOption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self{
            CombiningRuleOption::CR1 => write!(f,"CR1"),
            CombiningRuleOption::MCR1{kappa} => write!(f,"m-CR1(kappa={kappa}"),
            CombiningRuleOption::ECR => write!(f,"ECR"),
        }
    }
}

impl std::fmt::Display for InteractionParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self{
            InteractionParams::EpsilonKappa{epsilon, kappa} => write!(f,"epsilon={}, kappa={}",epsilon, kappa),
            InteractionParams::ECR => write!(f,"ECR"),
        }
    }
}

impl Default for CombiningRuleOption {
    fn default() -> Self {
        CombiningRuleOption::CR1
    }
}

impl Into<CombiningRuleOption> for Option<CombiningRuleOption> {

    fn into(self) -> CombiningRuleOption {

        self.unwrap_or_default()
    }
}

impl SiteInteraction {

    pub fn new_interactions(sites:&Vec<Site>, binary: HashMap<(usize,usize), AssociationBinaryRecord>)->Vec<Self>{

        let mut interactions = Vec::<SiteInteraction>::new();
        let s = sites.len();
        
        for j in 0..s{
            for l in j..s{
                
                let site_j = &sites[j];
                let site_l = &sites[l];
                let i = site_j.owner;
                let k = site_l.owner;

                let opt = binary.get(&(i, k)).or(binary.get(&(k, i)));                
    
                if site_j.cross_associate_with(site_l){
                    
                    interactions.push(
                        Self::from_sites(site_j, site_l, &opt.cloned().unwrap_or_default())
                    );

                } else if site_j.can_be_solvated_by(site_l) || site_l.can_be_solvated_by(site_j) {
                    
                    match opt {
                        Some(binary) => 
                            interactions.push(
                                Self::from_sites(site_j, site_l, binary)),
                        
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
        
        match &self.interaction_params {
            
            &InteractionParams::EpsilonKappa { epsilon, kappa } => {

                strength::cr1_association_strength_jl(t, f_ii, f_kk, epsilon, kappa)

            }

            InteractionParams::ECR => {

                let j = self.site_j;
                let l = self.site_l;

                let sj = &sites[j];
                let sl = &sites[l];

                let epsilon_j = sj.epsilon;
                let epsilon_l = sl.epsilon;
                
                let kappa_j = sj.kappa;
                let kappa_l = sl.kappa;

                strength::ecr_association_strength_jl(t, f_ii, f_kk, epsilon_j, epsilon_l, kappa_j, kappa_l)
                // // strength::cr1_association_strength_jl(t, f_ii, f_kk, self.epsilon, self.kappa)
                // let f = strength::elliot_factor_ik(f_ii, f_kk) ;
                // println!("caso1");
                // println!("f = {f}");

                // let d = strength::association_strength_jl(t, self.epsilon, self.kappa, f);

                // println!("DELTA = {d}");
                // d
            }
        }

    }

    pub fn association_strength_jl_dt(
        &self,
        t:f64,
        f_ii:f64,
        f_kk:f64,
        sites:&[Site])->f64{
        
        
        match &self.interaction_params {
            
            &InteractionParams::EpsilonKappa{epsilon, kappa} => {

                strength::cr1_association_strength_jl_dt(t, f_ii, f_kk, epsilon, kappa)

            }

            InteractionParams::ECR => {

                let j = self.site_j;
                let l = self.site_l;

                let sj = &sites[j];
                let sl = &sites[l];

                let epsilon_j = sj.epsilon;
                let epsilon_l = sl.epsilon;
                
                let kappa_j = sj.kappa;
                let kappa_l = sl.kappa;

                strength::ecr_association_strength_jl_dt(t, f_ii, f_kk, epsilon_j, epsilon_l, kappa_j, kappa_l)

            }
        }

    }
}

impl SiteInteraction{


    pub fn from_sites(site_j:&Site, site_l:&Site, binary: &AssociationBinaryRecord)->Self{

        match binary {
            
            AssociationBinaryRecord ::CombiningRule(combining_rule) => {
            
                
                match combining_rule {


                    CombiningRuleOption::CR1 => {
                        
                        let [epsilon_j, epsilon_l] = [site_j.epsilon, site_l.epsilon];
                        let [kappa_j, kappa_l] = [site_j.kappa, site_l.kappa];

                        let epsilon = 0.5 * (epsilon_j + epsilon_l);
                        let kappa = (kappa_j * kappa_l).sqrt();

                        let interaction_params = InteractionParams::EpsilonKappa { epsilon, kappa };
                        
                        Self { site_j: site_j.idx, site_l: site_l.idx, interaction_params }

                    }

                    &CombiningRuleOption::MCR1 { kappa } => {

                        let [epsilon_j, epsilon_l] = [site_j.epsilon, site_l.epsilon];
                        let epsilon = 0.5 * (epsilon_j + epsilon_l);

                        let interaction_params = InteractionParams::EpsilonKappa { epsilon, kappa };
                        
                        Self { site_j: site_j.idx, site_l: site_l.idx, interaction_params }
                    }

                    CombiningRuleOption::ECR => {

                        let interaction_params = InteractionParams::ECR;
                        
                        Self { site_j: site_j.idx, site_l: site_l.idx, interaction_params }
                        
                    }
                }
            }

            &AssociationBinaryRecord::Set { epsilon, kappa } => {

                let interaction_params = InteractionParams::EpsilonKappa { epsilon, kappa };
                Self { site_j: site_j.idx, site_l: site_l.idx, interaction_params }

            }

        }



    }


    // pub fn change_combining_rule(&mut self,combining_rule:CombiningRuleOption){

    //     self.combining_rule = combining_rule;
        
    // }
    // pub fn change_cross_parameters(&mut self,epsilon:Option<f64>,kappa:Option<f64>){

    //     self.epsilon = epsilon.unwrap_or(self.epsilon);

    //     self.kappa = kappa.unwrap_or(self.kappa);

    // } 
    // pub fn belongs_to(&self,i:usize,k:usize,map:&[usize])->bool{


    //     let j = self.site_j;
    //     let l = self.site_l;

    //     let owner_j = map[j];
    //     let owner_l = map[l];

    //     let owners = (owner_j,owner_l);

    //     if owners == (i,k){
    //         true
    //     } else if owners == (k,i){
    //         true
    //     } else {
    //         false  
    //     }
    // }
    
       
}

// impl From<(Site,Site)> for SiteInteraction {
//     fn from(value: (Site,Site)) -> Self {

//         Self::new(value.0, value.1,None,None,CombiningRuleOption::default())
//     }
// }

impl std::fmt::Display for Site {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Site(type={}, owner={}, idx={}, mul={}, epsilon={}, kappa={})",self.typ,self.owner,self.idx,self.mul,self.epsilon,self.kappa)
    }
}
impl std::fmt::Display for SiteInteraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SiteInteraction(j={}, l={}, params={})",self.site_j,self.site_l,self.interaction_params)
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


#[cfg(test)]
pub mod tests{

    use super::*;

    #[test]
    fn test_cross_association() {

        let [ew, bw] = [166.55e2, 0.0692];
        let [eacoh, bacoh] = [403.23e2, 4.5e-3];
        // let wco2 = 0.1836;

        let site1 = Site::new(SiteType::A, 0, 0,2., ew, bw);
        let site2 = Site::new(SiteType::B, 0, 1,2., ew, bw);
        let site3 = Site::new(SiteType::C, 1, 2,1., eacoh, bacoh);

        // let b1 = AssociationBinaryRecord::new(None, None, CombiningRuleOption::CR1);

        let sites = vec![site1, site2, site3];
        
        // let binary = HashMap::from([
        //     ((0, 1), b1)
        // ]);

        let interactions = ndarray::Array1::from_vec(SiteInteraction::new_interactions(&sites, HashMap::new()));

        
        let reff = ndarray::array![
            SiteInteraction::from_sites(&sites[0], &sites[1], &AssociationBinaryRecord::Set { epsilon: ew, kappa: bw }),
            SiteInteraction::from_sites(&sites[0], &sites[2], &AssociationBinaryRecord::Set { epsilon: 0.5 * (ew + eacoh), kappa: (bw * bacoh).sqrt() }),
            SiteInteraction::from_sites(&sites[1], &sites[2], &AssociationBinaryRecord::Set { epsilon: 0.5 * (ew + eacoh), kappa: (bw * bacoh).sqrt() }),
            SiteInteraction::from_sites(&sites[2], &sites[2], &AssociationBinaryRecord::Set { epsilon: eacoh, kappa: bacoh }),
        ];

        assert_eq!(interactions.len(), reff.len());
        assert_eq!(interactions, reff);

    }

    #[test]
    fn test_solvation(){
        
        let [ew, bw] = [166.55e2, 0.0692];

        let wco2 = 0.1836;

        let site1 = Site::new(SiteType::A, 0, 0, 2.,ew, bw);
        let site2 = Site::new(SiteType::B, 0, 1, 2.,ew, bw);
        let site3 = Site::new(SiteType::B, 1, 2, 1.,0.0, 0.0);

        let sites = vec![site1,site2,site3];
        
        let b = AssociationBinaryRecord::CombiningRule(CombiningRuleOption::MCR1 { kappa: wco2 });
        
        // let b1 = AssociationBinaryRecord::new(None, None, CombiningRuleOption::CR1);

        
        let binary = HashMap::from([
            ((0, 1), b)
        ]);
        let interactions = ndarray::Array1::from_vec(SiteInteraction::new_interactions(&sites,binary));

        assert_eq!(interactions.len(), 2);

        let reff = ndarray::array![
            SiteInteraction::from_sites(&sites[0], &sites[1], &AssociationBinaryRecord::Set { epsilon: ew, kappa: bw }),
            SiteInteraction::from_sites(&sites[0], &sites[2], &AssociationBinaryRecord::Set { epsilon: 0.5*ew, kappa: wco2 }),
        ];

        
        assert_eq!(interactions.len(), reff.len());
        assert_eq!(interactions, reff);

        // for inter in &interactions{

        //     let json = to_string_pretty(inter).unwrap();
        //     println!("{}",json)
        // }


    }
    
    #[test]
    fn test_ternary_with_solvation(){
        // let ew = 166.55e2;
        // let bw = 0.0692;
        // let eacoh = 403.23e2;
        // let bacoh = 4.5e-3;
        // let wco2 = 0.1836;

        // let site1 = Site::new(SiteType::A, 0, 0,2., ew, bw);
        // let site2 = Site::new(SiteType::B, 0, 1,2., ew, bw);
        // let site3 = Site::new(SiteType::C, 1, 2,1., eacoh, bacoh);
        // let site4 = Site::new(SiteType::B, 2, 3,1., 0.0, 0.0);
        // let b2 = AssociationBinaryRecord::new(None, Some(wco2), CombiningRuleOption::default());
        // let b1 = AssociationBinaryRecord::new(None, None, CombiningRuleOption::ECR);
        
        // let sites = vec![site1,site2,site3,site4];
        // let b = vec![BinaryParameter::new(b1, 0, 1) , BinaryParameter::new(b2, 0, 2)]
        //     .into_iter()
        //     .map(|r| ((r.id1,r.id2), r.model_record))
        //     .collect();

        // let interactions = SiteInteraction::new_interactions(&sites,b);

        // let i1 = &interactions[0];
        // let i2 = &interactions[1];
        // let i3 = &interactions[2];
        // let i4 = &interactions[3];

        // let water_params = InteractionParams::CR1 { epsilon: ew, kappa: bw };
        // let acoh_params = InteractionParams::CR1 { epsilon: eacoh, kappa: bacoh };
        // let water_acoh_params = InteractionParams::ECR;

        // let water_co2_params = InteractionParams::CR1 { epsilon: 0.5 * ew, kappa: wco2 };


        // for inter in &interactions{

        //     let json = to_string_pretty(inter).unwrap();
        //     println!("{},",json);

        // }
        // let n = interactions.len();
        // assert_eq!(n,5);

        // assert_eq!(i1.interaction_params, water_params);
        // assert_eq!(i4.epsilon, 0.5 * (ew + eacoh));
        // assert_eq!(i4.kappa,(bw * bacoh).sqrt());
        // assert_eq!(i3.epsilon,0.5 * ew);
        // assert_eq!(i3.kappa,wco2);
        // assert_eq!(i2.combining_rule,CombiningRuleOption::ECR);
        // assert_eq!(i4.combining_rule,CombiningRuleOption::ECR);


    }

}