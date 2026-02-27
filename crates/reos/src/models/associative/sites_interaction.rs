use std::collections::HashMap;

use super::sites::Site;
use super::combining_rule::CombiningRuleOption;
use super::parameters::AssociationBinaryRecord;

#[derive(Clone, PartialEq, Debug)]
pub enum InteractionParams{

    EpsilonKappa{epsilon:f64, kappa:f64},
    ECR
}


#[derive(Clone, Debug, PartialEq)]
pub struct SitesInteraction{
    pub site_j:usize, 
    pub site_l:usize,
    pub interaction_params:InteractionParams
}

impl SitesInteraction {

    pub fn new_interactions(sites:&Vec<Site>, binary: HashMap<(usize,usize), AssociationBinaryRecord>)->Vec<Self>{

        let mut interactions = Vec::<SitesInteraction>::new();
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

    pub fn association_strength_jl(&self, t:f64, f_ii:f64, f_kk:f64, sites:&[Site])->f64{
        
        match &self.interaction_params {
            
            &InteractionParams::EpsilonKappa { epsilon, kappa } => {

                super::strength::cr1_association_strength_jl(t, f_ii, f_kk, epsilon, kappa)
            }

            InteractionParams::ECR => {

                let [j,l] = [self.site_j, self.site_l];
                let [sj, sl] = [&sites[j], &sites[l]];

                let [epsilon_j, epsilon_l] = [sj.epsilon, sl.epsilon];
                let [kappa_j, kappa_l] = [sj.kappa, sl.kappa];

                super::strength::ecr_association_strength_jl(t, f_ii, f_kk, epsilon_j, epsilon_l, kappa_j, kappa_l)
            }
        }

    }

    pub fn association_strength_jl_dt(&self, t:f64, f_ii:f64, f_kk:f64, sites:&[Site])->f64{
        
        match &self.interaction_params {
            
            &InteractionParams::EpsilonKappa{epsilon, kappa} => {

                super::strength::cr1_association_strength_jl_dt(t, f_ii, f_kk, epsilon, kappa)

            }

            InteractionParams::ECR => {

                let [j,l] = [self.site_j, self.site_l];
                let [sj, sl] = [&sites[j], &sites[l]];

                let [epsilon_j, epsilon_l] = [sj.epsilon, sl.epsilon];
                let [kappa_j, kappa_l] = [sj.kappa, sl.kappa];

                super::strength::ecr_association_strength_jl_dt(t, f_ii, f_kk, epsilon_j, epsilon_l, kappa_j, kappa_l)

            }
        }

    }
}

impl SitesInteraction{

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

// impl From<(Site,Site)> for SitesInteraction {
//     fn from(value: (Site,Site)) -> Self {

//         Self::new(value.0, value.1,None,None,CombiningRuleOption::default())
//     }
// }

impl std::fmt::Display for InteractionParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self{
            InteractionParams::EpsilonKappa{epsilon, kappa} => write!(f,"epsilon={}, kappa={}",epsilon, kappa),
            InteractionParams::ECR => write!(f,"ECR"),
        }
    }
}


impl std::fmt::Display for Site {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Site(type={}, owner={}, idx={}, mul={}, epsilon={}, kappa={})",self.typ,self.owner,self.idx,self.mul,self.epsilon,self.kappa)
    }
}
impl std::fmt::Display for SitesInteraction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "SitesInteraction(j={}, l={}, params={})",self.site_j,self.site_l,self.interaction_params)
    }
}
