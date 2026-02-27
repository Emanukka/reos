
use serde::Serialize;
// allowed site interactions
const W: [[f64;3]; 3] = 
[
    [0.0,1.0,1.0],
    [1.0,0.0,1.0],
    [1.0,1.0,1.0],
];

pub const A: usize = 0;
pub const B: usize = 1;
pub const C: usize = 2;
pub const NS:  usize = 3;
pub const SITES:[usize;3] = [A,B,C];

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
    
    pub fn typ_idx(&self)->usize{
        match self.typ{
            SiteType::A => A,
            SiteType::B => B,
            SiteType::C => C,
        }
    }

    pub fn is_self_associative(&self)->bool{ (self.epsilon != 0.0) && (self.kappa != 0.0) }

    pub fn is_solvate(&self) -> bool { !self.is_self_associative() }

    pub fn cross_associate_with(&self,other:&Self)->bool{

        (W[self.typ_idx()][other.typ_idx()] == 1.0) 
        && (self.is_self_associative()) 
        && (other.is_self_associative()) 
    }

    pub fn can_be_solvated_by(&self,other:&Self)->bool{

        (W[self.typ_idx()][other.typ_idx()] == 1.0) 
        && (self.is_solvate() && other.is_self_associative())

    }

    pub fn new(typ:SiteType, owner:usize, idx:usize, mul:f64, epsilon:f64, kappa:f64)->Self{
        
        Self{typ, owner, idx, mul, epsilon, kappa}
    
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
mod tests{

    use super::super::combining_rule::*;
    use super::super::parameters::AssociationBinaryRecord;
    use super::super::sites_interaction::SitesInteraction;
    use super::super::sites::{Site, SiteType};
    use std::collections::HashMap;

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

        let interactions = ndarray::Array1::from_vec(SitesInteraction::new_interactions(&sites, HashMap::new()));

        
        let reff = ndarray::array![
            SitesInteraction::from_sites(&sites[0], &sites[1], &AssociationBinaryRecord::Set { epsilon: ew, kappa: bw }),
            SitesInteraction::from_sites(&sites[0], &sites[2], &AssociationBinaryRecord::Set { epsilon: 0.5 * (ew + eacoh), kappa: (bw * bacoh).sqrt() }),
            SitesInteraction::from_sites(&sites[1], &sites[2], &AssociationBinaryRecord::Set { epsilon: 0.5 * (ew + eacoh), kappa: (bw * bacoh).sqrt() }),
            SitesInteraction::from_sites(&sites[2], &sites[2], &AssociationBinaryRecord::Set { epsilon: eacoh, kappa: bacoh }),
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
        let interactions = ndarray::Array1::from_vec(SitesInteraction::new_interactions(&sites,binary));

        assert_eq!(interactions.len(), 2);

        let reff = ndarray::array![
            SitesInteraction::from_sites(&sites[0], &sites[1], &AssociationBinaryRecord::Set { epsilon: ew, kappa: bw }),
            SitesInteraction::from_sites(&sites[0], &sites[2], &AssociationBinaryRecord::Set { epsilon: 0.5*ew, kappa: wco2 }),
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

        // let interactions = SitesInteraction::new_interactions(&sites,b);

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