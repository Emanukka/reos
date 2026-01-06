
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize, ser::SerializeStruct};

use crate::{models::associative::sites::{ CombiningRule, NS, Site, SiteInteraction, SiteType}, parameters::{Parameters, records::{BinaryParameter, Properties}, }, state::eos::EosError};






#[derive(Clone)]
pub struct AssociativeParameters{
    pub multiplicity: Array1<f64>,
    pub sites:        Vec<Site>  ,
    pub interactions: Vec<SiteInteraction>,
    pub na:           usize,
    pub nb:           usize,
    pub nc:           usize,
    pub map:          Array1<usize>
}

impl AssociativeParameters {
    
    /// Change the all sites pair combining rule that belongs to 2 associative components in the mixture  
    pub fn change_combining_rule(&mut self,i:usize,k:usize,combining_rule:super::sites::CombiningRule){

        let interactions = &mut self.interactions;

        for inter in interactions{

            if inter.belongs_to(i, k, self.map.as_slice().unwrap()){

                inter.change_combining_rule(combining_rule);
            }
        }
    }

    // pub fn set_solvation_interaction(&mut self,i:usize,k:usize,epsilon:Option<f64>,kappa:f64){

    //     let interactions = &mut self.interactions;
    //     let sites = &self.sites;
    //     let mut solvation = SiteInteraction::solvations_from_ik(i, k, sites, epsilon, kappa);

    //     interactions.append(&mut solvation);
        
    //     // for inter in interactions{

            
    //     // }
    // }
    // pub fn set_binary_from_owners(&mut self,i:usize,k:usize,epsilon:Option<f64>,kappa:Option<f64>) {
        
    //     let mut hit = 0;
    //     for (j,site_j) in self.sites.iter().enumerate(){
    //         for (l,site_l) in self.sites.iter().enumerate(){
    //             hit+=1;
    //             if (site_j.has_owner(i)) && (site_l.has_owner(k)){

    //                 match epsilon {
    //                     Some(eps)=>{
    //                         self.epsilon[(j,l)] = eps;
    //                         self.epsilon[(l,j)] = eps
    //                     }
    //                     None => {}
    //                 }

    //                 match kappa {
    //                     Some(bet)=>{
    //                         self.kappa[(j,l)] = bet;
    //                         self.kappa[(l,j)] = bet;
    //                     }
    //                     None => {}
    //                 }
    //             }
    //         }
    //     }

    //     if hit == 0{
    //         println!("no site-pair was match")
    //     }
    // }
}

impl Parameters<AssociationPureRecord,AssociationBinaryRecord> for AssociativeParameters {

    fn raw(records: Vec<AssociationPureRecord>, binary: Vec<BinaryParameter<AssociationBinaryRecord>>, _properties: Option<Properties>)->Self {
        
        let ncomp = records.len();
        let mut mapping:Vec<usize>=Vec::with_capacity(ncomp);

        for (i,record) in records.iter().enumerate(){

            match record.get_type() {

                AssociationType::Associative=>{mapping.push(i)}
                _ => continue

            }
        }

        let n = mapping.len();

        let mut n_a = 0;
        let mut n_b = 0; 
        let mut n_c = 0; 

        let mut sites:Vec<Site> = Vec::with_capacity(NS * n); 
        let mut m:Vec<f64> = Vec::with_capacity(NS); 
        let mut s = 0;
        let mut total_sites_per_component = Vec::with_capacity(n);
        let mut map:Vec<usize> = Vec::with_capacity( NS * n);

        for k in 0..n{

            let i = mapping[k]; 
            let record=&records[i];
            
            let na= record.na as f64;
            let nb= record.nb as f64;
            let nc= record.nc as f64;
            
            let mut total_sites_k = 0;

            if na!=0.0{
                sites.push(Site::new(SiteType::A, i, s,record.epsilon,record.kappa));
                map.push(i);
                n_a += 1;
                m.push(na);
                s+=1;
                total_sites_k+=1;
            }
            if nb!=0.0{
                sites.push(Site::new(SiteType::B, i, s,record.epsilon,record.kappa));
                map.push(i);
                n_b += 1;
                m.push(nb);
                s+=1;
                total_sites_k+=1;
            }
            if nc!=0.0{
                sites.push(Site::new(SiteType::C, i, s,record.epsilon,record.kappa));
                map.push(i);
                n_c += 1;
                m.push(nc);
                s+=1;
                total_sites_k+=1;
            }
            total_sites_per_component.push(total_sites_k);            
        }

        let multiplicity = Array1::from_vec(m);
        let map = Array1::from_vec(map);
        let interactions = SiteInteraction::interactions_from_sites(&sites, binary); 

        AssociativeParameters{
            na: n_a,
            nb: n_b,
            nc: n_c,
            map,
            interactions,
            multiplicity,
            sites,
        }
    }

    fn get_properties(&self)->&Properties { panic!() }
    

}

fn array2_to_vec(a: &Array2<f64>) -> Vec<Vec<f64>> {
    a.rows()
     .into_iter()
     .map(|r| r.to_vec())
     .collect()
}

fn array2_to_string(a: &Array2<f64>) -> String {
    let rows: Vec<String> = a
        .rows()
        .into_iter()
        .map(|r| {
            let elems: Vec<String> =
                r.iter().map(|x| format!("{}", x)).collect();
            format!("[{}]", elems.join(", "))
        })
        .collect();

    format!("[{}]", rows.join(",\n"))
}

impl Serialize for AssociativeParameters {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: serde::Serializer {
            
        let fields = 7;
        let mut state = serializer.serialize_struct("AssociativeParameters", fields)?;
        
        state.serialize_field("sites", &self.sites)?;

        state.serialize_field("multiplicity", &self.multiplicity.to_vec())?;

        state.serialize_field("interactions", &self.interactions)?;


        // state.serialize_field("P", &array2_to_vec(&self.p))?;
        // state.serialize_field("lambdat", &array2_to_vec(&self.lambda_t))?;
        // state.serialize_field("gammat", &array2_to_vec(&self.gamma_t))?;
        // state.serialize_field("epsilon", &array2_to_vec(&self.epsilon))?;
        // state.serialize_field("kappa", &array2_to_vec(&self.kappa))?;

        state.end()        
    }
}


impl std::fmt::Display for AssociativeParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(
            f, 
            "sites:{:?}\nmultiplicity:{:?}\ninteractions:{:?}" ,
            &self.sites,
            &self.multiplicity.to_vec(),
            &self.interactions
            // &array2_to_string(&self.p),
            // &array2_to_string(&self.lambda_t),
            // &array2_to_string(&self.gamma_t))
            // &array2_to_string(&self.epsilon),
            // &array2_to_string(&self.kappa),)
        )

        // write!(f, "multiplicity:{:#?}",&self.multiplicity.to_vec());
        // write!(f, "P:{:#?}",&array2_to_string(&self.p))


        // match self{
        //     SiteType::A => write!(f,"A"),
        //     SiteType::B => write!(f,"B"),
        //     SiteType::C => write!(f,"C"),
        // }
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

#[derive(Clone,Serialize,Deserialize,Debug)]
pub struct AssociationBinaryRecord{
    #[serde(default)]
    pub epsilon:Option<f64>,
    #[serde(default)]
    pub kappa:Option<f64>,
    #[serde(default)]
    pub combining_rule:Option<CombiningRule>
}


impl AssociationBinaryRecord {
    
    pub fn new(epsilon:Option<f64>,kappa:Option<f64>,combining_rule:Option<CombiningRule>)->Self{

        Self { epsilon, kappa, combining_rule }
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

    use serde_json::{from_str, to_string, to_string_pretty};
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

        let p = AssociativeParameters::new(vec![pr1,pr2], vec![br]);

        let string = to_string(&p).unwrap();

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