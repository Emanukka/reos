

use std::{collections::HashMap, fmt::write};

use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize, ser::SerializeStruct};

use crate::{models::associative::sites::{ CombiningRule, NS, Site, SiteInteraction, SiteType}, parameters::{Parameters, PureRecord}, state::eos::EosError};

pub const W:[[f64;3];3]=[[0.0,1.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0]];





#[derive(Clone)]
pub struct AssociativeParameters{
    pub multiplicity: Array1<f64>,
    pub sites:        Vec<Site>  ,
    pub interactions: Vec<SiteInteraction>,
    pub id:           Array2<f64>,
    pub lambda_t:       Array2<f64>,
    pub gamma_t:        Array2<f64>,
}

impl AssociativeParameters {
    
    /// Change the all sites pair combining rule that belongs to 2 associative components in the mixture  
    pub fn change_combining_rule(&mut self,i:usize,k:usize,combining_rule:super::sites::CombiningRule){

        let interactions = &mut self.interactions;

        for inter in interactions{

            if inter.belongs_to(i, k){

                inter.change_combining_rule(combining_rule);
            }
        }
    }

    pub fn set_solvation_interaction(&mut self,i:usize,k:usize,epsilon:Option<f64>,beta:f64){

        let interactions = &mut self.interactions;
        let sites = &self.sites;
        let mut solvation = SiteInteraction::solvations_from_ik(i, k, sites, epsilon, beta);

        interactions.append(&mut solvation);
        
        // for inter in interactions{

            
        // }
    }
    // pub fn set_binary_from_owners(&mut self,i:usize,k:usize,epsilon:Option<f64>,beta:Option<f64>) {
        
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

    //                 match beta {
    //                     Some(bet)=>{
    //                         self.beta[(j,l)] = bet;
    //                         self.beta[(l,j)] = bet;
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

impl Parameters <AssociationPureRecord> for AssociativeParameters {
    
    fn from_records(records:Vec<AssociationPureRecord>)->Self{

        let ncomp = records.len();
        let mut mapping:Vec<usize>=Vec::with_capacity(ncomp);

        for (i,record) in records.iter().enumerate(){
            match record.get_type() {
                AssociationType::Inert=>{continue;}
                AssociationType::Associative=>{mapping.push(i)}

            }
        }

        let n = mapping.len();

        let mut sites:Vec<Site> = Vec::with_capacity(NS*n); 
        let mut m:Vec<f64> = Vec::with_capacity(NS); 
        #[allow(non_snake_case)]
        let mut S = 0; //Site counter
        let mut total_sites_per_component = Vec::with_capacity(n);
        for k in 0..n{

            let i = mapping[k]; 
            let record=&records[i];
            
            let na=record.na as f64;
            let nb=record.nb as f64;
            let nc=record.nc as f64;
            
            let mut total_sites_k = 0;

            if na!=0.0{
                sites.push(Site::new(SiteType::A, i, S,record.epsilon,record.beta));

                m.push(na);
                S+=1;
                total_sites_k+=1;
            }
            if nb!=0.0{
                sites.push(Site::new(SiteType::B, i, S,record.epsilon,record.beta));

                m.push(nb);
                S+=1;
                total_sites_k+=1;
            }
            if nc!=0.0{
                sites.push(Site::new(SiteType::C, i, S,record.epsilon,record.beta));


                m.push(nc);
                S+=1;
                total_sites_k+=1;
            }
            total_sites_per_component.push(total_sites_k);            
        }

        let multiplicity = Array1::from_vec(m);
        let ns=sites.len();

        let interactions = SiteInteraction::interactions_from_sites(&sites,todo!()); 


        let mut lambda:Array2<f64>=Array2::zeros((n,S));
        let mut gamma:Array2<f64>=Array2::zeros((ncomp,n));
        let id = Array2::eye(S);

        let mut sumj = 0;
        for k in 0..n{
            
            let jk = total_sites_per_component[k];
            gamma[(mapping[k],k)] = 1.0;

            for j in sumj..sumj+jk{
                lambda[(k,j)] = 1.0
            }
            sumj+=jk

        }
        let lambda_t = lambda.t().to_owned();
        let gamma_t = gamma.t().to_owned();

        AssociativeParameters{
            interactions,
            lambda_t,
            gamma_t,
            multiplicity,
            id,
            sites,
        }
    }
}

// #[derive(Serialize)]
// struct Mat {
//     a: Vec<Vec<f64>>,
// }

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

        // state.serialize_field("P", &array2_to_vec(&self.p))?;
        state.serialize_field("lambdat", &array2_to_vec(&self.lambda_t))?;
        state.serialize_field("gammat", &array2_to_vec(&self.gamma_t))?;
        // state.serialize_field("epsilon", &array2_to_vec(&self.epsilon))?;
        // state.serialize_field("beta", &array2_to_vec(&self.beta))?;

        state.end()        
    }
}


impl std::fmt::Display for AssociativeParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(
            f, 
            "sites:{:?}\nmultiplicity:{:?}\nlambdat:\n{}\ngammat:\n{}\n" ,
            &self.sites,
            &self.multiplicity.to_vec(),
            // &array2_to_string(&self.p),
            &array2_to_string(&self.lambda_t),
            &array2_to_string(&self.gamma_t))
            // &array2_to_string(&self.epsilon),
            // &array2_to_string(&self.beta),)


        // write!(f, "multiplicity:{:#?}",&self.multiplicity.to_vec());
        // write!(f, "P:{:#?}",&array2_to_string(&self.p))


        // match self{
        //     SiteType::A => write!(f,"A"),
        //     SiteType::B => write!(f,"B"),
        //     SiteType::C => write!(f,"C"),
        // }
    }
}
#[derive(Serialize, Deserialize,Debug,Clone,PartialEq)]
#[serde(rename_all = "lowercase")]

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
    pub beta:f64,
    #[serde(default)]
    pub na:usize,
    #[serde(default)]
    pub nb:usize,
    #[serde(default)]
    pub nc:usize,
    // #[serde(default)]
    // pub typ:AssociationType,

}
impl PureRecord for AssociationPureRecord {}

#[derive(Clone)]
pub struct AssociationBinaryRecord{
    pub comp1:usize,
    pub comp2:usize,
    pub epsilon:Option<f64>,
    pub kappa:Option<f64>,
    pub combining_rule:Option<CombiningRule>
}


impl AssociationBinaryRecord {
    
    pub fn new(comp1:usize,comp2:usize,epsilon:Option<f64>,kappa:Option<f64>,combining_rule:Option<CombiningRule>)->Self{

        Self { comp1, comp2, epsilon, kappa, combining_rule }
    }

    pub fn map_from_vec(binary:Vec<Self>)->HashMap<(usize,usize),Self>{

        let mut map:HashMap<(usize,usize), Self> = HashMap::new();

        binary.iter().map(
            |v| {
                map.insert((v.comp1,v.comp2), v.clone());
            }
        );
        map
    }
}
impl AssociationPureRecord {
    
    pub fn new(
        epsilon:f64,
        beta:f64,
        na:usize,
        nb:usize,
        nc:usize,
    )->Self{
        Self{
        epsilon,
        beta,
        na,
        nb,
        nc,
        // typ,
        }
    }
    pub fn inert() -> Self {
        AssociationPureRecord{
        epsilon:0.0,
        beta:0.0,
        na:0,
        nb:0,
        nc:0,
        }
    }

    pub fn solvate(sites:[usize;NS])-> Self{
        AssociationPureRecord{
        epsilon:0.0,
        beta:0.0,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        }
    }

    pub fn associative(epsilon:f64,beta:f64,sites:[usize;NS])->Self{
        AssociationPureRecord{
        epsilon,
        beta,
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


pub mod readyto{
    use crate::models::associative::parameters::{AssociationBinaryRecord, AssociationPureRecord, AssociativeParameters};

    pub fn water4c_cpa()->AssociationPureRecord{

        AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0])
    }
    pub fn acetic1a_cpa()->AssociationPureRecord{
        AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1])
    }

    // pub fn binary_water_acetic()->AssociationBinaryRecord{

    // }
//     pub fn water_acetic_acid()->E<SCPA>{
//                 //Records
//             //1:Water, 2:Acetic Acid
//             let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
//             let c2=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);

//             let a1=AssociationPureRecord::associative(
//                 166.55e2, 
//                 0.0692, 
//                 [2,2,0],
//             );
//             let a2=AssociationPureRecord::associative(
//                 403.23e2, 
//                 4.5e-3, 
//                 [0,0,1],
//             );


//             let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
//             let mut parameters = CPAParameters::from_records(records);
//             parameters.assoc.change_combining_rule(0, 1, sites::CombiningRule::ECR);
//             parameters.cubic.set_kij(0, 1, -0.222);
//             let cpa = SCPA::from_parameters(parameters);
//             E::from_residual(cpa)

//     } 
// // } 
// pub fn acetic_acid_water()->E<SCPA>{
//             //Records
//         //1:Water, 2:Acetic Acid
//         let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
//         let c2=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);

//         let a1=AssociationPureRecord::associative(
//             166.55e2, 
//             0.0692, 
//             [2,2,0],
//         );
//         let a2=AssociationPureRecord::associative(
//             403.23e2, 
//             4.5e-3, 
//             [0,0,1],
//         );

//         let records = vec![CPAPureRecord::new(c2, a2),CPAPureRecord::new(c1, a1)];
//         let mut parameters = CPAParameters::from_records(records);
//         parameters.assoc.change_combining_rule(0, 1, sites::CombiningRule::ECR);
//         parameters.cubic.set_kij(0, 1, -0.222);
//         let cpa = SCPA::from_parameters(parameters);
//         E::from_residual(cpa)

// } 

// pub fn water_co2()->E<SCPA>{
//             //Records
//         let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
//         let c2=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

//         let a1=AssociationPureRecord::associative(
//             166.55e2, 
//             0.0692, 
//             [2,2,0],
//         );
//         let a2=AssociationPureRecord::solvate(
//             [0,1,0]);

//         let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
//         let mut parameters = CPAParameters::from_records(records);
//         // parameters.cubic.set_kij(0, 1, -0.222);

//         parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
//         parameters.assoc.set_solvation_interaction(0, 1, None, 0.1836);

//         let cpa = SCPA::from_parameters(parameters);
//         //Create new State
//         E::from_residual(cpa)
// } 
// pub fn co2_water()->E<SCPA>{
//             //Records
//         //1:Water, 2:Acetic Acid
//         let c1=CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);
//         let c2=CubicPureRecord::new_set1(0.35079, 0.0272e-3, 0.7602, 304.12);

//         let a1=AssociationPureRecord::associative(
//             166.55e2, 
//             0.0692, 
//             [2,2,0],
//         );
//         let a2=AssociationPureRecord::solvate(
//             [0,1,0],
//         );

//         let records = vec![CPAPureRecord::new(c2, a2),CPAPureRecord::new(c1, a1)];
//         let mut parameters = CPAParameters::from_records(records);

//         parameters.cubic.set_kij_temperature_dependent(0, 1, -0.15508, 0.000877);
//         // parameters.assoc.set_binary_from_owners(0, 1, None, Some(0.1836));

//         parameters.assoc.set_solvation_interaction(0, 1, None, 0.1836);

//         let cpa = SCPA::from_parameters(parameters);
//         //Create new State
//         E::from_residual(cpa)

// } 


// pub fn methanol_2b()->E<SCPA>{
//             //Records
//         //1:metoh, 2:oct
//         let c1=CubicPureRecord::new_set1(0.40531, 0.0000309, 0.4310, 513.);

//         let a1=AssociationPureRecord::associative(
//             24591.0, 
//             0.01610, 
//             [1,1,0],
//         );
//         let records = vec![CPAPureRecord::new(c1, a1)];
//         let parameters = CPAParameters::from_records(records);


//         let cpa = SCPA::from_parameters(parameters);
//         //Create new State
//         E::from_residual(cpa)

        
// } 
// pub fn methanol_3b()->E<SCPA>{
//             //Records
//         //1:metoh, 2:oct
//         let c1=
//         CubicPureRecord::new_set1(
//             4.5897e-1, 
//             0.0334e-3, 
//             1.0068,
//             513.);

//         let a1=AssociationPureRecord::associative(
//             160.70e2, 
//             34.4e-3, 
//             [2,1,0],
//         );

//         let records = vec![CPAPureRecord::new(c1, a1)];
//         let parameters = CPAParameters::from_records(records);


//         let cpa = SCPA::from_parameters(parameters);
//         E::from_residual(cpa)
// } 


// pub fn acoh_octane()->E<SCPA>{
//             //Records
//         //1:Acetic Acid, 2:Octane
//         let c1=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);
//         let c2=CubicPureRecord::new_set1(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

//         let a1=AssociationPureRecord::associative(
//             403.23e2, 
//             4.5e-3, 
//             [0,0,1],
//         );
//         let a2=AssociationPureRecord::inert();


//         let records = vec![CPAPureRecord::new(c1, a1),CPAPureRecord::new(c2, a2)];
//         let parameters = CPAParameters::from_records(records);


//         let mut cpa = SCPA::from_parameters(parameters);
//         cpa.cubic.parameters.set_kij(0, 1, 0.064);

//         E::from_residual(cpa)

        
// } 
// pub fn octane_acoh()->E<SCPA>{
//             //Records
//         //1:Acetic Acid, 2:Octane
//         let c1=CubicPureRecord::new_set1(0.91196, 0.0468e-3, 0.4644, 594.8);
//         let c2=CubicPureRecord::new_set1(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

//         let a1=AssociationPureRecord::associative(
//             403.23e2, 
//             4.5e-3, 
//             [0,0,1],
//         );
//         let a2=AssociationPureRecord::inert();


//         let records = vec![CPAPureRecord::new(c2, a2),CPAPureRecord::new(c1, a1)];
//         let parameters = CPAParameters::from_records(records);


//         let mut cpa = SCPA::from_parameters(parameters);

//         cpa.cubic.parameters.set_kij(0, 1, 0.064);

//         E::from_residual(cpa)

        
// } 


}


#[cfg(test)]
mod tests{

    use serde_json::{from_str, to_string, to_string_pretty};
    use super::*;


    #[test]
    fn test1(){

        let data1 = r#"
        {   
            "epsilon": 166.55e2, 
            "beta": 0.0692,
            "na":   2,
            "nb":   2
        }
        "#;
        
        let data2 = r#"
        {
            "nb":   1
        }
        "#;

        let c1: AssociationPureRecord = from_str(data1).unwrap();
        let c2: AssociationPureRecord = from_str(data2).unwrap();

        // let c1_string = serde_json::to_string_pretty(&c1).unwrap();

        // println!("{}",c1_string)

        let p = AssociativeParameters::from_records(vec![c1,c2]);

        // let string = to_string(&p).unwrap();
        let string = p.to_string();


        println!("{}",string)

    }
}