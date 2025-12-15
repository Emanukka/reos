

use std::fmt::write;

use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize, ser::SerializeStruct};

use crate::{models::{ associative::sites::{NS, Site, SiteType}}, parameters::{Parameters, PureRecord}, state::eos::EosError};

pub const W:[[f64;3];3]=[[0.0,1.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0]];



pub trait CombiningRule{

    fn which()->Self where  Self: Sized;

    fn to_string()->String;

    fn association_strength_jl(
        t:f64,
        g_contact:f64,
        bij:&Array2<f64>,
        epsilon:&Array2<f64>,
        beta:&Array2<f64>,
    )->Array2<f64>;

}



#[derive(Clone)]
pub struct AssociativeParameters{
    pub multiplicity: Array1<f64>,
    pub p:            Array2<f64>,
    pub epsilon:      Array2<f64>,
    pub beta:         Array2<f64>,
    pub sites:        Vec<Site>  ,
    pub lambda:       Array2<f64>,
    pub gamma:        Array2<f64>,
    pub lambda_t:     Array2<f64>,
    pub gamma_t:      Array2<f64>,
}

impl AssociativeParameters {
    
    pub fn set_binary_from_owners(&mut self,i:usize,k:usize,epsilon:Option<f64>,beta:Option<f64>) {
        
        let mut hit = 0;
        for (j,site_j) in self.sites.iter().enumerate(){
            for (l,site_l) in self.sites.iter().enumerate(){
                hit+=1;
                if (site_j.has_owner(i)) && (site_l.has_owner(k)){

                    match epsilon {
                        Some(eps)=>{
                            self.epsilon[(j,l)] = eps;
                            self.epsilon[(l,j)] = eps
                        }
                        None => {}
                    }

                    match beta {
                        Some(bet)=>{
                            self.beta[(j,l)] = bet;
                            self.beta[(l,j)] = bet;
                        }
                        None => {}
                    }
                }
            }
        }

        if hit == 0{
            println!("no site-pair was match")
        }
    }
}

impl Parameters <AssociationPureRecord> for AssociativeParameters {
    
    fn from_records(records:Vec<AssociationPureRecord>)->Self{

        let ncomp = records.len();
        // let mut nsolv:Vec<usize>=Vec::with_capacity(ncomp);
        let mut nassoc:Vec<usize>=Vec::with_capacity(ncomp);
        for (i,record) in records.iter().enumerate(){
            // vb[i]=record.b.clone();
            match record.get_type() {
                AssociationType::Inert=>{continue;}
                AssociationType::Associative=>{nassoc.push(i)}

                // AssociationType::Solvate=>{nassoc.push(i);}
            }
        }
        let n = nassoc.len();
        let mut veps= Vec::<f64>::new();
        let mut vbeta= Vec::<f64>::new();


        let mut sites:Vec<Site> = Vec::with_capacity(NS*n); 
        let mut m:Vec<f64> = Vec::with_capacity(NS); 
        #[allow(non_snake_case)]
        let mut S = 0; //Site counter
        let mut total_sites_per_component = Vec::with_capacity(n);
        for k in 0..n{

            let i = nassoc[k]; 
            let record=&records[i];
            
            let na=record.na as f64;
            let nb=record.nb as f64;
            let nc=record.nc as f64;
            
            let mut total_sites_k = 0;

            if na!=0.0{
                sites.push(Site::new(SiteType::A, i, S,record.epsilon,record.beta));
                veps.push(record.epsilon);
                vbeta.push(record.beta);
                m.push(na);
                S+=1;
                total_sites_k+=1;
            }
            if nb!=0.0{
                sites.push(Site::new(SiteType::B, i, S,record.epsilon,record.beta));
                veps.push(record.epsilon);
                vbeta.push(record.beta);
                m.push(nb);
                S+=1;
                total_sites_k+=1;
            }
            if nc!=0.0{
                sites.push(Site::new(SiteType::C, i, S,record.epsilon,record.beta));
                veps.push(record.epsilon);
                vbeta.push(record.beta);

                m.push(nc);
                S+=1;
                total_sites_k+=1;
            }

            total_sites_per_component.push(total_sites_k);            

        }

        let veps = Array2::from_shape_vec((S,1), veps).unwrap();
        let vbeta = Array2::from_shape_vec((S,1), vbeta).unwrap();

        let multiplicity = Array1::from_vec(m);
        let ns=sites.len();
        assert_eq!(ns,S);

        let epsilon:Array2<f64> = 0.5*(&veps + &veps.t());
        let beta:Array2<f64> = (&vbeta * &vbeta.t()).sqrt();

        let mut p:Array2<f64>=Array2::zeros((S,S));


        for (j,site_j) in sites.iter().enumerate(){
            for (l,site_l) in sites.iter().enumerate(){

                if W[site_j.t()][site_l.t()] == 1.0{
                    // && (epsilon[(l,l)] != 0.0){

                    p[(j,l)]=1.0;
                    p[(l,j)]=1.0;
                }
            }
        }
        let mut lambda:Array2<f64>=Array2::zeros((n,S));
        let mut gamma:Array2<f64>=Array2::zeros((ncomp,n));
        let mut sumj = 0;
        for k in 0..n{
            
            let jk = total_sites_per_component[k];
            gamma[(nassoc[k],k)] = 1.0;

            for j in sumj..sumj+jk{
                lambda[(k,j)] = 1.0
            }
            sumj+=jk

        }
        let lambda_t=lambda.t().to_owned();
        let gamma_t=gamma.t().to_owned();

        AssociativeParameters{
            p,
            lambda,
            gamma,
            lambda_t,
            gamma_t,
            multiplicity,
            epsilon,
            beta,
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

        state.serialize_field("P", &array2_to_vec(&self.p))?;
        state.serialize_field("lambda", &array2_to_vec(&self.lambda))?;
        state.serialize_field("gamma", &array2_to_vec(&self.gamma))?;
        state.serialize_field("epsilon", &array2_to_vec(&self.epsilon))?;
        state.serialize_field("beta", &array2_to_vec(&self.beta))?;

        state.end()        
    }
}


impl std::fmt::Display for AssociativeParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        write!(
            f, 
            "sites:{:?}\nmultiplicity:{:?}\nP:\n{}\nlambda:\n{}\ngamma:\n{}\nepsilon:\n{}\nbeta:\n{}" ,
            &self.sites,
            &self.multiplicity.to_vec(),
            &array2_to_string(&self.p),
            &array2_to_string(&self.lambda),
            &array2_to_string(&self.gamma),
            &array2_to_string(&self.epsilon),
            &array2_to_string(&self.beta),)


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
impl PureRecord for AssociationPureRecord {
    
}
impl AssociationPureRecord {
    
    pub fn new(
        epsilon:f64,
        beta:f64,
        na:usize,
        nb:usize,
        nc:usize,
        // typ:AssociationType,
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
        // typ:AssociationType::Inert,
        }
    }

    pub fn solvate(sites:[usize;NS])-> Self{
        AssociationPureRecord{
        epsilon:0.0,
        beta:0.0,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        // typ:AssociationType::Associative,
        }
    }

    pub fn associative(epsilon:f64,beta:f64,sites:[usize;NS])->Self{
        AssociationPureRecord{
        epsilon,
        beta,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        // typ:AssociationType::Associative,
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