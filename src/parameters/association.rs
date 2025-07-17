
use core::panic;
use ndarray::{ Array1, Array2};
use serde::{Deserialize, Serialize};
use crate::models::{Site, A, B, C, NS, SITES};
use crate::state::eos::{EosError};

#[derive(PartialEq,Debug,Clone,Copy,Serialize, Deserialize)]
#[serde(rename_all="lowercase")]
pub enum AssociationRule {
    ///Eps^{ij} and Beta^{ij} from Experimental data
    EXP,
    ///Combining Rule for Solvation Interaction (Only Beta^{ij} from Experimental data )
    MCR1,
    ///Elliot's Combining Rule
    ECR,
    ///Konteogeorgis' Combining Rule
    CR1
}

impl TryFrom<String,> for AssociationRule {
    type Error =EosError;

    fn try_from(value: String) -> Result<Self, Self::Error> {
        
        match value.to_lowercase().as_str() {
            "cr1"=>{Ok(Self::CR1)},
            "ecr"=>{Ok(Self::ECR)},
            "mcr1"=>{Ok(Self::MCR1)},
            "exp"=>{Ok(Self::EXP)},
            &_=>{
                Err(
                    EosError::RuleError
                )
            }
        }
    }
}
#[derive(PartialEq,Debug,Clone,Copy,Serialize, Deserialize)]
#[serde(rename_all="lowercase")]
pub enum EpsilonRule{
    Art,
    Geo,

}
impl Default for EpsilonRule {
    fn default() -> Self {
        Self::Art
    }
}
impl EpsilonRule {
    
    fn cross(&self,i:f64,j:f64)->f64{

        match self {
            Self::Art=>{0.5*(i+j)},
            Self::Geo=>{(i*j).sqrt()}
        }
    }
}


#[derive(PartialEq,Debug,Clone,Copy)]
pub struct AssocBin{
    pub lij:f64,
    pub rule:AssociationRule,
}

impl Default for AssocBin {
    
    fn default() -> Self {
        Self::new(0.0, AssociationRule::CR1)
    }
}
impl AssocBin
{

    pub fn new(lij:f64,rule:AssociationRule)->Self{
        Self{lij,rule}
    }
}



#[derive(Clone,Debug)]
pub struct ASCParameters{
    //Pure
    pub ncomp:    usize,
    pub nself: Vec<usize>,
    pub nsolv: Vec<usize>, //redundante
    pub nassoc:Vec<usize>,
    pub vb:       Array1<f64>,
    //Binary
    pub binary: Array2<AssocBin>,
    pub site_multiplicity: Array2<Site>,
    pub eps_cross_mat: Array2<f64>,
    pub beta_cross_mat: Array2<f64>,

}
#[derive(Serialize, Deserialize,Debug,Clone,PartialEq)]
#[serde(rename_all = "lowercase")]

pub enum AssociationType{
    Inert,
    Solvate,
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
    #[serde(default)]
    pub typ:AssociationType,
    pub b:f64

}

impl AssociationPureRecord {
    
    pub fn inert(b:f64) -> Self {
        AssociationPureRecord{
        epsilon:0.0,
        beta:0.0,
        na:0,
        nb:0,
        nc:0,
        typ:AssociationType::Inert,
        b
        }
    }

    pub fn solvate(sites:[usize;NS],b:f64)-> Self{
        AssociationPureRecord{
        epsilon:0.0,
        beta:0.0,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        typ:AssociationType::Solvate,
        b
        }
    }

    pub fn associative(epsilon:f64,beta:f64,sites:[usize;NS],b:f64)->Self{
        AssociationPureRecord{
        epsilon,
        beta,
        na:sites[0],
        nb:sites[1],
        nc:sites[2],
        typ:AssociationType::Associative,
        b
        }
    }
}

// set_cr1_self-assoc-components_interactions (ou seja: pra todos componentes associativos
// será preenchida os indices i diferente de j da matrix eps e beta)
// daí dps vc poderia setar uma interação especifica pra eles 2 (mudando pra ecr; mesmo que
// tenha preenchido a matriz eps,beta, o cálculo de delta chama cr1 calculado a partir de 
// cr1)
// 
impl ASCParameters {
    
    /// Todas interações binárias são CR1; para modificar, usar 'set_binary( )'
    pub fn from_records(records:Vec<AssociationPureRecord>)->Self{

        //Pure
        let ncomp = records.len();
        let mut vb = Array1::<f64>::zeros(ncomp);
        #[allow(non_snake_case)]
        // let mut S=Array2::<f64>::zeros((NS,ncomp));
        let mut S=Array2::<Site>::default((NS,ncomp));
        //Binary (não são preenchidos agora)
        let mut mepsilon_cross=Array2::<f64>::zeros((ncomp,ncomp));
        let mut mbeta_cross=Array2::<f64>::zeros((ncomp,ncomp));
        let mut nself:Vec<usize>=Vec::with_capacity(ncomp);
        let mut nsolv:Vec<usize>=Vec::with_capacity(ncomp);
        // apesar de ter default,
        // tenho que declarar que i interage com j, mesmo que ambos sejam componenentes associativos
        // 
        let binary = Array2::<AssocBin>::default((ncomp,ncomp));
        for (i,record) in records.into_iter().enumerate(){
            vb[i]=record.b;
            match record.typ {
                AssociationType::Inert=>{continue;}
                AssociationType::Solvate=>{nsolv.push(i)}
                AssociationType::Associative=>{nself.push(i)}
                }

            S[(A,i)]=Site::A(record.na as f64);
            S[(B,i)]=Site::B(record.nb as f64);
            S[(C,i)]=Site::C(record.nc as f64);
            mepsilon_cross[(i,i)]= record.epsilon;
            mbeta_cross[(i,i)]=record.beta;

        }

        // cálculo das interações base (associativo-associativo CR1)
        for &i in &nself{
        for &j in &nself{
            if i==j{continue;}
            else{
                mepsilon_cross[(i,j)]=(mepsilon_cross[(i,i)]+mepsilon_cross[(j,j)])*0.5;
                mbeta_cross[(i,j)]=(mbeta_cross[(i,i)]*mbeta_cross[(j,j)]).sqrt();
                }
            }
        }
        let nassoc=([nself.clone(),nsolv.clone()]).concat();

        ASCParameters{
            ncomp,
            nself,
            nassoc,
            nsolv,
            binary,
            site_multiplicity:S,
            eps_cross_mat:mepsilon_cross,
            beta_cross_mat:mbeta_cross,
            vb
        }
    }
    //
    pub fn set_binary_interaction(
        &mut self,
        i:usize,
        j:usize,
        rule:AssociationRule,
        eps:Option<f64>,
        beta:Option<f64>,
        
    ){
        assert_ne!(i,j);
        let me=&mut self.eps_cross_mat;
        let mb=&mut self.beta_cross_mat;
        let binary=&mut self.binary;
        let mut k=i;
        let mut l=j;

        let nsolv=&self.nsolv;

        //k é solvente
        //l é soluto 

        match rule {
            AssociationRule::MCR1=>{
                if nsolv.contains(&j){}
                else if nsolv.contains(&i) {
                    k=j;
                    l=i;
                }else {
                    panic!("MCR1 implies 'i' is Self-Assoc., and 'j' is Solvate (v.v), but 
                    was found from pure record that they aren't. Verify the pure records.")
                }
                
                assert_ne!(me[(k,k)],0.0,"both componentes are solvates! error");
                me[(k,l)]=(me[(k,k)])*0.5;
                mb[(k,l)]=beta.expect("Alert! for solv. interaction betw.'i' and 'j',it's necessary BetaCross' value.");
            }
            AssociationRule::EXP=>{
                me[(k,l)]=eps.expect("Alert! Missing  EpsCross from Experimental Data");
                mb[(k,l)]=beta.expect("Alert! Missing BetaCross from Experimental Data");
            }
            _=>{}

        }
        
        me[(l,k)]=me[(k,l)];
        mb[(l,k)]=mb[(k,l)];

        binary[(k,l)]=AssocBin::new(0.0, rule);
        binary[(l,k)]=binary[(k,l)];

 
    }

}