
use core::panic;
use std::default;
use ndarray::{ iter, Array, Array1, Array2};
use serde::{Deserialize, Serialize};
use crate::models::cpa::CPA;
use crate::models::{Site, A, B, C, IDEAL_GAS_CONST, NS, SITES};
use crate::parameters::cubic::CubicPureRecord;
use crate::state::eos::{EosError};
use crate::state::E;
pub const W:[[f64;3];3]=[[0.0,1.0,1.0],[1.0,0.0,1.0],[1.0,1.0,1.0]];

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

// type ΔfnPtr = fn(i:usize,k:usize,t:f64,g:f64,b:&Array1<f64>,eps:&Array2<f64>,beta:&Array2<f64>)->f64; 
// #[derive(Clone,Debug)]
// pub struct Δfn(pub ΔfnPtr);

// impl Default for Δfn {
//     fn default() -> Self {
//         Δfn(cr1)
//     }
// }
// // (i,k,t,gmix,vb,eps,beta)
//         // ((p.vb[i]+p.vb[k])*0.5)*beta_cross * gmix *(
//             // (( 
//                 // ( eps_cross) / (IDEAL_GAS_CONST * t) ) ).exp() -1.0 )
// pub fn cr1(i:usize,k:usize,t:f64,g:f64,bi:f64,bk:f64,eps:&Array2<f64>,beta:&Array2<f64>)->f64{

//     // ((b[i]+b[k])*0.5)*beta[(i,k)]*g *((( eps[(i,k)] / (IDEAL_GAS_CONST * t) ) ).exp() -1.0 )
//             ((bi+bk)*0.5)*beta[(i,k)]* g *(
//             (( 
//                 ( eps[(i,k)]) / (IDEAL_GAS_CONST * t) ) ).exp() -1.0 )
// }
// pub fn ecr(i:usize,k:usize,t:f64,g:f64,b:&Array1<f64>,eps:&Array2<f64>,beta:&Array2<f64>)->f64{

//     (cr1(i, i, t, g, b, eps, beta)*cr1(k, k, t, g, b, eps, beta)).sqrt()


// }


// #[allow(non_snake_case)]
// pub fn ecr(&self,t:f64,gmix:f64,i:usize,k:usize)->f64
// {
//     (cr1(t, gmix, i, i)*self.CR1(t, gmix, k,k)).sqrt()
// }

#[derive(Clone,Debug)]
pub struct ASCParameters{
    //Pure
    pub ncomp:    usize,
    pub nsolv:    Vec<usize>, 
    pub nassoc:   Vec<usize>,
    pub vb:       Array1<f64>,
    //Binary
    pub binary: Array2<AssocBin>,
    pub eps_cross_mat: Array2<f64>,
    pub beta_cross_mat: Array2<f64>,
    pub site_multiplicity: Array2<f64>,
    pub f: Array1<Site>,
    pub s1: Array1<f64>,
    pub m:Array1<f64>,
    // pub delta:Array2<Δfn>,
    pub map:Array1<usize>,
    pub pmat:Array2<f64>
    


}

impl ASCParameters {
    
    /// Todas interações binárias são CR1; para modificar, usar 'set_binary( )'
    pub fn from_records(records:Vec<AssociationPureRecord>)->Self{

        let ncomp = records.len();
        let mut vb = Array1::<f64>::zeros(ncomp);
        let mut nself:Vec<usize>=Vec::with_capacity(ncomp);
        let mut nsolv:Vec<usize>=Vec::with_capacity(ncomp);

        //i=0,...,Ntotal

        // Ntotal=n+n'
        // n'==Nao assoc
        // n== Assoc´

        // n e n' sao vetores diferentes (tamanhos dif)
        // pra nao ter um vetor/matriz NxN (o que causa redundancia para componentes nao associatrivos)
        // cria-se um mapa (i_base_total,i_base_associativa), de modo que, quando for necessário votlar 
        // para a base total - cálculo de vetor potencial químico -, pode-se apenas consultar o mapa

        //Apenas verifica quem é associativo,solvato ou nao-associativo
        for (i,record) in records.iter().enumerate(){
            vb[i]=record.b.clone();
            match record.typ {
                AssociationType::Inert=>{continue;}
                AssociationType::Solvate=>{nsolv.push(i)}
                AssociationType::Associative=>{nself.push(i)}
            }
        }

        let nassoc=([nself.clone(),nsolv.clone()]).concat();
        let n =nassoc.len();

        //acessar map[idx i do componente associativo no espaço dos compoenntes associativos]--->index i no espaço de todos comps. da mistura

        let mut map:Array1<usize>=Array1::zeros(n);
        for i in 0..ncomp{
            for (idx_assoc,&j) in nassoc.iter().enumerate(){
                if i==j{
                    map[idx_assoc]=i;
                    break;
                }
            }
        }

        let mut f:Vec<Site>=Vec::with_capacity(NS*n); //nao necessariamente vai ser preenchido, pois i pode não conter j

        let mut s=Array2::<f64>::zeros((NS,n));
        let mut s1:Vec<f64>=Vec::with_capacity(NS); //nao necessariamente vai ser preenchido, pois i pode não conter j

        // let mut S=Array2::<Site>::default((NS,ncomp));
        let mut mepsilon_cross=Array2::<f64>::zeros((n,n));
        let mut mbeta_cross=Array2::<f64>::zeros((n,n));
        // let delta=Array2::<Δfn>::default((n,n));
        let binary = Array2::<AssocBin>::default((n,n));

        for i in 0..n{
            let j=map[i];
            let record=&records[j];

            let na=record.na as f64;
            let nb=record.nb as f64;
            let nc=record.nc as f64;
            //Não importa a ordem de f
            // mas é necessário que a ordem de f e s batam
            if na!=0.0{
                f.push(Site::A(i));
                s1.push(na);
            }
            if nb!=0.0{
                f.push(Site::B(i));
                s1.push(nb);

            }
            if nc!=0.0{
                f.push(Site::C(i));
                s1.push(nc);
            }

            s[(A,i)]=record.na as f64;
            s[(B,i)]=record.nb as f64;
            s[(C,i)]=record.nc as f64;

            mepsilon_cross[(i,i)]= record.epsilon;
            mbeta_cross[(i,i)]=record.beta;
        }
        let f=Array1::from_vec(f);
        let s1=Array1::from_vec(s1);
        let ns=f.len();
        let mut pmat:Array2<f64>=Array2::zeros((ns,ns));

        assert_eq!(f.len(),s1.len());

        for alpha in 0..f.len(){
            for beta in 0..f.len(){
                
                let j=f[alpha].t();
                let l=f[beta].t();

                if W[j][l]==1.0{
                    pmat[(alpha,beta)]=1.0;
                    pmat[(beta,alpha)]=1.0;
                }
            }
        }
        for &i in &nself{
            for &j in &nself{
                if i==j{continue;}
                else{
                    mepsilon_cross[(i,j)]=(mepsilon_cross[(i,i)]+mepsilon_cross[(j,j)])*0.5;
                    mbeta_cross[(i,j)]=(mbeta_cross[(i,i)]*mbeta_cross[(j,j)]).sqrt();
                }
            }
        }


        // println!("map={}",&map);

        ASCParameters{
            ncomp,
            nassoc,
            nsolv,
            binary,
            pmat,
            site_multiplicity:s,
            eps_cross_mat:mepsilon_cross,
            beta_cross_mat:mbeta_cross,
            map,
            f,
            s1,
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
        // let delta: &mut ndarray::ArrayBase<ndarray::OwnedRepr<Δfn>, ndarray::Dim<[usize; 2]>>=&mut self.delta;


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
            AssociationRule::ECR=>{
                // delta[(i,j)]=Δfn(ecr);
                // delta[(j,i)]=Δfn(ecr);
            }
            _=>{}

        }
        
        me[(l,k)]=me[(k,l)];
        mb[(l,k)]=mb[(k,l)];

        binary[(k,l)]=AssocBin::new(0.0, rule);
        binary[(l,k)]=binary[(k,l)];

        // println!("Delta={:#?}",delta);

 
    }

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



impl std::fmt::Display for ASCParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, " === Associative Parameters ===")?;
        writeln!(f, "  Associative Components: {:?}", self.nassoc)?;
        writeln!(f, "  Hard-sphere Volumes (vb):")?;
        writeln!(f, "    {:?}", self.vb.to_vec())?;

        writeln!(f, "  Epsilon Cross Matrix (εᵢⱼ):")?;
        for row in self.eps_cross_mat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }

        writeln!(f, "  Beta Cross Matrix (βᵢⱼ):")?;
        for row in self.beta_cross_mat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }

        writeln!(f, "  Association Site Matrix (Sⱼᵢ):")?;
        for row in self.site_multiplicity.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }


        writeln!(f, "  Binary Parameters (rule and lᵢⱼ=aT+b):")?;
        for row in self.binary.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        // writeln!(f, "  εᵢⱼ Rule): {:?}", self.)?;


        // writeln!(f, "  Solvated Components: {:?}", self.non_assoc_comps)?;

        // writeln!(f, "  Association Interactions: {:?}", self.interaction)?;
        Ok(())
    }
}



//ready-to parameters

pub fn water_acetic_acid()->E<CPA>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
            0.0145e-3
        );
        let a2=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
            0.0468e-3
        );


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c1,c2],
            vec![a1,a2]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.0),-0.222 );
        cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
        //Create new State
        E::from_residual(cpa)

} 

pub fn water_co2()->E<CPA>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new(0.35079, 0.0272e-3, 0.7602, 304.12);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [2,2,0],
            0.0145e-3
        );
        let a2=AssociationPureRecord::solvate(
            [0,1,0],
            0.0272e-3);


        //CPA eos
        
        let mut cpa=CPA::from_records(
            vec![c1,c2],
            vec![a1,a2]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.000877),-0.15508 );
        cpa.assoc.set_binary(0, 1, AssociationRule::MCR1, None, Some(0.1836));

        //Create new State
        E::from_residual(cpa)

} 
pub fn multiple_water(n:usize)->E<CPA>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);

        let a1=AssociationPureRecord::associative(
            166.55e2, 
            0.0692, 
            [1,1,1],
            0.0145e-3
        );
        let mut vec_asc: Vec<AssociationPureRecord>=vec![];
        let mut vec_cub: Vec<CubicPureRecord>=vec![];
        for _ in 0..n{
            vec_asc.push(a1.clone());
            vec_cub.push(c1.clone());
        }

        let cpa=CPA::from_records(vec_cub,vec_asc);


        //Create new State
        E::from_residual(cpa)

} 

pub fn methanol_2b()->E<CPA>{
            //Records
        //1:metoh, 2:oct
        let c1=CubicPureRecord::new(0.40531, 0.0000309, 0.4310, 513.);

        let a1=AssociationPureRecord::associative(
            24591.0, 
            0.01610, 
            [1,1,0],
            0.0000309
        );

        let  cpa=CPA::from_records(
            vec![c1],
            vec![a1]);

        
        //Create new State
        E::from_residual(cpa)
} 
pub fn methanol_3b()->E<CPA>{
            //Records
        //1:metoh, 2:oct
        let c1=
        CubicPureRecord::
        new(
            4.5897e-1, 
            0.0334e-3, 
            1.0068,
            513.);

        let a1=AssociationPureRecord::associative(
            160.70e2, 
            34.4e-3, 
            [2,1,0],
            0.0334e-3
        );

        let  cpa=CPA::from_records(
            vec![c1],
            vec![a1]);

        
        //Create new State
        E::from_residual(cpa)
} 


pub fn acoh_octane()->E<CPA>{
            //Records
        //1:Acetic Acid, 2:Octane
        let c1=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c2=CubicPureRecord::new(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

        let a1=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
            0.0468e-3
        );
        let a2=AssociationPureRecord::inert(0.1424e-3);


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c1,c2],
            vec![a1,a2]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.0), 0.064 );
        
        //Create new State
        E::from_residual(cpa)

        
} 
pub fn octane_acoh()->E<CPA>{
            //Records
        //1:Acetic Acid, 2:Octane
        let c1=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c2=CubicPureRecord::new(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

        let a1=AssociationPureRecord::associative(
            403.23e2, 
            4.5e-3, 
            [0,0,1],
            0.0468e-3
        );
        let a2=AssociationPureRecord::inert(0.1424e-3);


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c2,c1],
            vec![a2,a1]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.0), 0.064 );
        
        //Create new State
        E::from_residual(cpa)

        
} 

#[cfg(test)]
pub mod tests{
    use ndarray::array;

    use crate::parameters::association::{acoh_octane, methanol_2b, methanol_3b, octane_acoh, water_acetic_acid, water_co2};

    pub fn map_assoc(){
        
        let eos = octane_acoh();
        let map=&eos.residual.assoc.parameters.map;
        // println!("Assoc=\n{}",eos.residual.assoc.parameters);
        let x1=0.4;
        let x=array![x1,1.-x1];
        // println!("x1={}",x[map[0]]);
        assert_eq!(x[map[0]],0.6);

    }
    pub fn map_f_water_co2(){
        println!("---WATER & CO2---\n");

        let eos = water_co2();
        let map=&eos.residual.assoc.parameters.map;
        // println!("Assoc=\n{}",eos.residual.assoc.parameters);
        let x1=0.4;
        let x=array![x1,1.-x1];
        // println!("x1={}",x[map[0]]);
        println!("S={}",&eos.residual.assoc.parameters.s1);
        println!("T ( (sitio_j,comp_i)---->sitio_alpha )={}",&eos.residual.assoc.parameters.f);
        println!("P=\n{}",&eos.residual.assoc.parameters.pmat);


    }
    pub fn map_f_water_acoh(){
        
        println!("---WATER & ACETIC ACID---\n");
        let eos = water_acetic_acid();
        let map=&eos.residual.assoc.parameters.map;
        // println!("Assoc=\n{}",eos.residual.assoc.parameters);
        let x1=0.4;
        let x=array![x1,1.-x1];
        // println!("x1={}",x[map[0]]);
        println!("S={}",&eos.residual.assoc.parameters.s1);
        println!("T ( (sitio_j,comp_i)---->sitio_alpha )={}",&eos.residual.assoc.parameters.f);
        println!("P=\n{}",&eos.residual.assoc.parameters.pmat);


    }
    pub fn map_f_acoh_octane(){
        println!("---AcOH & Octane---\n");

        let eos = acoh_octane();
        let map=&eos.residual.assoc.parameters.map;
        // println!("Assoc=\n{}",eos.residual.assoc.parameters);
        let x1=0.4;
        let x=array![x1,1.-x1];
        // println!("x1={}",x[map[0]]);
        println!("S={}",&eos.residual.assoc.parameters.s1);
        println!("T ( (sitio_j,comp_i)---->sitio_alpha )={}",&eos.residual.assoc.parameters.f);
        println!("P=\n{}",&eos.residual.assoc.parameters.pmat);

    }
    pub fn map_f_metoh_3b(){

        println!("---MeOH 3B---\n");

        let eos = methanol_3b();
        println!("S={}",&eos.residual.assoc.parameters.s1);
        println!("T ( (sitio_j,comp_i)---->sitio_alpha )={}",&eos.residual.assoc.parameters.f);
        println!("P=\n{}",&eos.residual.assoc.parameters.pmat);

    }
    pub fn map_f_metoh_2b(){
                println!("---MeOH 2B---\n");

        let eos = methanol_2b();
        println!("S={}",&eos.residual.assoc.parameters.s1);
        println!("T ( (sitio_j,comp_i)---->sitio_alpha )={}",&eos.residual.assoc.parameters.f);
        println!("P=\n{}",&eos.residual.assoc.parameters.pmat);
    }

    #[test]
    pub fn map(){

        map_f_water_co2();
        println!("\n");

        map_f_water_acoh();
        println!("\n");

        map_f_acoh_octane();
        println!("\n");

        map_f_metoh_2b();
        println!("\n");

        map_f_metoh_3b();
        println!("\n");

    }
}