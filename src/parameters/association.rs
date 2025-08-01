
use core::panic;
use std::fmt::Debug;
use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use crate::models::cpa::CPA;
use crate::models::{Site,NS};
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



#[derive(Clone)]
pub struct ASCParameters{
    //Pure
    pub nassoc:   Vec<usize>,
    pub nsolv:   Vec<usize>,
    pub vb:       Array1<f64>,
    ///Indicates if compi and compj are CR1 (alpha_ij=0) or ECR (alpha_ij=1) 
    pub alphamat: Array2<f64>,
    pub epsmat: Array2<f64>,
    pub betamat: Array2<f64>,
    pub asc_into_global: Array1<usize>,
    ///f map ( site_alpha--->(type_of_Site,component) )
    pub f: Array1<Site>,
    pub s: Array1<f64>,
    /// ns x ns (allowed interactions between sites)
    pub pmat:Array2<f64>,
    /// nAssoc x ns (m=Tx'^T))
    pub tmat:Array2<f64>,
    /// N x nAssoc (x'=Hx^T)
    pub hmat:Array2<f64>,
    


}

impl ASCParameters {
    
    /// Todas interações binárias são CR1; para modificar, usar 'set_binary( )'
    pub fn from_records(records:Vec<AssociationPureRecord>)->Self{

        let ncomp = records.len();
        let mut vb = Array1::<f64>::zeros(ncomp);
        let mut nsolv:Vec<usize>=Vec::with_capacity(ncomp);
        let mut nassoc:Vec<usize>=Vec::with_capacity(ncomp);
        for (i,record) in records.iter().enumerate(){
            vb[i]=record.b.clone();
            match record.typ {
                AssociationType::Inert=>{continue;}
                AssociationType::Solvate=>{nsolv.push(i);nassoc.push(i);}
                AssociationType::Associative=>{nassoc.push(i)}
            }
        }

        // let nassoc=([nself.clone(),nsolv.clone()]).concat(); // altera a ordem
        let n =nassoc.len();

        let mut veps=Array1::<f64>::zeros(n);
        let mut vbeta=Array1::<f64>::zeros(n);

        let mut hmat:Array2<f64>=Array2::zeros((ncomp,n)); // N x nAssoc

        let mut get_global_from_associative:Array1<usize>=Array1::zeros(n); //só vou utilizar aqui; quando for realziar conta no associative.rs,
        //vou apenas fazert as operações matriciais
        
        for i in 0..ncomp{
            for (idx_assoc,&j) in nassoc.iter().enumerate(){
                if i==j{
                    get_global_from_associative[idx_assoc]=i;
                    hmat[(i,idx_assoc)]=1.0;

                    break;
                }
            }
        }

        let mut f:Vec<Site>=Vec::with_capacity(NS*n); //transformação F
        let mut s1:Vec<f64>=Vec::with_capacity(NS); //nao necessariamente vai ser preenchido, pois i pode não conter j
        
        // let mut S=Array2::<Site>::default((NS,ncomp));
        // let mut mepsilon_cross=Array2::<f64>::zeros((n,n));
        // let mut mbeta_cross=Array2::<f64>::zeros((n,n));


        for i in 0..n{
            //j : indice do componente associativo i no vetor global de tamanho Ntotal 
            //necessário pois records é len=Ntotal

            let j=get_global_from_associative[i]; 
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

            //veps tem tamanho n 
            veps[i]= record.epsilon;
            vbeta[i]=record.beta;
        }

        let f=Array1::from_vec(f);
        let ns=f.len();
        let s=Array1::from_vec(s1);
        let mut mepsilon_cross=Array2::<f64>::zeros((ns,ns));
        let mut mbeta_cross=Array2::<f64>::zeros((ns,ns));
        let alpha_mat: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 2]>>=Array2::default((ns,ns));

        let mut pmat:Array2<f64>=Array2::zeros((ns,ns));
        let mut tmat:Array2<f64>= Array2::zeros((n,ns));

        assert_eq!(f.len(),s.len());

        //1. Criando Matriz P de interçaões permitidas
        //2. Criando Matriz T de existencia do sítio ativo alpha no componente i //nsXn
        //3. Criando matriz eps e beta (ns x ns)
        // sitios->associativo->mistura
        for (alpha,site_alpha) in f.iter().enumerate(){

            let compi=site_alpha.i();
            tmat[(compi,alpha)]=1.0;
            
            for (beta,site_beta) in f.iter().enumerate(){

                let compk=site_beta.i();

                //Pega os tipos do sítio alpha e beta
                let j=site_alpha.t();
                let l=site_beta.t();

                if W[j][l]==1.0{
                    pmat[(alpha,beta)]=1.0;
                    pmat[(beta,alpha)]=1.0;

                }
                //Nao utilizar isso lá dentro, porque vai afetar no caso de solvatação
                // deixa a matriz P agir na correção

                //apesar de ser nao fisico, é muito útil para obter eps do componente
                let ecross=(veps[compi]+veps[compk])*0.5;
                let bcross=((vbeta[compi]*vbeta[compk])).sqrt();

                mepsilon_cross[(alpha,beta)]=ecross;
                mbeta_cross[(alpha,beta)]=bcross;


  
            }

        }

        ASCParameters{
            nassoc,
            nsolv,
            alphamat:alpha_mat,
            pmat,
            asc_into_global:get_global_from_associative,
            epsmat:mepsilon_cross,
            betamat:mbeta_cross,
            hmat,
            f,
            tmat,
            s,
            vb
        }

    }
    //
    pub fn set_binary_interaction(
        &mut self,
        compi:usize,
        compj:usize,
        rule:AssociationRule,
        eps:Option<f64>,
        beta:Option<f64>,
    ){
        assert_ne!(compi,compj);
        let me=&mut self.epsmat;
        let mb=&mut self.betamat;
        let alphamat=&mut self.alphamat;
        let f=&self.f;
        let nsolv=&self.nsolv;
        // dbg!(f);

        for (idx_alpha,site_alpha) in f.iter().enumerate(){
            let compi_from_site_alpha=site_alpha.i();
            if compi_from_site_alpha!=compi{continue;}

            for (idx_beta,site_beta) in f.iter().enumerate(){
                let compj_from_site_beta=site_beta.i();
                if compj_from_site_beta!=compj{continue;}       

                match rule {
                    //Apenas altera a matriz eps e beta
                    AssociationRule::MCR1=>{

                        let mut solvent_site=idx_alpha;
                        let mut solvate_site=idx_beta;

                        if nsolv.contains(&compj_from_site_beta){}
                        else if nsolv.contains(&compi_from_site_alpha){
                            //change sites
                            // println!("changing sites");
                            solvate_site= idx_alpha;
                            solvent_site= idx_beta;
                            // solvent=compj_from_site_beta;
                            // solvate=compi_from_site_alpha;
                        }else {
                            panic!("MCR1 implies 'i' is Self-Assoc., and 'j' is Solvate (v.v), but 
                            was found from pure record that they aren't. Verify the pure records.")
                        }
                        
                        me[(solvent_site,solvate_site)]=(me[(solvent_site,solvent_site)])*0.5;
                        mb[(solvent_site,solvate_site)]=beta.expect("Alert! for solv. interaction betw.'i' and 'j',it's necessary BetaCross' value.");
                        
                        me[(solvate_site,solvent_site)]=me[(solvent_site,solvate_site)];
                        mb[(solvate_site,solvent_site)]=mb[(solvent_site,solvate_site)];
                    }
                    //Apenas altera a matriz eps e beta
                    AssociationRule::EXP=>{
                        me[(idx_alpha,idx_beta)]=eps.expect("Alert! Missing  EpsCross from Experimental Data");
                        mb[(idx_alpha,idx_beta)]=beta.expect("Alert! Missing BetaCross from Experimental Data");
                        
                        me[(idx_beta,idx_alpha)]=me[(idx_alpha,idx_beta)];
                        mb[(idx_beta,idx_alpha)]=mb[(idx_alpha,idx_beta)];
                        
                    }
                    AssociationRule::ECR=>{
                        alphamat[(idx_alpha,idx_beta)]=1.0;
                        alphamat[(idx_beta,idx_alpha)]=1.0;

                    }
                    _=>{
                    }
                }
            }
        }

        // println!("alpha_ij=\n{}",&alphamat);
        // println!("eps=\n{}",&self.epsmat);
        // println!("beta=\n{}",&self.betamat);
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


impl Debug for ASCParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "  Associative Components: {:?}", self.nassoc)?;
        writeln!(f, "  Hard-sphere Volumes (vb):")?;
        writeln!(f, "    {:?}", self.vb.to_vec())?;

        writeln!(f, "F:")?;
        writeln!(f,"{:?}", self.f.to_vec())?;

        writeln!(f, "S:")?;
        writeln!(f,"{:?}", self.s.to_vec())?;
        
        writeln!(f, "T:")?;
        for row in self.tmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "H:")?;
        for row in self.hmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "P:")?;
        for row in self.pmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "E:")?;
        for row in self.epsmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "B:")?;
        for row in self.betamat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }



        Ok(())
    }
    
}
impl std::fmt::Display for ASCParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "nassoc: {:?}", self.nassoc)?;
        writeln!(f, "nsolv : {:?}", self.nsolv)?;
        writeln!(f, "b     : {}",self.vb)?;

        writeln!(f, "F: dim=({:#?}) {}",self.f.dim(),self.f)?;

        writeln!(f, "S: dim=({}) {}",self.s.dim(),self.s)?;
        
        writeln!(f, "T: dim=({}×{})",self.tmat.dim().0,self.tmat.dim().1)?;
        for row in self.tmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "H: dim=({}×{})",self.hmat.dim().0,self.hmat.dim().1)?;
        for row in self.hmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "P: dim=({}×{})",self.pmat.dim().0,self.pmat.dim().1)?;
        for row in self.pmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "E:")?;
        for row in self.epsmat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        writeln!(f, "B:")?;
        for row in self.betamat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }



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
pub fn water_octane_acetic_acid()->E<CPA>{
            //Records
        //1:Water, 2:Acetic Acid
        let c1=CubicPureRecord::new(0.12277, 0.0145e-3, 0.6736, 647.14);
        let c2=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);
        let c3=CubicPureRecord::new(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

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

        let a3=AssociationPureRecord::inert(0.1424e-3);


        //CPA eos
        let mut cpa=CPA::from_records(
            vec![c1,c3,c2],
            vec![a1,a3,a2]);

        //Set binary parameters 
        cpa.cubic.set_binary(0,1, Some(0.0),-0.222 );
        cpa.assoc.set_binary(0, 1, AssociationRule::ECR, None, None);
        
        //Create new State
        E::from_residual(cpa)

} 
pub fn acetic_acid_water()->E<CPA>{
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
            vec![c2,c1],
            vec![a2,a1]);

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
pub fn co2_water()->E<CPA>{
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
            vec![c2,c1],
            vec![a2,a1]);

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
    use std::sync::Arc;

    use approx::assert_relative_eq;
    use ndarray::{array, linalg::Dot, Array1};

    use crate::{parameters::association::{acetic_acid_water, acoh_octane, co2_water, methanol_2b, methanol_3b, octane_acoh, water_acetic_acid, water_co2, water_octane_acetic_acid}, state::{density_solver::DensityInitialization, S}};
    
    pub fn map_assoc(){
        let eos = octane_acoh();
        let h=&eos.residual.assoc.parameters.hmat;
        // println!("Assoc=\n{}",eos.residual.assoc.parameters);
        let x1=0.4;
        let x=array![x1,1.-x1];
        let x_= h.dot(&x.t());
        assert_eq!(x_[0],0.6);
    }
    // #[test]
    pub fn associative_solvate(){
        println!("---WATER & CO2---\n");
        let p=500e5;
        let t=298.15;
        let x=array![0.8,0.2];
        let eos = Arc::new(water_co2());
        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        println!{"{}",format!("{}",state)};

    }
    // #[test]
    pub fn solvate_associative(){
        println!("---CO2 & WATER---\n");
        let eos = co2_water().into();
        let p=500e5;
        let t=298.15;
        let x=array![0.2,0.8];

        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        println!{"{}",format!("{}",state)};

    }
    // #[test]
    pub fn associative_4c_associative_1a(){
        
        println!("---WATER & ACETIC ACID---\n");
        let eos = water_acetic_acid().into();
        let p=500e5;
        let t=298.15;
        let x=array![0.2,0.8];
        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        println!{"{}",format!("{}",state)};
    }
    pub fn associative_4c_inert_associative_1a(){
        
        println!("---WATER & OCTANE & ACETIC ACID---\n");
        let eos = water_octane_acetic_acid().into();
        let p=500e5;
        let t=298.15;
        let xtotal=1.+1e-20;
        let x=array![0.2/xtotal,1e-20/xtotal,0.8/xtotal];

        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();


        println!{"{}",format!("{}",state)};

        
    }
    pub fn associative_1a_associative_4c(){
        
        println!("---ACETIC ACID & WATER---\n");
        let eos = acetic_acid_water().into();
        let p=500e5;
        let t=298.15;
        let x=array![0.8,0.2];

        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        println!{"{}",format!("{}",state)};

        
    }
    #[test]
    pub fn associative_1a_inert(){
        println!("---AcOH & Octane---\n");

        let eos = acoh_octane().into();
        let t=298.15;
        let p=500e5;
        
        let x=array![0.2,0.8];
        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        println!{"{}",format!("{}",state)};

    }
    pub fn associative_3b(){

        println!("---MeOH 3B---\n");
        let eos=Arc::new(methanol_3b());
        let x=array![1.0];

        let t=298.15;
        let p=500e5;
        let cmp= [0.0007471714606619553];

        let state=S::new_tpx(&eos,
             t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        let rho=state.rho;

        let phi=state.lnphi().unwrap().exp();
        assert_relative_eq!(phi[0],cmp[0],epsilon=1e-8);

        let xassoc=eos.residual.assoc.X_tan(t, rho, &x).unwrap();
        let xa=xassoc[0];
        let xb=xassoc[1];
        // let gmix=eos.residual.assoc.g_func(rho, &x);
        assert_relative_eq!(2.*xa-1.0,xb,epsilon=1e-10);

        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();

        println!{"{}",format!("{}",state)};

    }
    pub fn associative_2b(){
        println!("---MeOH 2B---\n");

        let eos = methanol_2b().into();
        let t=298.15;
        let p=500e5;
        let x=array![1.0];

        let state=S::new_tpx(&eos, t, p, x.clone(), DensityInitialization::Vapor).unwrap();
        println!{"{}",format!("{}",state)};

    }

    #[test]
    pub fn dbg_assoc_p(){

        associative_solvate();
        println!("\n");
        solvate_associative();
        println!("\n");

        associative_4c_associative_1a();
        println!("\n");
        associative_4c_inert_associative_1a();
        println!("\n");

        associative_1a_associative_4c();
        println!("\n");

        associative_1a_inert();
        println!("\n");

        associative_2b();
        println!("\n");

        associative_3b();
        println!("\n");

        // cargo test dbg_assoc_p -- --nocapture >src/parameters/dbg/dbg_assoc_p.txt

    }
}