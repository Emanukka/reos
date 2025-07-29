use std::cell::RefCell;

use ndarray::{Array1, Array2, Axis};
use crate::models::{A, B, C, SITES};
use crate::parameters::association::{ ASCParameters, AssociationRule, W};
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult};
use super::{IDEAL_GAS_CONST,NS};




pub struct Associative{
    pub parameters:ASCParameters
}
// pub struct Associative{
//     pub parameters:ASCParameters,
//     pub X:RefCell< X>,
// }

// struct X{
//     pub data:Array2<f64>,
//     pub rho:f64,

// }

// impl X {
    
//     fn default(rho:f64,ncomp:usize)->Self{
//         Self{
//             data:Array2::<f64>::from_elem((NS,ncomp),0.2),
//             rho
//         }
//     }
//     pub fn change_coordinate(&mut self,rho:f64,vb:&Array1<f64>,x:&Array1<f64>){
        
//         let rhomax=1.0/vb.dot(x);
//         let ncomp=x.len();
//         let s0=self.rho/rhomax;
//         let s1=rho/rhomax;
//         let eps=0.4;
//         let dif = (s1-s0).abs();
//         self.rho=rho;
//         if dif<=eps{
//         }else {
//             // println!("Mudança de coordenada");
//             self.data=Array2::<f64>::from_elem((NS,ncomp),0.2)
//         }
//     }
// }

impl Associative {
    
    // pub fn new(p:ASCParameters)->Self{
    //     let ncomp=p.ncomp;
    //     Self{parameters:p,X:
    //         RefCell::new(X::default(0.0, ncomp))}
    // }
    pub fn new(p:ASCParameters)->Self{
        Self{parameters:p}
    }

    pub fn set_binary(
        &mut self,
        i:usize,
        j:usize,
        rule:AssociationRule,
        eps:Option<f64>,
        beta:Option<f64>){self.parameters.set_binary_interaction(i, j, rule, eps, beta);}
            



}
impl Associative {
    pub fn g_func(&self,rho:f64,vx:&Array1<f64>)->f64{
        1.0 / (1.0 - 1.9 * (rho * self.parameters.vb.dot(vx) / 4.0))
    }

    pub fn dlngdrho(&self,rho:f64,vx:&Array1<f64>)->f64{

        let bm = self.parameters.vb.dot(vx);
        let gm =self.g_func(rho,vx);

        let result = (1.0 / gm) * (-1.0) * (gm*gm) * (-1.9 * bm / 4.0);
        // println!("dlngdrho={result}");
        result

    }

    // !!!
    pub fn dlngdni(&self,rho:f64,vx:&Array1<f64>)->Array1<f64>{

        let gmix = self.g_func(rho,vx);
        gmix * 1.9 *(rho*&self.parameters.vb/4.0)
    }

    pub fn h_func(&self,vx:&Array1<f64>,xassoc:&Array2<f64>)->f64{
        let s_matrix = &self.parameters.site_multiplicity;
        let mut sxx=0.0;
        
        let map=&self.parameters.map;
        let n=map.len();
        // for &i in &self.parameters.nassoc{
        for i in 0..n{
            let mut sx=0.0;
            for j in 0..NS{
                sx+=(1.0-xassoc[(j,i)])*s_matrix[(j,i)]
            }

            sxx+=vx[map[i]]*sx;
        }
        sxx
    }
}

impl Associative {
    #[allow(non_snake_case)]
    
    pub fn CR1(&self,t:f64,gmix:f64,i:usize,k:usize)->f64{
        let p = &self.parameters;
        let map=&p.map;

        let eps_cross= p.eps_cross_mat[(i,k)];
        let beta_cross= p.beta_cross_mat[(i,k)];
      
        ((p.vb[map[i]]+p.vb[map[k]])*0.5)*beta_cross * gmix *(
            (( 
                ( eps_cross) / (IDEAL_GAS_CONST * t) ) ).exp() -1.0 )
    }

    #[allow(non_snake_case)]
    pub fn ECR(&self,t:f64,gmix:f64,i:usize,k:usize)->f64
    {
        (self.CR1(t, gmix, i, i)*self.CR1(t, gmix, k,k)).sqrt()
    }

    pub fn association_constants(&self,t:f64,rho:f64,vx:&Array1<f64>,gmix: f64)->Array2<f64>{
        let nassoc=&self.parameters.nassoc;
        let n=nassoc.len();

        let f=&self.parameters.f;
        let s1=&self.parameters.s1;
        let ns=self.parameters.s1.len();
        // let n=self.parameters.ncomp;


        // let vb=&self.parameters.vb;
        // let eps=&self.parameters.eps_cross_mat;
        // let beta=&self.parameters.beta_cross_mat;

        let mut matk =Array2::zeros((ns,ns));
        // let delta=&self.parameters.delta;
        let map=&self.parameters.map;

        for alpha in 0..ns{
            for beta in 0..ns{
                if alpha>beta{continue;}
                else{
                    // println!("vi={}",vx[map[i]]);
                    let comp_i=f[alpha].i();
                    let comp_j=f[beta].i();

                    match self.parameters.binary[(comp_i,comp_j)].rule{

                    AssociationRule::CR1=>{
                        matk[(alpha,beta)]=self.CR1(t, gmix, comp_i, comp_j)
                        // matk[(i,k)]=cr1(i,k,t,gmix,vb,eps,beta)*vx[i]*vx[k]*rho
                    }
                    AssociationRule::ECR=>{
                        matk[(alpha,beta)]=self.ECR(t, gmix, comp_i, comp_j)
                        // matk[(i,k)]=ecr(i,k,t,gmix,vb,eps,beta)*vx[i]*vx[k]*rho
                    }
                    //mCR1 e exp
                    _=>{
                        matk[(alpha,beta)]=self.CR1(t, gmix, comp_i, comp_j)
                        }
                    }

                    // matk[(i,k)]=(delta[(i,k)].0)(i,k,t,gmix,vb,eps,beta)*vx[map[i]]*vx[map[k]]*rho;
                    // matk[(i,k)]=(delta[(i,k)].0)(i,k,t,gmix,vb,eps,beta)*vx[i]*vx[k]*rho;

                    matk[(beta,alpha)]=matk[(alpha,beta)]
                }
            }
        }

        matk
    }
    // pub fn association_constants(&self,t:f64,rho:f64,vx:&Array1<f64>,gmix: f64)->Array2<f64>{
    //     let nassoc=&self.parameters.nassoc;
    //     let n=nassoc.len();
    //     // let n=self.parameters.ncomp;


    //     // let vb=&self.parameters.vb;
    //     // let eps=&self.parameters.eps_cross_mat;
    //     // let beta=&self.parameters.beta_cross_mat;

    //     let mut matk =Array2::zeros((n,n));
    //     // let delta=&self.parameters.delta;
    //     let map=&self.parameters.map;
    //     for i in 0..n{
    //         for k in 0..n{
    //             if i>k{continue;}
    //             else{
    //                 // println!("vi={}",vx[map[i]]);
    //                 match self.parameters.binary[(i,k)].rule{
    //                 AssociationRule::CR1=>{
    //                     matk[(i,k)]=self.CR1(t, gmix, i, k)*vx[map[i]]*vx[map[k]]*rho
    //                     // matk[(i,k)]=cr1(i,k,t,gmix,vb,eps,beta)*vx[i]*vx[k]*rho
    //                 }
    //                 AssociationRule::ECR=>{
    //                     matk[(i,k)]=self.ECR(t, gmix, i, k)*vx[map[i]]*vx[map[k]]*rho
    //                     // matk[(i,k)]=ecr(i,k,t,gmix,vb,eps,beta)*vx[i]*vx[k]*rho
    //                 }
    //                 //mCR1 e exp
    //                 _=>{
    //                     matk[(i,k)]=self.CR1(t, gmix, i, k)*vx[map[i]]*vx[map[k]]*rho
    //                     }
    //                 }

    //                 // matk[(i,k)]=(delta[(i,k)].0)(i,k,t,gmix,vb,eps,beta)*vx[map[i]]*vx[map[k]]*rho;
    //                 // matk[(i,k)]=(delta[(i,k)].0)(i,k,t,gmix,vb,eps,beta)*vx[i]*vx[k]*rho;

    //                 matk[(k,i)]=matk[(i,k)]
    //             }
    //         }
    //     }

    //     matk
    // }

    #[allow(non_snake_case)]
    // pub fn X_tan(&self,t:f64,rho:f64,vx:&Array1<f64>)-> Result<Array2<f64>,EosError>{
    pub fn X_tan(&self,t:f64,rho:f64,vx:&Array1<f64>)-> Result<Array2<f64>,EosError>{

        let ncomp = self.parameters.ncomp;
        // let vb=&self.parameters.vb;

        // let mut a=self.X.borrow_mut();
        // a.change_coordinate(rho, vb, vx);

        let assoc_comps = &self.parameters.nassoc;
        let n=assoc_comps.len();

        let S = &self.parameters.site_multiplicity;
        let gmix=self.g_func(rho, vx);
        let matk=self.association_constants(t, rho,vx, gmix);

        let mut x_assoc =  Array2::from_elem((NS,n),0.2);
        // let mut x_assoc =  Array2::from_elem((NS,ncomp),0.2);
        let map=&self.parameters.map;
        // let mut x_assoc =  a.data.clone();
        let omega=0.25;
        // let mut mat_error:Array2<f64> = Array2::zeros((NS,ncomp));
        let mut mat_error:Array2<f64> = Array2::zeros((NS,n));
        let mut x_old: Array2<f64>;
        let mut res = 1.0;
        let mut it:i32 = 0; 
        const TOL:f64 = 1e-11; 
        const MAX:i32 = 10000;

        while (res>TOL) & (it<MAX) {

            it+=1;
            x_old= x_assoc.clone();
            // for &i in assoc_comps{
            for i in 0..n{
            for s1 in SITES {

                if S[(s1,i)]==0.0{x_assoc[(s1,i)]=1.0;continue;}
                let mut sum1: f64 = 0.0;

                for k in 0..n{
                // for &k in assoc_comps{
                for s2 in SITES{
                        sum1+= matk[(i,k)]*x_old[(s2,k)]*S[(s2,k)]*W[s1][s2]
                    }
                }
                x_assoc[(s1,i)] = (1.0-omega)*(vx[map[i]]/(vx[map[i]]+sum1)) + omega*x_assoc[(s1,i)] ;
                // x_assoc[(s1,i)] = (1.0-omega)*(vx[i]/(vx[i]+sum1)) + omega*x_assoc[(s1,i)] ;
                mat_error[(s1,i)] = ( (x_assoc[(s1,i)]-x_old[(s1,i)])/x_assoc[(s1,i)] ).abs();
                }
            }
            res = *mat_error.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
        }
        if it == MAX{return Err(EosError::NotConverged("X_tan".to_string()));}
        else {
            // x_grande.data=x_assoc.clone();
            // println!("it={it}");
            // println!("omega={omega}");
            // println!("K=\n{}",&matk);
            // println!("X={}\n",&x_assoc);
            // a.data=x_assoc.clone();

            return Ok(x_assoc);
        }
    }   

    #[allow(non_snake_case)]
    pub fn X_tan_(&self,t:f64,rho:f64,vx:&Array1<f64>)-> Result<Array2<f64>,EosError>{
        let assoc_comps = &self.parameters.nassoc;
        let n_assoc=assoc_comps.len();

        let S = &self.parameters.site_multiplicity;
        let gmix=self.g_func(rho, vx);
        let matk=self.association_constants(t, rho,vx, gmix);

        let mut x_assoc =  Array2::from_elem((NS,n_assoc),0.2);
        // let mut x_assoc =  Array2::from_elem((NS,ncomp),0.2);
        let map=&self.parameters.map;


        // let mut x_assoc =  a.data.clone();

        let omega=0.25;
        // let mut mat_error:Array2<f64> = Array2::zeros((NS,ncomp));
        let mut mat_error:Array2<f64> = Array2::zeros((NS,n_assoc));
        let mut x_old: Array2<f64>;
        let mut res = 1.0;
        let mut it:i32 = 0; 
        const TOL:f64 = 1e-11; 
        const MAX:i32 = 10000;
        let n=NS*n_assoc;

        while (res>TOL) & (it<MAX) {

            it+=1;
            x_old= x_assoc.clone();
            
            for r in 0..n{

                let i=r/NS;
                let s1=r%NS;
                let mut sum1= 0.0;

                // if S[(s1,i)]==0.0{x_assoc[(s1,i)]=1.0;continue;}
                for t in 0..n {

                    let k=t/NS;
                    let s2=t%NS;
                    // if S[(s2,k)]==0.0{continue;}
                    
                    sum1+= matk[(i,k)]*x_old[(s2,k)]*S[(s2,k)]*W[s1][s2];

                    // sum1+= delta[(i,k)]*vx[k]*x_old[(s2,k)]*S[(s2,k)].n()*(&S[(s1,i)]*&S[(s2,k)]);
                    // println!("order={}",1.0/(1.0 + rho*sum1));
                }

                x_assoc[(s1,i)] = (1.0-omega)*(vx[map[i]]/(vx[map[i]]+sum1)) + omega*x_assoc[(s1,i)] ;
                // x_assoc[(s1,i)] = (1.0-omega)*(vx[i]/(vx[i]+sum1)) + omega*x_assoc[(s1,i)] ;
                mat_error[(s1,i)] = ( (x_assoc[(s1,i)]-x_old[(s1,i)])/x_assoc[(s1,i)] ).abs();
            }
            // println!("{}  {}  {} {}",x_assoc[(0,0)], x_assoc[(1,0)], x_assoc[(2,1)], it);
            res = *mat_error.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
            //res = mat_error.sum();
            // println!("{}",res);
        }

    // dbg!(&x_assoc);
    // dbg!(res);
    // println!("materror={mat_error}");
    if it == MAX{
        return Err(EosError::NotConverged("X_tan".to_string()));
        
    }

    else {
        
        // x_grande.data=x_assoc.clone();

        // println!("X={}",&x_assoc);
        return Ok(x_assoc);
    }
    
}   

//Xmichelssen-NR
pub fn xmichelssen(){
    
}
    #[allow(non_snake_case)]
    pub fn Q(&self,t:f64,rho:f64,x:&Array1<f64>){

        let xassoc=x;


    }

}



impl Residual for Associative {

    fn components(&self)->usize {
        self.parameters.ncomp
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.parameters.vb.dot(x)
    }
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        // let delta=&self.calc_delta_mat(t, rho, x);
        let xassoc=self.X_tan(t,rho, x)?;

        let h = self.h_func(x,&xassoc);
        let dlngdrho = self.dlngdrho(rho,x);

        Ok(
        -IDEAL_GAS_CONST*t*((1.0/(2.*(1./rho)))*((h)*(1.+(1./(1./rho))*dlngdrho)))

        )
    }
    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

        // let delta=&self.calc_delta_mat(t, rho, x);
        let xassoc=self.X_tan(t,rho, x)?;
        let h = self.h_func(x,&xassoc);
        let nassoc=&self.parameters.nassoc;
        let n =nassoc.len();

        let mut mu = -0.5*h*self.dlngdni(rho,x); // Vetor Ntotal

        let map=&self.parameters.map;
        let ln_xassoc = xassoc.ln();
        let sites=&self.parameters.site_multiplicity;
        // let lnx_mult_by_smat:Array2<f64> = ln_xassoc*&self.parameters.site_multiplicity;
        // let lnx_mult_by_smat:Array2<f64> = ln_xassoc;
        // let mut mu=

        //varre tbm components inertes , pois o residual deles é não-nulo
        for i in 0..n{
            let mut slogx=0.0;
            for j in SITES{
                slogx+=ln_xassoc[(j,i)]*sites[(j,i)]
            }
            mu[map[i]]+=slogx 
        } 
        // let slogx = lnx_mult_by_smat.sum_axis(ndarray::Axis(0));
        // println!("mu={}",&mu);
        Ok(
        mu
        )

        // Ok(
        // slogx - (h/2.0)*dlngdni
        // )
    }
    
}


