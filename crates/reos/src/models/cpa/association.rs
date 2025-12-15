use std::marker::PhantomData;

use ndarray::{Array1, Array2};
use crate::models::associative::parameters::{CombiningRule, AssociativeParameters};
use crate::models::cpa::rdf::{Rdf, RdfModel};
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult};


#[derive(Clone)]
pub struct Associative<T:CombiningRule,R:RdfModel>{
    pub parameters:AssociativeParameters,
    pub rdf:Rdf<R>,
    cr:PhantomData<T>,

}


impl<T:CombiningRule,R:RdfModel> Associative<T,R> {
    
    pub fn from_parameters(
        parameters:AssociativeParameters,
        rdf:Rdf<R>
        )->Self{
        Self
        {
            parameters,
            rdf,
            cr:PhantomData
        }
    }
}
impl<T:CombiningRule,R:RdfModel> Associative<T,R> {

    pub fn m(
        &self,
        x:&Array1<f64>
        )->Array2<f64>{
        
        let x:Array2<f64> = x.to_shape((x.len(),1)).unwrap().to_owned();
        let l_t:&Array2<f64> = &self.parameters.lambda_t;
        let g_t:&Array2<f64>=&self.parameters.gamma_t;
        let x_asc:Array2<f64> = g_t.dot(&x);
        let m: Array2<f64> = l_t.dot(&x_asc);
        m
    }

    pub fn h(
        &self,
        x:&Array1<f64>,
        unbonded:&Array1<f64>
        )->f64{
        
        let mult = &self.parameters.multiplicity;
        let m= self.m(x);
        (m.t().dot(&( (1.0 - unbonded) * mult ) ))[0]
    }


    pub fn k(
        &self,
        t:f64,
        rho:f64,
        x:&Array1<f64>,
        m:&Array2<f64>,
        )->Array2<f64>{
        
        let p=&self.parameters.p;
        let g_contact = self.rdf.g(rho, x);
        let epsilon = &self.parameters.epsilon;
        let beta = &self.parameters.beta;
        let bij = &self.rdf.bij;
        let delta = T::association_strength_jl(t,g_contact,bij,epsilon,beta);
        let mm = m*&m.t();
        let kmat=delta*mm*p*rho;
        kmat
    }

    pub fn x_tan(&self,t:f64,rho:f64,x:&Array1<f64>)-> Result<Array1<f64>,EosError>{
        
        let ns = self.parameters.sites.len();
        let mult= &self.parameters.multiplicity;
        let m= self.m(x);
        let kmat=self.k(t, rho, x,&m);
        let m = m.flatten();
        let mut x_new =  Array1::from_elem(ns,0.2);

        let omega=0.25;
        let mut mat_error:Array1<f64> = Array1::zeros(ns);
        let mut x_old: Array1<f64>;
        let mut res = 1.0;
        let mut it:i32 = 0; 
        const TOL:f64 = 1e-12;
        const MAX:i32 = 10000;  
        // let n=NS*n_assoc;

        while (res>TOL) & (it<MAX) {
            it+=1;
            x_old= x_new.clone();
            let z=kmat.dot(&(&x_old*mult));

            let correc=&m/(&m+z);
            x_new = (1.0-omega)*&correc + omega*&x_old ;

            mat_error=((&x_new-&x_old)/&x_new).abs();
            res = *mat_error.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
            //res = mat_error.sum();
            // println!("{}",res);
        }
        
    if it == MAX{
        return Err(EosError::NotConverged("x_tan".to_string()));
    }
    else {
        return Ok(x_new);
    }
    
    }   


    // pub fn grad(&self,t:f64,rho:f64,x:&Array1<f64>,x_assoc:&Array1<f64>,)->Array1<f64>{

    //     let m=self.m(x);
    //     let ns=self.parameters.s.len();
    //     let multplicity:&Array1<f64> = &self.parameters.s;
    //     let mmat=m.to_shape((ns,1)).unwrap();
    //     let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
    //     let pmat=&self.parameters.pmat;

    //     let gmix = self.parameters.rdf.rdf(rho, x, &self.parameters.vb);
    //     let kmat=self.delta_mat(t, gmix)*mm_mat*pmat*rho;

    //     m*(1.0/x_assoc-1.0) -kmat.dot(&(x_assoc*multplicity))
    // }


    // pub fn hessian(&self,m:&Array1<f64>,k:&Array2<f64>,x_assoc:&Array1<f64>)->Array2<f64>{

    //     let ns=self.parameters.s.len();
    //     let id=Array2::<f64>::eye(ns);
    //     let div=m/x_assoc.pow2();
    //     let d=div.to_shape((ns,1)).unwrap();

    //     -(id.dot(&d)+k)

    // }

}

// impl<T:RDF> Associative<T> {

//     pub fn cr1_mat(&self,t:f64,gmix:f64)->Array2<f64>{

//         let p = &self.parameters;
//         let ns= p.s.len();

//         let eps_mat= &p.epsmat;
//         let beta_mat= &p.betamat;

//         let vb=&p.vb;

        
//         let b_asc=p.hmat.t().dot(vb);
        
//         let b_sites=p.tmat.t().dot(&b_asc);
//         let b_mat=b_sites.to_shape((ns,1)).unwrap();

//         let bij_mat=(&b_mat+&b_mat.t())*0.5;

        
//         bij_mat*beta_mat*gmix*((eps_mat/(IDEAL_GAS_CONST*t)).exp()-1.0)



                
//     }

//     pub fn ecr_mat(&self,cr1:&Array2<f64>)->Array2<f64>{

//         let ns=cr1.dim().0;
//         let diag=cr1.diag();
//         let diag=diag.to_shape((ns,1)).unwrap();
//         let diag_mult=diag.dot(&diag.t());

//         (diag_mult).sqrt()
//     }



//     pub fn delta_mat(&self,t:f64,gmix: f64)->Array2<f64>{
        
//         let alphamat=&self.parameters.alphamat;


//         let cr1=self.cr1_mat(t, gmix);
//         let ecr=self.ecr_mat(&cr1);

//         let delta=cr1*(1.0-alphamat)+ecr*alphamat;
//         delta
//     }
//     pub fn association_strength(&self,t:f64,rho:f64,x:&Array1<f64>)->Array2<f64>{
//         let m=self.m(x);
//         let ns=self.parameters.f.len();

//         let mmat=m.to_shape((ns,1)).unwrap();
//         let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
//         let pmat=&self.parameters.pmat;

//         // let gmix=self.g_func(rho, x);
//         let gmix = self.parameters.rdf.rdf(rho, x, &self.parameters.vb);

//         // let delta=self.delta(t, gmix);
//         let delta=self.delta_mat(t, gmix);

//         let kmat=delta*mm_mat*pmat*rho;
//         kmat
//     }
//     #[allow(non_snake_case)]
//     pub fn x_tan(&self,t:f64,rho:f64,x:&Array1<f64>)-> Result<Array1<f64>,EosError>{
        
//         let ns=self.parameters.f.len();
//         let s:&Array1<f64> = &self.parameters.s;
//         let m=self.m(x);
//         let kmat=self.association_strength(t, rho, x);

//         let mut x_new =  Array1::from_elem(ns,0.2);

//         let omega=0.25;
//         // let mut mat_error:Array2<f64> = Array2::zeros((NS,ncomp));
//         let mut mat_error:Array1<f64> = Array1::zeros(ns);
//         // let mut mat_error:Array2<f64> = Array2::zeros((ns,1));
//         let mut x_old: Array1<f64>;
//         let mut res = 1.0;
//         let mut it:i32 = 0; 
//         const TOL:f64 = 1e-12;
//         const MAX:i32 = 10000;
//         // let n=NS*n_assoc;

//         while (res>TOL) & (it<MAX) {
//             it+=1;
//             x_old= x_new.clone();
//             let z=kmat.dot(&(&x_old*s));

//             let correc=&m/(&m+z);
//             x_new = (1.0-omega)*&correc + omega*&x_old ;

//             mat_error=((&x_new-&x_old)/&x_new).abs();
//             res = *mat_error.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
//             //res = mat_error.sum();
//             // println!("{}",res);
//         }
//     // dbg!(&x_assoc);
//     // dbg!(res);
//     // println!("materror={mat_error}");
//     if it == MAX{
//         return Err(EosError::NotConverged("x_tan".to_string()));
//     }
//     else {
//         // x_grande.data=x_assoc.clone();

//         // println!("it={}",it);
//         // self.x=x_new.clone();
//         return Ok(x_new);
//     }
    
//     }   


//     pub fn grad(&self,t:f64,rho:f64,x:&Array1<f64>,x_assoc:&Array1<f64>,)->Array1<f64>{

//         let m=self.m(x);
//         let ns=self.parameters.s.len();
//         let multplicity:&Array1<f64> = &self.parameters.s;
//         let mmat=m.to_shape((ns,1)).unwrap();
//         let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
//         let pmat=&self.parameters.pmat;

//         let gmix = self.parameters.rdf.rdf(rho, x, &self.parameters.vb);
//         let kmat=self.delta_mat(t, gmix)*mm_mat*pmat*rho;

//         m*(1.0/x_assoc-1.0) -kmat.dot(&(x_assoc*multplicity))
//     }


//     pub fn hessian(&self,m:&Array1<f64>,k:&Array2<f64>,x_assoc:&Array1<f64>)->Array2<f64>{

//         let ns=self.parameters.s.len();
//         let id=Array2::<f64>::eye(ns);
//         let div=m/x_assoc.pow2();
//         let d=div.to_shape((ns,1)).unwrap();

//         -(id.dot(&d)+k)

//     }
// }



impl<T:CombiningRule,R:RdfModel> Residual for Associative<T,R> {

    fn components(&self)->usize {
        self.parameters.gamma.nrows()
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        todo!()
    }
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        // let delta=&self.calc_delta_mat(t, rho, x);

        let xassoc=self.x_tan(t,rho, x)?;
        // dbg!("aqui");

        let h = self.h(x,&xassoc);
        let dlngdrho = self.rdf.dlngdrho(rho, x);

        Ok(
        // -IDEAL_GAS_CONST*t*((1.0/(2.*(1./rho)))*((h)*(1.+(1./(1./rho))*dlngdrho)))
        -((1.0/(2.*(1./rho)))*((h)*(1.+(1./(1./rho))*dlngdrho)))


        )
    }

    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

        let unbonded=self.x_tan(t,rho, x).unwrap();
        let mult=&self.parameters.multiplicity;

        let lnx=(unbonded.ln())*mult;
        // let f=&self.parameters.f;
        let lambda=&self.parameters.lambda;
        let mu1=lambda.dot(&lnx); //nX1
        let gamma=&self.parameters.gamma;
        let h = self.h(x,&unbonded);
        let mu2 = -0.5*h*self.rdf.ndlngdni(rho, x); 

        let mu = gamma.dot(&mu1) + mu2;
        Ok(
            mu
        )
    }
    
    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        let xassoc=self.x_tan(t,rho, x).unwrap();

        let m=self.m(x);

        let mult=&self.parameters.multiplicity;

        let v= (xassoc.ln() - 0.5*xassoc + 0.5)*mult;
        Ok(
        (m.dot(&v))[0]
        )
    }
}


