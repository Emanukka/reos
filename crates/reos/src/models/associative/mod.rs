pub mod sites;
pub mod parameters;
use crate::models::associative::parameters::AssociativeParameters;
use crate::state::eos::{EosError};

use ndarray::{Array1, Array2};


#[derive(Clone)]
pub struct Associative{
  pub parameters:AssociativeParameters,
}


impl Associative {
    
  pub fn from_parameters(parameters:AssociativeParameters)->Self{
    Self{
      parameters
    }
  }
}

impl Associative {

    pub fn m(
        &self,
        x:&Array1<f64>
        )->Array1<f64>{
        
        let l_t:&Array2<f64> = &self.parameters.lambda_t;
        let g_t:&Array2<f64>=&self.parameters.gamma_t;
        let x_asc:Array1<f64> = g_t.dot(x);
        let m= l_t.dot(&x_asc);
        m
    }

    pub fn h(
        &self,
        x:&Array1<f64>,
        unbonded:&Array1<f64>
        )->f64{
        
        let mult = &self.parameters.multiplicity;
        let m= self.m(x);
        m.dot(&( (1.0 - unbonded) * mult ) )
    }

    pub fn association_strength(
        &self,
        t:f64,
        volf:Array2<f64>)->Array2<f64>{

        let interactions = &self.parameters.interactions;
        let s = self.parameters.sites.len();
        let mut delta = Array2::zeros((s,s));

        for interaction in interactions{
            let j = interaction.site_j.idx;
            let l = interaction.site_l.idx;
            let i = interaction.site_j.c;
            let k = interaction.site_l.c;
            let f_ii = volf[(i,i)];
            let f_kk = volf[(k,k)];

            delta[(j,l)] = interaction.association_strength_jl(t, f_ii, f_kk);
            delta[(l,j)] = delta[(j,l)]
            
        }
        delta 
        
    }

    /// $K_{jl} = \rho m_{j} m_{l} \Delta_{jl}$
    pub fn association_constants(
        &self,
        t:f64,
        rho:f64,
        x:&Array1<f64>,
        volf:Array2<f64>
        )->Array2<f64>{
        
        let interactions = &self.parameters.interactions;
        let s = self.parameters.sites.len();
        let mut kmat = Array2::zeros((s,s));
        for interaction in interactions{
            let j = interaction.site_j.idx;
            let l = interaction.site_l.idx;
            let i = interaction.site_j.c;
            let k = interaction.site_l.c;
            let f_ii = volf[(i,i)];
            let f_kk = volf[(k,k)];

            let delta_jl = interaction.association_strength_jl(t, f_ii, f_kk);
            kmat[(j,l)] = delta_jl * x[i] * x[k] * rho;
            kmat[(l,j)] = kmat[(j,l)]
            
        }
        kmat

    }

    pub fn x_tan(
      &self,
      x:&Array1<f64>,
      kmat:&Array2<f64>)-> Result<Array1<f64>,EosError>{
        
      let s = self.parameters.sites.len();
      let mult= &self.parameters.multiplicity;
      let m= self.m(x);
      let mut x_new =  Array1::from_elem(s,0.2);
      let omega= 0.25;
      let mut x_old: Array1<f64>;
      let mut res = 1.0;
      let mut it:i32 = 0; 
      let max = 500;
      while ( res > 1e-12) && ( it < max) {
          it+=1;
          x_old = x_new.clone();
          let z= kmat.dot(&(&x_old*mult));
          let correc= &m/(&m+z);
          x_new = (1.0 - omega) * &correc + omega * &x_old ;
          let err = ((&x_new-&x_old)/&x_new).abs();
          res = *err.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
      }
      
      if it == max{
          return Err(
              EosError::NotConverged("Unbonded Sites fraction".into()) )
      } else {
          return Ok(x_new) 
      }
    
    }   

    pub fn r_pressure(
      &self,
      rho:f64,
      dlng_drho:f64,
      x:&Array1<f64>,
      unbonded:&Array1<f64>
      )->f64 {

      let h = self.h(x,&unbonded);
      - rho * (1. + rho * dlng_drho) * h / 2.

    }

    pub fn r_chemical_potential(
      &self,
      x:&Array1<f64>,
      ndlng_dni:&Array1<f64>,
      unbonded:&Array1<f64>,
      )->Array1<f64>{

      let lnx= unbonded.ln() * &self.parameters.multiplicity;
      let sites= &self.parameters.sites;
      let h = self.h(x,&unbonded);
      let mut mu = - 0.5 * h * ndlng_dni;

      for site_j in sites{
          let i = site_j.c;
          let j = site_j.idx;
          mu[i] += lnx[j] 
      }

      mu
   

    }
    
    pub fn r_helmholtz(
      &self,
      x:&Array1<f64>,
      unbonded:&Array1<f64>
      )->f64{

      let m=self.m(x);
      let v= (unbonded.ln() - 0.5 * unbonded + 0.5) * &self.parameters.multiplicity;
      m.dot(&v)

    }

}


#[cfg(test)]
mod tests{
  
}
//   def Q_michelsen(self,m,X,K):

//     XM=X*self.M
//     return np.sum(m*(np.log(X)-X+1)*self.M) - 0.5*XM.T@(K@XM)

//     pub fn q_michelsen(
//         &self,
//         m:&Array1<f64>,
//         unbonded:&Array1<f64>,
//         k:&Array2<f64>) -> f64{
            
//         let mult = &self.parameters.multiplicity;
//         let xm = & (unbonded * mult);
//         (m * (unbonded.ln() - unbonded + 1.) * mult).sum() - 0.5 * xm.dot( &k.dot(xm) )

//     }
// // -self.Id*( 1/X*( m+K@(X*self.M) )) - K
//     pub fn grad(&self,
//         m:&Array1<f64>,
//         unbonded:&Array1<f64>,
//         k:&Array2<f64>) -> Array1<f64>{
        
//         let mult = &self.parameters.multiplicity;
//         let xm = & (unbonded * mult);

//         m * (1. / unbonded - 1.) - k.dot(xm)

//     }
//     pub fn hessian(
//         &self,
//         m:&Array1<f64>,
//         unbonded:&Array1<f64>,
//         k:&Array2<f64>) -> Array2<f64>{
//         let mult = &self.parameters.multiplicity;
        
//         let s = mult.len();
//         let h: DMatrix<f64> = DMatrix::zeros(s, s);
        
//         let aaaa= DVector::from(m);
//         let id = &self.parameters.id;
//         let xm = unbonded * mult;
//         - id * (
//             1. / unbonded * (
//                 m + k.dot(&xm)
//             )
//         ) - k

//     }

    // pub fn x_michelsen(
    //     &self,
    //     rho:f64,
    //     t:f64,
    //     x:&Array1<f64>) -> Array1<f64>{

    //     let s = self.parameters.sites.len();
    //     let mult= &self.parameters.multiplicity;
    //     let m= self.m(x).flatten();
    //     let k= self.k(t, rho, x);

    //     let mut x_new =  Array1::from_elem(s,0.2);

    //     let w= 0.25;
    //     let mut x_old: Array1<f64>;
    //     let mut res = 1.0;
    //     let mut it:i32 = 0; 
    //     const TOL:f64 = 1e-7;
    //     const MAX:i32 = 100;  

    //     while nsb < 5 {
            
    //         let correc = m / (k.dot(& x_old * mult));
    //         x_old = (1. - w) * correc + w * x_old;

    //     }

    //     while (res > TOL) && it < MAX {

    //         let h = self.hessian(m, &x_old, &k);
    //         let g = self.grad(m,&x_old,&k);

    //         // let hinv = 
    //         let new = LU::new(h);

    //         // x_new = solve_;
    //     } 

    //     k


    // }
    // pub fn grad(&self,t:f64,rho:f64,x:&Array1<f64>,x_assoc:&Array1<f64>,)->Array1<f64>{

    //     let m=self.m(x);
    //     let ns=self.parameters.s.len();
    //     let multplicity:&Array1<f64> = &self.parameters.s;
    //     let mmat=m.to_shape((ns,1)).unwrap();
    //     let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
    //     let pmat=&self.parameters.pmat;

    //     let gmix = self.parameters.rdf.rdf(rho, x, &self.parameters.vb);
    //     let kmat=self.delta_mat(t, gmix)*mm_mat*pmat*rho;

    //     m*(1.0/x_assoc-1.0) - kmat.dot(&(x_assoc*multplicity))
    // }


    // pub fn hessian(&self,m:&Array1<f64>,k:&Array2<f64>,x_assoc:&Array1<f64>)->Array2<f64>{

    //     let ns=self.parameters.s.len();
    //     let id=Array2::<f64>::eye(ns);
    //     let div=m/x_assoc.pow2();
    //     let d=div.to_shape((ns,1)).unwrap();

    //     -(id.dot(&d)+k)

    // }


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



// impl<R:RdfModel> Residual for Associative<R> {

//     fn components(&self)->usize {
//         self.parameters.gamma.nrows()
//     }
//     fn bmix(&self,_x:&Array1<f64>)->f64 {
//         todo!()
//     }
//     fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

//         let unbonded= self.x_tan(t,rho, x)?;
//         let h = self.h(x,&unbonded);
//         let dlngdrho = self.rdf.dlngdrho(rho, x);

//         Ok(
//         - rho * (1. + rho * dlngdrho) * h / 2.
//         )
//     }

//     fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

//         let unbonded=self.x_tan(t,rho, x).unwrap();
//         let mult=&self.parameters.multiplicity;

//         let lnx= unbonded.ln() * mult;
//         let sites=&self.parameters.sites;
//         // let lambda=&self.parameters.lambda;
//         // let mu1=lambda.dot(&lnx); //nX1
//         // let gamma=&self.parameters.gamma;
//         let h = self.h(x,&unbonded);
//         // let mu2 = -0.5 * h * self.rdf.ndlngdni(rho, x); 
//         let mut mu = -0.5 * h * self.rdf.ndlngdni(rho, x); 

//         for site_j in sites{
//             let i = site_j.c;
//             let j = site_j.j;

//             mu[i] += lnx[j] 
//         }
//         // let mu = gamma.dot(&mu1) + mu2;
//         Ok(
//             mu
//         )
//     }
    
//     fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
//         let xassoc=self.x_tan(t,rho, x).unwrap();

//         let m=self.m(x);

//         let mult=&self.parameters.multiplicity;

//         let v= (xassoc.ln() - 0.5*xassoc + 0.5)*mult;
//         Ok(
//         m.dot(&v)
//         )
//     }
// }


