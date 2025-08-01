
use ndarray::{Array1, Array2};
use crate::parameters::association::{ ASCParameters, AssociationRule};
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult};
use super::{IDEAL_GAS_CONST};




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
//     }a
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
    pub fn g_func(&self,rho:f64,x:&Array1<f64>)->f64{
        1.0 / (1.0 - 1.9 * (rho * self.parameters.vb.dot(x) / 4.0))
    }

    pub fn dlngdrho(&self,rho:f64,x:&Array1<f64>)->f64{

        let bm = self.parameters.vb.dot(x);
        let gm =self.g_func(rho,x);

        let result = (1.0 / gm) * (-1.0) * (gm*gm) * (-1.9 * bm / 4.0);
        // println!("dlngdrho={result}");
        result

    }

    // !!!
    pub fn dlngdni(&self,rho:f64,x:&Array1<f64>)->Array1<f64>{

        let gmix = self.g_func(rho,x);
        gmix * 1.9 *(rho*&self.parameters.vb/4.0)
    }
    pub fn get_m(&self,x:&Array1<f64>)->Array1<f64>{
        let t=&self.parameters.tmat;
        let h=&self.parameters.hmat;

        let x_asc=h.t().dot(x);
        let m=t.t().dot(&x_asc);
        m

    }
    pub fn h_func(&self,x:&Array1<f64>,xassoc:&Array1<f64>)->f64{
        
        let s = &self.parameters.s;
        let m=self.get_m(x);
        let yasc=1.0-xassoc;
        
        m.dot(&(yasc*s))

    }
}

impl Associative {
    #[allow(non_snake_case)]
    
    pub fn cr1(&self,t:f64,gmix:f64,alpha:usize,beta:usize)->f64{

        let p = &self.parameters;
        let f=&p.f;
        let map=&p.asc_into_global;
        let eps_cross= p.epsmat[(alpha,beta)];
        let beta_cross= p.betamat[(alpha,beta)];
        // let vb: ArrayBase<ndarray::OwnedRepr<&[f64]>, Dim<[usize; 2]>>= array![[vb_sclice.in]];
        // let b=(vb+vb.t())*0.5;
        let vb=&p.vb;
        let i=f[alpha].i();
        let j=f[beta].i();

        ((vb[map[i]]+vb[map[j]])*0.5)*beta_cross * gmix *(
            (( 
                ( eps_cross) / (IDEAL_GAS_CONST * t) ) ).exp() -1.0 )
    }

    pub fn ecr(&self,t:f64,gmix:f64,alpha:usize,beta:usize)->f64{
        (self.cr1(t, gmix, alpha, alpha)*self.cr1(t, gmix, beta,beta)).sqrt()
    }

    pub fn delta(&self,t:f64,gmix: f64)->Array2<f64>{
        
        let ns=self.parameters.s.len();
        let alphamat=&self.parameters.alphamat;
        let mut delta =Array2::zeros((ns,ns));
        for alpha in 0..ns{
            for beta in alpha..ns{

                let dcr1=self.cr1(t, gmix, alpha, beta);
                let decr=self.ecr(t, gmix, alpha, beta);
                let alpha_ij=alphamat[(alpha,beta)];
                delta[(alpha,beta)]=dcr1*(1.0-alpha_ij)+decr*alpha_ij;
                delta[(beta,alpha)]=delta[(alpha,beta)]

            }
        }
        delta
    }

    #[allow(non_snake_case)]
    pub fn X_tan(&self,t:f64,rho:f64,x:&Array1<f64>)-> Result<Array1<f64>,EosError>{
        let f=&self.parameters.f;
        let ns=f.len();
        let s:&Array1<f64> = &self.parameters.s;
        let m=self.get_m(x);

        let mmat=m.to_shape((ns,1)).unwrap();
        let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
        let pmat=&self.parameters.pmat;

        let gmix=self.g_func(rho, x);
        let kmat=self.delta(t, gmix)*mm_mat*pmat*rho;

        // println!("K=\n{}",&kmat);
        //nsX1
        // let mut x_novo =  Array2::from_elem((ns,1),0.2);
        let mut x_novo =  Array1::from_elem(ns,0.2);
        // let mut x_assoc =  Array2::from_elem((NS,ncomp),0.2);


        // dbg!(&x_novo);
        // dbg!(&(&x_novo*s));
        // dbg!(&(&x_n);
        // let mut x_assoc =  a.data.clone();

        let omega=0.25;
        // let mut mat_error:Array2<f64> = Array2::zeros((NS,ncomp));
        let mut mat_error:Array1<f64> = Array1::zeros(ns);
        // let mut mat_error:Array2<f64> = Array2::zeros((ns,1));
        let mut x_old: Array1<f64>;
        let mut res = 1.0;
        let mut it:i32 = 0; 
        const TOL:f64 = 1e-11; 
        const MAX:i32 = 10000;
        // let n=NS*n_assoc;

        while (res>TOL) & (it<MAX) {
            it+=1;
            x_old= x_novo.clone();
            let z=kmat.dot(&(&x_old*s));

            let correc=&m/(&m+z);
            x_novo = (1.0-omega)*&correc + omega*&x_old ;

            mat_error=((&x_novo-&x_old)/&x_novo).abs();
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

        // println!("it={}",it);
        return Ok(x_novo);
    }
    
    }   

}



impl Residual for Associative {

    fn components(&self)->usize {
        self.parameters.vb.len()
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.parameters.vb.dot(x)
    }
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        // let delta=&self.calc_delta_mat(t, rho, x);

        let xassoc=self.X_tan(t,rho, x)?;
        // dbg!("aqui");

        let h = self.h_func(x,&xassoc);
        let dlngdrho = self.dlngdrho(rho,x);

        Ok(
        -IDEAL_GAS_CONST*t*((1.0/(2.*(1./rho)))*((h)*(1.+(1./(1./rho))*dlngdrho)))

        )
    }

    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

        let xassoc=self.X_tan(t,rho, x).unwrap();
        let nassoc=&self.parameters.nassoc;

        let s=&self.parameters.s;
        let n =nassoc.len();

        let lnx_s=(xassoc.ln())*s;
        // let f=&self.parameters.f;
        let tmat=&self.parameters.tmat;
        let mu1=tmat.dot(&lnx_s); //nX1
        let hmat=&self.parameters.hmat;
        // println!("T@(ln(X)*S)=\n{}",mu1);
        let h = self.h_func(x,&xassoc);
        let mu2 = -0.5*h*self.dlngdni(rho,x); // Vetor Ntotal
        let mu=hmat.dot(&mu1)+mu2;
        Ok(
            mu
        )

        // Ok(
        // slogx - (h/2.0)*dlngdni
        // )
    }
    
}


