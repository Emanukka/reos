
use ndarray::{Array1, Array2};
use crate::models::IDEAL_GAS_CONST;
use crate::models::cpa::parameters::{ ASCParameters, AssociationRule};
use crate::models::cpa::rdf::RDF;
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult};


#[derive(Clone)]
pub struct Associative<T:RDF>{
    pub parameters:ASCParameters,
    pub rdf: T
}


impl<T:RDF> Associative<T> {
    
    pub fn new(p:ASCParameters, rdf: T)->Self{
        Self{parameters:p, rdf}
        // Self{x:Array1::<f64>::default(p.f.len()),parameters:p}
    }

    pub fn set_binary(
        &mut self,
        i:usize,
        j:usize,
        rule:AssociationRule,
        eps:Option<f64>,
        beta:Option<f64>){self.parameters.set_binary_interaction(i, j, rule, eps, beta);}
            



}
impl<T:RDF> Associative<T> {
    // pub fn g_func(&self,rho:f64,x:&Array1<f64>)->f64{
    //     1.0 / (1.0 - 1.9 * (rho * self.parameters.vb.dot(x) / 4.0))
    // }

    // pub fn dlngdrho(&self,rho:f64,x:&Array1<f64>)->f64{

    //     let bm = self.parameters.vb.dot(x);
    //     let gm =self.g_func(rho,x);

    //     let result = (1.0 / gm) * (-1.0) * (gm*gm) * (-1.9 * bm / 4.0);
    //     // println!("dlngdrho={result}");
    //     result

    // }

    // // !!!
    // pub fn dlngdni(&self,rho:f64,x:&Array1<f64>)->Array1<f64>{

    //     let gmix = self.g_func(rho,x);
    //     gmix * 1.9 *(rho*&self.parameters.vb/4.0)
    // }
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

impl<T:RDF> Associative<T> {

    pub fn cr1_mat(&self,t:f64,gmix:f64)->Array2<f64>{

        let p = &self.parameters;
        let ns= p.s.len();

        let eps_mat= &p.epsmat;
        let beta_mat= &p.betamat;

        let vb=&p.vb;

        
        let b_asc=p.hmat.t().dot(vb);
        
        let b_sites=p.tmat.t().dot(&b_asc);
        let b_mat=b_sites.to_shape((ns,1)).unwrap();

        let bij_mat=(&b_mat+&b_mat.t())*0.5;

        
        bij_mat*beta_mat*gmix*((eps_mat/(IDEAL_GAS_CONST*t)).exp()-1.0)



                
    }

    pub fn ecr_mat(&self,cr1:&Array2<f64>)->Array2<f64>{

        let ns=cr1.dim().0;
        let diag=cr1.diag();
        let diag=diag.to_shape((ns,1)).unwrap();
        let diag_mult=diag.dot(&diag.t());

        (diag_mult).sqrt()
    }



    pub fn delta_mat(&self,t:f64,gmix: f64)->Array2<f64>{
        
        let alphamat=&self.parameters.alphamat;


        let cr1=self.cr1_mat(t, gmix);
        let ecr=self.ecr_mat(&cr1);

        let delta=cr1*(1.0-alphamat)+ecr*alphamat;
        delta
    }
    pub fn association_strength(&self,t:f64,rho:f64,x:&Array1<f64>)->Array2<f64>{
        let m=self.get_m(x);
        let ns=self.parameters.f.len();

        let mmat=m.to_shape((ns,1)).unwrap();
        let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
        let pmat=&self.parameters.pmat;

        // let gmix=self.g_func(rho, x);
        let gmix = self.rdf.rdf(rho, x, &self.parameters.vb);

        // let delta=self.delta(t, gmix);
        let delta=self.delta_mat(t, gmix);

        let kmat=delta*mm_mat*pmat*rho;
        kmat
    }
    #[allow(non_snake_case)]
    pub fn X_tan(&self,t:f64,rho:f64,x:&Array1<f64>)-> Result<Array1<f64>,EosError>{
        
        let ns=self.parameters.f.len();
        let s:&Array1<f64> = &self.parameters.s;
        let m=self.get_m(x);
        let kmat=self.association_strength(t, rho, x);

        let mut x_novo =  Array1::from_elem(ns,0.2);

        let omega=0.25;
        // let mut mat_error:Array2<f64> = Array2::zeros((NS,ncomp));
        let mut mat_error:Array1<f64> = Array1::zeros(ns);
        // let mut mat_error:Array2<f64> = Array2::zeros((ns,1));
        let mut x_old: Array1<f64>;
        let mut res = 1.0;
        let mut it:i32 = 0; 
        const TOL:f64 = 1e-12;
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
        // self.x=x_novo.clone();
        return Ok(x_novo);
    }
    
    }   


    pub fn grad(&self,t:f64,rho:f64,x:&Array1<f64>,x_assoc:&Array1<f64>,)->Array1<f64>{

        let m=self.get_m(x);
        let ns=self.parameters.s.len();
        let multplicity:&Array1<f64> = &self.parameters.s;
        let mmat=m.to_shape((ns,1)).unwrap();
        let mm_mat:Array2<f64>=mmat.dot(&mmat.t());
        let pmat=&self.parameters.pmat;

        let gmix = self.rdf.rdf(rho, x, &self.parameters.vb);
        let kmat=self.delta_mat(t, gmix)*mm_mat*pmat*rho;

        m*(1.0/x_assoc-1.0) -kmat.dot(&(x_assoc*multplicity))
    }


    pub fn hessian(&self,m:&Array1<f64>,k:&Array2<f64>,x_assoc:&Array1<f64>)->Array2<f64>{

        let ns=self.parameters.s.len();
        let id=Array2::<f64>::eye(ns);
        let div=m/x_assoc.pow2();
        let d=div.to_shape((ns,1)).unwrap();

        -(id.dot(&d)+k)

    }
}



impl<T: RDF> Residual for Associative<T> {

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
        let dlngdrho = self.rdf.dlngdrho(rho,x,&self.parameters.vb);

        Ok(
        -IDEAL_GAS_CONST*t*((1.0/(2.*(1./rho)))*((h)*(1.+(1./(1./rho))*dlngdrho)))

        )
    }

    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

        let xassoc=self.X_tan(t,rho, x).unwrap();
        let nassoc=&self.parameters.nassoc;

        let s=&self.parameters.s;

        let lnx_s=(xassoc.ln())*s;
        // let f=&self.parameters.f;
        let tmat=&self.parameters.tmat;
        let mu1=tmat.dot(&lnx_s); //nX1
        let hmat=&self.parameters.hmat;
        // println!("T@(ln(X)*S)=\n{}",mu1);
        let h = self.h_func(x,&xassoc);
        // let mu2 = -0.5*h*self.dlngdni(rho,x); // Vetor Ntotal
        let mu2 = -0.5*h*self.rdf.ndlngdni(rho,x,&self.parameters.vb); // Vetor Ntotal

        let mu=hmat.dot(&mu1)+mu2;
        Ok(
            mu
        )

        // Ok(
        // slogx - (h/2.0)*dlngdni
        // )
    }
    
    fn residual_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        let xassoc=self.X_tan(t,rho, x).unwrap();

        let m=self.get_m(x);

        let multiplicity=&self.parameters.s;

        let v=xassoc.ln()-0.5*xassoc+0.5;
        let w=multiplicity*m;
        Ok(
        IDEAL_GAS_CONST*t*(w.dot(&v))
        )

    }
}


