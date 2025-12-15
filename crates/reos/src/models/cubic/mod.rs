use crate::models::cubic::parameters::CubicParameters;
use crate::state::eos::EosResult;
use crate::residual::Residual;
use ndarray::Array1;
use ndarray::Array2;
use ndarray::ArrayView1;
use crate::models::IDEAL_GAS_CONST as R;


pub const SRK_KAPPA_FACTORS:  [f64; 3] = [0.480000, 1.57400, -0.17600];
pub const PR76_KAPPA_FACTORS: [f64; 3] = [0.374640, 1.54226, -0.26992];
pub const PR78_KAPPA_FACTORS: [f64; 4] = [0.374642, 1.48503, -0.164423, 0.016666];

pub mod parameters;

pub trait CubicModel{
    fn model() -> Self where Self: Sized;
    fn eps(&self)->f64;
    fn sig(&self)->f64;
    fn omega(&self,)->f64;
    fn psi(&self)->f64;
    fn kappa_from_w(&self,w:f64)->f64;

    fn acrit(&self,tc:f64,pc:f64)->f64{
        self.psi()*(R*tc).powf(2.0)/pc
    }
    fn bcrit(&self,tc:f64,pc:f64)->f64{
        self.omega()*R*tc/pc
    }
    fn which(&self)-> String;
}
// #[derive(Clone)]
pub struct Cubic<T:CubicModel>{
    pub parameters:CubicParameters<T>,
}

impl<T:CubicModel> Cubic<T> {
    
    pub fn from_parameters(parameters:CubicParameters<T>)->Self{
        Self { parameters }
    }

}
pub struct SRK;
#[derive(serde::Serialize)]
pub struct PR76;
#[derive(serde::Serialize)]
pub struct PR78;


impl CubicModel for SRK{
    
    fn model()->Self where Self: Sized {
        SRK
    }
    fn omega(&self,)->f64 {
        0.08664
    }
    fn psi(&self)->f64 {
        0.42748
    }
    fn eps(&self)->f64 {
        0.0
    }
    fn sig(&self)->f64 {
        1.0
    }
    fn which(&self)-> String {
        "SRK".to_string()
    }

    fn kappa_from_w(&self,w:f64)->f64 {
        let factors = &SRK_KAPPA_FACTORS;
        factors[0] + w*factors[1] + w.powi(2)*factors[2]

    }

}


impl CubicModel for PR76 {
    
    fn model()->Self where Self: Sized {
        PR76
    }
    fn omega(&self,)->f64 {
        0.07780
    }
    fn psi(&self)->f64 {
        0.45724
    }
    fn eps(&self)->f64 {
        1.0-2.0_f64.sqrt()
    }

    fn sig(&self)->f64 {
        1.0+2.0_f64.sqrt()
    }

    fn which(&self)-> String {
        "PR76".to_string()    
    }

    fn kappa_from_w(&self,w:f64)->f64 {
        let factors = &PR76_KAPPA_FACTORS;
        factors[0] + w*factors[1] + w.powi(2)*factors[2]

    }
}

impl CubicModel for PR78 {
    
    fn model()->Self where Self: Sized {
        PR78
    }
    fn omega(&self,)->f64 {
        PR76.omega()
    }
    fn psi(&self)->f64 {
        PR76.psi()
    }
    fn eps(&self)->f64 {
        PR76.eps()
    }

    fn sig(&self)->f64 {
        PR76.sig()
    }

    fn which(&self)-> String {
        "PR78".to_string()    
    }

    fn kappa_from_w(&self,w:f64)->f64 {
        let factors76 = &PR76_KAPPA_FACTORS;
        let factors78 = &PR78_KAPPA_FACTORS;

        if w<0.491 {
            factors76[0] + w*factors76[1] + w.powi(2)*factors76[2]

        } else{
            factors78[0] + w*factors78[1] + w.powi(2)*factors78[2] + w.powi(3)*factors78[3]
        }
        
    }
}

impl <T:CubicModel> Cubic<T>{

    fn falpha(&self,
        t:f64,
        )->Array2<f64>{
        let p = &self.parameters;
        let tr =  t/&p.tc;
        let alpha =  (1.0 + &p.kappa * (1.0 - tr.sqrt())).pow2();
        alpha
    }
    
    fn faij(&self,t:f64,alpha:&Array2<f64>)->Array2<f64>{
        
        let p = &self.parameters;
        let ai = &p.a0*alpha;
        let ait = ai.t();
        let sqrt_ai_ait = (ai.dot(&ait)).sqrt();
        let kij = &p.aij + &p.bij*t; //kij = Aij + Bij*T
        (1. - kij)*sqrt_ai_ait
    }
    fn famix(&self,x:&Array1<f64>,aij:&Array2<f64>)->f64{
        x.dot(&aij.dot(x))
    }
    fn fq(&self,amix:f64,bmix:f64,t:f64)->f64{
        amix/bmix/R/t
    }
    fn fdadni(&self,amix:f64,x:&Array1<f64>,aij:&Array2<f64>)->Array1<f64>{
        2.*aij.dot(x) - amix
    }
    fn fdbdni(&self,)->ArrayView1<f64>{
        self.parameters.b.column(0)
    }
    fn fdqdni(&self,q:f64,amix:f64,bmix:f64,dadni:&Array1<f64>,dbdni:&ArrayView1<f64>)->Array1<f64>{
        q*(1. + dadni/amix - dbdni/bmix)
    }
    fn fi(&self,bmix:f64,vm:f64)->f64{
        let sig = self.parameters.model.sig();
        let eps = self.parameters.model.eps();

        (1.0/(sig-eps))*f64::ln((vm + sig*bmix )/(vm + eps*bmix))
    }
    

}

impl<T:CubicModel> Residual for Cubic<T> {
    
    fn components(&self)->usize {
        self.parameters.ncomp
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.parameters.b.t().dot(x)[0] //1xN Nx1
    }

    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        let bmix=self.bmix(x);
        let alpha = self.falpha(t);
        let aij = self.faij(t, &alpha);
        let amix=self.famix(x, &aij);
        let vm=1.0/rho;
        let delta= (vm+self.parameters.model.sig()*bmix)*(vm+self.parameters.model.eps()*bmix);

        // let repulsive= (R*t*bmix)/(vm*(vm-bmix));
        let repulsive= bmix/(vm*(vm-bmix));

        let attractive = amix/delta/R/t;

        Ok(
            repulsive-attractive
        )

    }

    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        let vm=1./rho;
        let bmix=self.bmix(x);
        let alpha = self.falpha(t);
        let aij = self.faij(t, &alpha);
        let amix = self.famix(x, &aij);
        let dadni = &self.fdadni(amix, x, &aij);
        let dbdni=&self.fdbdni();
        let ln=(vm/(vm-bmix)).ln();
        let q = self.fq(amix, bmix, t);
        let dqdni = self.fdqdni(q, amix, bmix, dadni, dbdni);

        let i = self.fi(bmix, vm);
        // let z_residual=(self.pressure(t,rho,x)?)/(rho*R*t);
        let z_residual=self.r_pressure(t,rho,x)?/rho;


        Ok(
        (dbdni/bmix)*z_residual +ln - dqdni*i
        )
    }
    
    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        let vm=1./rho;
        let bmix=self.bmix(x);
        let alpha = self.falpha(t);
        let aij = self.faij(t, &alpha);
        let amix = self.famix(x, &aij);
        let ln=(vm/(vm-bmix)).ln();
        let q= self.fq(amix, bmix, t);
        let i = self.fi(bmix, vm);        

        Ok(
            ln-q*i
        )     
    }
}


mod tests{
    use crate::{models::cubic::{Cubic, SRK, parameters::{CubicParameters, CubicPureRecord}}, parameters::Parameters, residual::Residual};


    // #[test]
    fn it_works() {
        //1:Acetic Acid, 2:Octane
        // let c1=CubicPureRecord::new(0.91196, 0.0468e-3, 0.4644, 594.8);
        // let c2=CubicPureRecord::new(34.8750e-1, 0.1424e-3, 0.99415, 568.7);

        // let mut c = Cubic::from_parameters(CubicParameters::<SRK>::from_records(vec![c1,c2]));
        // c.parameters.set_kij(0, 1,  0.064 );
        // //Set binary parameters 
        // // cpa.cubic.set_binary(0,1, Some(0.0), 0.064 );
        // let t = 353.15;
        // let alpha = c.falpha(t);
        // let faij = c.faij(t, &alpha); 
        // let x = ndarray::arr1(&[0.5,0.5]);
        // let bmix = c.bmix(&x);
        // let amix = c.famix(&x, &faij);
        // println!("amix: {}", amix);
        // println!("bmix: {}", bmix);
        // println!("fq: {}", c.fq(amix, bmix, t));

        // // let p = c.pressure(353.15, 200.0, &x).unwrap();
        // // println!("P: {}", p);       
        // // let mu = c.residual_chemical_potential(353.15, 200.0, &x).unwrap();
        // // println!("mu: {}", mu);       
        // // let a = c.residual_helmholtz(353.15, 200.0, &x).unwrap();
        // // println!("A: {}", a);       
        // //Create new State
        // // E::from_residual(cpa)
    }

}