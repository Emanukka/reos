use std::f64::consts::SQRT_2;
use crate::models::cubic::parameters::CubicParameters;
use crate::residual::Residual;
use ndarray::Array1;
use ndarray::Array2;
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
        1.0 - SQRT_2
    }

    fn sig(&self)->f64 {
        1.0 + SQRT_2
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


    fn alpha_i(&self, i:usize, t:f64) -> f64 {

        let p = &self.parameters;
        let tr = t / p.tc[i];

        (1.0 + p.kappa[i] * (1.0 - tr.sqrt())).powi(2)

    }

    fn dai_dt(&self, i:usize, t:f64, alpha:f64) -> f64 {

        let p = &self.parameters;
        let a0 = p.a0[i];
        let kappa = p.kappa[i];
        let tr = t / p.tc[i];

        - a0 * kappa * (alpha * tr).sqrt() / t

    }


    fn faij(&self,t:f64)->Array2<f64>{
        
        let p = &self.parameters;
        let n = p.ncomp;
        let a0 = &p.a0;

        Array2::from_shape_fn((n,n), |(i,j)| { 

            let kij = p.aij[(i,j)] + p.bij[(i,j)] * t;
            let ai = a0[i] * self.alpha_i(i, t);
            let aj = a0[j] * self.alpha_i(j, t);

            (1. - kij) * ( ai * aj ).sqrt()
        
        })

    }

    fn da_dt(&self, t:f64, x:&Array1<f64>) -> f64 {

        let p = &self.parameters;
        let a0 = &p.a0; 

        let daij = Array2::from_shape_fn((p.ncomp,p.ncomp), |(i,j)| {

            let kij = p.aij[(i,j)] + p.bij[(i,j)] * t;
            let dkij = p.bij[(i,j)];

            let alpha_i = self.alpha_i(i, t);
            let alpha_j = self.alpha_i(j, t);

            let dai = self.dai_dt(i, t, alpha_i); 
            let daj = self.dai_dt(j, t, alpha_j); 

            let ai = a0[i] * alpha_i;
            let aj = a0[j] * alpha_j;
            let sqrt = (ai * aj).sqrt();

            0.5 * (1. - kij) * ( dai * aj + daj * ai) / sqrt - sqrt * dkij
        
        });

        x.dot(&daij.dot(x))

    }

    fn bmix(&self,x:&Array1<f64>)->f64{

        self.parameters.b.dot(x) 

    }

    fn famix(&self, x:&Array1<f64>, ax:&Array1<f64>)->f64{

        x.dot(ax)

    }

    fn fq(&self,amix:f64,bmix:f64,t:f64)->f64{
        amix / bmix / R / t
    }

    fn fdadni(&self, amix:f64, ax:&Array1<f64>)->Array1<f64>{
        
        let n = self.parameters.ncomp;

        Array1::from_shape_fn(n, |i| {

            2.0 * ax[i] - amix
            
        })


    }

    fn fi(&self, bmix:f64, vm:f64)->f64{

        let sig = self.parameters.model.sig();
        let eps = self.parameters.model.eps();

        (1.0/(sig-eps))*f64::ln((vm + sig*bmix )/(vm + eps*bmix))
    }

    fn delta(&self,bmix:f64, vm:f64) -> f64 {

        let sig = self.parameters.model.sig();
        let eps = self.parameters.model.eps();

        (vm + sig * bmix) * (vm + eps * bmix)
    }

    fn ln_z_rep(&self, bmix:f64, vm:f64) -> f64 {

        (vm / (vm - bmix) ).ln()

    }
}

impl<T:CubicModel> Residual for Cubic<T> {
    
    fn molar_weight(&self)->&Array1<f64> {
        &self.parameters.properties.molar_weight
    }
    fn components(&self)->usize {
        self.parameters.ncomp
    }

    fn max_density(&self,x:&Array1<f64>)->f64 {

        1.0 / self.bmix(x) 

    }


    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->f64 {

        let vm= 1.0 / rho;

        let aij = self.faij(t);

        let amix = self.famix(x, &aij.dot(x));
        let bmix = self.bmix(x);

        let delta = self.delta(bmix, vm);

        let repulsive= bmix / vm / (vm - bmix);
        let attractive = amix / delta / R / t;

        repulsive - attractive

    }

    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->Array1<f64> {

        let p = &self.parameters;
        let n = p.ncomp;
        let vm = 1./rho;
        let bmix= self.bmix(x);

        let aij = self.faij(t);
        let ax = aij.dot(x);
        let amix = self.famix(x, &ax);
        let dadni = &self.fdadni(amix, &ax);
        
        let ln_z_rep= self.ln_z_rep(bmix, vm);
        let q = self.fq(amix, bmix, t);
        let i_upper = self.fi(bmix, vm);

        let delta = self.delta(bmix, vm);
        let repulsive = bmix / vm / (vm - bmix);
        let attractive = amix / delta / R / t;
        let pressure = repulsive - attractive;

        let z_residual = pressure / rho;


        Array1::from_shape_fn(n, |i| {

            let db_dni = p.b[i];
            let da_dni = dadni[i];

            let dq_dni = q * (1.0 + da_dni / amix - db_dni / bmix);

            db_dni * z_residual / bmix + ln_z_rep - dq_dni * i_upper

        })



    }
    
    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>) -> f64 {

        let vm = 1./rho;
        let bmix=self.bmix(x);

        let aij = self.faij(t);
        let amix = self.famix(x, &aij.dot(x));
        
        let ln_z_rep= self.ln_z_rep(bmix, vm);
        let q = self.fq(amix, bmix, t);
        let i = self.fi(bmix, vm);       

        ln_z_rep - q * i

    }

    fn r_entropy(&self, t:f64, d:f64, x:&Array1<f64>)->f64 {

        let vm= 1./d;
        let bmix=self.bmix(x);
        let aij = self.faij(t);
        let amix = self.famix(x, &aij.dot(x));
        let q = self.fq(amix, bmix, t);
        let da_dt = self.da_dt(t, x);
        let i = self.fi(bmix, vm);       
        let ln_z_rep= self.ln_z_rep(bmix, vm);
        

        q * i * t * da_dt / amix - ln_z_rep

    }
}


mod tests{
    use approx::assert_relative_eq;
    use ndarray::array;

    use crate::{models::cubic::{Cubic, SRK, parameters::{CubicParameters, CubicPureRecord}}, parameters::{Parameters, records::PureRecord}, residual::Residual};


    fn water()->PureRecord<CubicPureRecord>{

        let m = CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14);

        PureRecord::new(0.0, "water".to_string(), m)

    }

    fn cub()->Cubic<SRK>{

        let w = water();
        let p = CubicParameters::new(vec![w], vec![]);
        Cubic::from_parameters(p)
    }

    #[test]
    fn test_srk_helmholtz() {
        
        let srk = self::cub();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = srk.r_helmholtz(t, d, &x);

        assert_relative_eq!(val, -0.058144295861,epsilon = 1e-10)

    }

    #[test]
    fn test_srk_entropy() {
        
        let srk = cub();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = srk.r_entropy(t, d, &x);

        assert_relative_eq!(val, -0.041951593945, epsilon = 1e-10)

    }

    #[test]
    fn test_srk_chem_pot() {
        
        let srk = cub();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = srk.r_chemical_potential(t, d, &x);

        assert_relative_eq!(val[0], -0.115660251059, epsilon = 1e-10)

    }
    
    #[test]
    fn test_srk_pressure() {
        
        let srk = cub();
        let t = 298.15;
        let d= 1_000.;
        let x = array![1.0];

        let val = srk.r_pressure(t, d, &x);

        assert_relative_eq!(val, -57.5159551979349, epsilon = 1e-10)

    }
    

}