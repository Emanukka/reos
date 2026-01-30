use core::f64;
use std::ops::Div;

use crate::models::cubic::parameters::CubicParameters;
use crate::residual::Residual;
use ndarray::Array1;
use ndarray::Array2;
use crate::models::IDEAL_GAS_CONST as R;


pub const SRK_KAPPA_FACTORS:  [f64; 3] = [0.480000, 1.57400, -0.17600];
pub const PR76_KAPPA_FACTORS: [f64; 3] = [0.374640, 1.54226, -0.26992];
pub const PR78_KAPPA_FACTORS: [f64; 4] = [0.374642, 1.48503, -0.164423, 0.016666];

pub mod models;
pub mod parameters;
pub mod alpha;
// #[derive(Clone)]
pub struct Cubic{
    pub parameters:CubicParameters,
}

impl Cubic {
    
    pub fn from_parameters(parameters:CubicParameters)->Self{
        Self { parameters }
    }

}



impl Cubic{


    pub fn transvol(&self, v:f64, x:&Array1<f64>) -> f64 {

        let c = self.parameters.vvolt.dot(x);
        v + c

    }
    pub fn transb(&self, b:f64, x:&Array1<f64>) -> f64 {

        let c = self.parameters.vvolt.dot(x);
        b + c

    }
    
    fn faij(&self,t:f64)->Array2<f64>{
        
        let p = &self.parameters;
        let n = p.ncomp;
        let a0 = &p.a0;
        let tc = &p.tc;
        Array2::from_shape_fn((n,n), |(i,j)| { 

            let kij = p.aij[(i,j)] + p.bij[(i,j)] * t;

            let ai = a0[i] * p.alpha.alpha(i, t / tc[i]);
            let aj = a0[j] * p.alpha.alpha(j, t / tc[j]);

            (1. - kij) * ( ai * aj ).sqrt()
        
        })

    }

    fn da_dt(&self, t:f64, x:&Array1<f64>) -> f64 {

        let p = &self.parameters;
        let a0 = &p.a0; 
        let tc = &p.tc; 

        let daij = Array2::from_shape_fn((p.ncomp,p.ncomp), |(i,j)| {

            let kij = p.aij[(i,j)] + p.bij[(i,j)] * t;
            let dkij = p.bij[(i,j)];

            let ai = a0[i] * p.alpha.alpha(i, t / tc[i]);
            let aj = a0[j] * p.alpha.alpha(j, t / tc[j]);

            let dai= a0[i] * p.alpha.dalpha_dt(i, t, tc[i]); 
            let daj= a0[j] * p.alpha.dalpha_dt(j, t, tc[j]); 

            let sqrt = (ai * aj).sqrt();

            0.5 * (1. - kij) * ( dai * aj + daj * ai) / sqrt - sqrt * dkij
        
        });

        x.dot(&daij.dot(x))

    }

    fn bmix(&self, x:&Array1<f64>)->f64{

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

    fn fi(&self, b:f64, vt:f64)->f64{

        let sig = self.parameters.sigma;
        let eps = self.parameters.epsilon;
        
        (1.0 / (sig - eps))*f64::ln((vt + sig * b )/(vt + eps * b))
    }

    fn delta(&self, bmix:f64, vm:f64) -> f64 {


        let sig = self.parameters.sigma;
        let eps = self.parameters.epsilon;

        (vm + sig * bmix) * (vm + eps * bmix)
    }

    fn tdelta(&self, b:f64, vt:f64) -> f64 {


        let sig = self.parameters.sigma;
        let eps = self.parameters.epsilon;

        (vt + sig * b) * (vt + eps * b)
    }

    fn ln_z_rep(&self, b:f64, v:f64, vt:f64) -> f64 {

        (v / (vt - b) ).ln()

    }
}

impl Residual for Cubic  {
    
    fn get_properties(&self)->&crate::parameters::Properties {
        &self.parameters.properties
    }
    fn molar_weight(&self)->&Array1<f64> {
        &self.parameters.properties.molar_weight
    }
    fn components(&self)->usize {
        self.parameters.ncomp
    }

    fn max_density(&self,x:&Array1<f64>)->f64 {

        1.0 / self.bmix(x) 

    }


    fn r_pressure(&self, t:f64, rho:f64, x:&Array1<f64>)->f64 {

        let v = 1.0 / rho;
        let b = self.bmix(x);
        let c = self.parameters.vvolt.dot(x);
        // pode transladar apenas 1?
        
        // let vtrans = self.transvol(v, x);
        // let btrans = self.transb(b, x);
        
        let aij = self.faij(t);
        let amix = self.famix(x, &aij.dot(x));

        let delta = self.tdelta(b, v + c);

        let repulsive = (b - c) / v / (v + c - b);
        let attractive = amix / delta / R / t;

        repulsive - attractive

    }

    fn r_chemical_potential(&self,t:f64, rho:f64, x:&Array1<f64>)->Array1<f64> {

        let p = &self.parameters;
        let n = p.ncomp;
        let vm = 1. / rho;
        let bmix= self.bmix(x);

        let aij = self.faij(t);
        let ax = aij.dot(x);
        let amix = self.famix(x, &ax);
        let dadni = &self.fdadni(amix, &ax);
        
        let ln_z_rep= self.ln_z_rep(bmix, vm, vm);
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

        // let vm = 1./rho;
        let vm = 1. / rho;
        let bmix = self.bmix(x);

        let aij = self.faij(t);
        let amix = self.famix(x, &aij.dot(x));
        
        let ln_z_rep= self.ln_z_rep(bmix, vm, vm);
        let q = self.fq(amix, bmix, t);
        let i = self.fi(bmix, vm);       

        ln_z_rep - q * i

    }

    fn r_entropy(&self, t:f64, d:f64, x:&Array1<f64>)->f64 {

        let v = 1.0 / d;
        let b = self.bmix(x);
        let c = self.parameters.vvolt.dot(x);
        let vt = v + c;
        
        let aij = self.faij(t);
        let amix = self.famix(x, &aij.dot(x));
        let q = self.fq(amix, b, t);
        let da_dt = self.da_dt(t, x);
        let i = self.fi(b, vt);       
        let ln_z_rep = self.ln_z_rep(b, v, vt);
        

        q * i * t * da_dt / amix - ln_z_rep

    }
}

#[cfg(test)]
pub mod utilis {
    use crate::{models::cubic::{Cubic, models::PR78, parameters::{CubicParameters, CubicPureRecord}}, parameters::{Parameters, PureRecord}};


    pub fn nhexane_dehlouz() -> Cubic {
           
        let pr = CubicPureRecord::Twu91 { 
            tc: 507.60, 
            pc: 30.25 * 1e5, 
            l: 0.28726, 
            n: 2.01991, 
            m: 0.83405, 
            volt: Some(0.81628 / 1e6)
        };
        let pr = PureRecord::new(86.17848, "n-hexane".into(), pr);
        Cubic::from_parameters(CubicParameters::new(vec![pr], vec![], PR78.into()))
    }
}

#[cfg(test)]
pub mod tests{
    
    

    use std::{fs::OpenOptions, io::BufWriter};

    use crate::{models::cubic::{Cubic, models::{PR78, SRK}, parameters::{CubicParameters, CubicPureRecord}}, parameters::{ Parameters, records::PureRecord}, residual::Residual};
    use approx::assert_relative_eq;
    use ndarray::array;
    use crate::models::IDEAL_GAS_CONST as R;
    fn water()->PureRecord<CubicPureRecord>{

        let m = CubicPureRecord::new_set1(0.12277, 0.0145e-3, 0.6736, 647.14, None);

        PureRecord::new(0.0, "water".to_string(), m)

    }


    fn cub()->Cubic{

        let w = water();
        let p = CubicParameters::new(vec![w], vec![], SRK.into());
        Cubic::from_parameters(p)
    }

    fn nhexane_yohann() -> Cubic {
           
        let pr = CubicPureRecord::Twu91 { 
            tc: 507.60, 
            pc: 30.25 * 1e5, 
            l: 0.2557, 
            n: 2.1871, 
            m: 0.8377, 
            volt: Some(0.7925 / 1e6)
        };
        let pr = PureRecord::new(0.0, "n-hexane".into(), pr);
        Cubic::from_parameters(CubicParameters::new(vec![pr], vec![], PR78.into()))
    }


    fn pressure_dehlouz(t:f64, v:f64, a0:f64, alpha:f64, b:f64, c:f64,)->f64{
        
        let rep = R * t / (v + c - b);
        let att = a0 * alpha / ( (v + c) * (v + c + b) + b * (v + c - b) );

        rep - att
        
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
    
    #[test]
    fn test_volume_translation_tcpr78() {
        
        let cub = super::utilis::nhexane_dehlouz();

        let param = &cub.parameters;

        let t = 198.15;
        let d = 8499.433742;
        let v = 1. / d;
        let x = array![1.];

        let c = param.vvolt[0];
        let a0 = param.a0[0];
        let b = param.b[0];

        let tr = t / param.tc[0];
        let alpha = param.alpha.alpha(0, tr);

        let p_dehlouz = pressure_dehlouz(t, v, a0, alpha, b, c);
        

        let p = R * t * (d + cub.r_pressure(t, d, &x));

        assert_relative_eq!(p / 1e5, p_dehlouz / 1e5, epsilon = 1e-9);
        
        assert_relative_eq!(p / 1e5, 2.0070345678300554, epsilon = 1e-9);
        
        // let content = format!("R = {},\nerr_P = {:.6} %", R, (p - 1e5).abs()/1e5 * 100.);
        // let arquivo = OpenOptions::new()
        //     .create(true)   // cria se não existir
        //     .append(true)   // escreve no final se existir
        //     .open("errp.txt")
        //     .expect("Erro ao abrir/criar o arquivo");

        // use std::io::Write;
        // let mut writer = BufWriter::new(arquivo);
        // writeln!(writer,"{}\n",content).expect("Erro ao escrever");
        // assert_relative_eq!(p_dehlouz, 1e5, epsilon = 1e0);
        
    }
}