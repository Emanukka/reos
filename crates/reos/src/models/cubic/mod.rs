use core::f64;

use crate::models::cubic::mixing_rule::MixingRuleModel;
use crate::models::cubic::models::CubicModel;
use crate::models::cubic::parameters::CubicParameters;
use crate::residual::Residual;
use ndarray::Array1;
use crate::models::IDEAL_GAS_CONST as R;


pub const SRK_KAPPA_FACTORS:  [f64; 3] = [0.480000, 1.57400, -0.17600];
pub const PR76_KAPPA_FACTORS: [f64; 3] = [0.374640, 1.54226, -0.26992];
pub const PR78_KAPPA_FACTORS: [f64; 4] = [0.374642, 1.48503, -0.164423, 0.016666];

pub mod models;
pub mod parameters;
pub mod alpha;
pub mod combining_rule;
pub mod mixing_rule;
pub mod options;

#[cfg(test)]
pub mod tests;
// #[derive(Clone)]
pub struct Cubic{
    pub parameters:CubicParameters,
}

impl Cubic {
    
    pub fn from_parameters(parameters:CubicParameters) -> Self {
        Self { parameters }
    }

}


/// w = (B, D, d1, d2), where B = f(N), D = f(N, T) and d1,d2 are parameters of a generic cubic.
/// such that  w = u ∩ x
pub struct W {
    pub b:f64,
    pub d:f64,
    pub d1:f64,
    pub d2:f64
}


pub struct DwDni{
    pub b:Vec<f64>,
    pub d:Vec<f64>
}

pub struct DwDt{
    pub d:f64,
}

struct F;

impl F {

    fn me(t:f64, v:f64, w:&W) -> f64 {

        - Self::g(v, w) - w.d * Self::f(v, w) / t
    }

    fn g(v:f64, w:&W) -> f64 {

        let inner = (v - w.b) / v;
        // let inner = (v + c - b) / v;
        inner.ln()
    }
    fn f(v:f64, w:&W) -> f64 {

        let d3 = w.d2 - w.d1;
        let r = (v + w.d2 * w.b ) / (v + w.d1 * w.b); 
        // let r = (v + c + d2 * b ) / (v + c + d1 * b); 
        
        r.ln() / R / w.b / d3
    }
}
/// Interface for derivatives of F in u = (n, T, V, B, D),
struct DFu ;

// Derivatives of aux. functions g(V, B) and f(V, B) 
impl DFu  {

    fn gv(v:f64, w:&W) -> f64{

        w.b / v / (v - w.b)
    }

    fn gb(v:f64, w:&W) -> f64{

        - 1. / (v - w.b)
    }


    fn fv(v:f64, w:&W) -> f64 {

        // let d3 = w.d2 - w.d1;
        - 1. / R / (v + w.d1 * w.b) / (v + w.d2 * w.b)
    }

    fn fb(v:f64, f:f64, fv:f64, w:&W) -> f64 {

        - (f + fv * v) / w.b

    }



}

// Derivaties of F in respect to variables in u = (n, T, V, B, D)
impl DFu {

    fn dt(t:f64, f:f64, w:&W) -> f64 {
        w.d * f / t.powi(2)
    }
    
    fn db(t:f64,v:f64, w:&W) -> f64 {

        let gb = Self::gb(v, w);
        let f = F::f(v, w);
        let fv = Self::fv(v, w);
        let fb = Self::fb(v, f, fv, w);
        
        - gb - w.d / t * fb
        // let df_dd = - f / t;
    }

    fn dd(t:f64,v:f64, w:&W) -> f64 {

        let f = F::f(v, w);
        - f / t
    }

    fn dv(t:f64, v:f64, w:&W) -> f64 {

        let gv = Self::gv(v, w);
        let fv = Self::fv(v, w);

        - gv - fv * w.d / t
    }

}

/// Interface for derivaties of F in respect to x = (n, T, V)
pub struct DFx; 

impl DFx {
    fn dni(t:f64, v:f64, w:&W, dw_dni:&DwDni) -> Array1<f64> {

        let n = dw_dni.b.len();
        let df_dn = - F::g(v, w);
        let df_db = DFu::db(t, v, w);
        let df_dd = DFu::dd(t, v, w);
        
        let mut df_dni = Array1::from_elem(n, df_dn);
        
        for i in 0..n {

            df_dni[i] += {
                
                df_db * dw_dni.b[i] + 
                df_dd * dw_dni.d[i]
            }
        }

        df_dni
    }

    fn dt(t:f64, v:f64, w:&W, dw_dt:&DwDt) -> f64{

        let f = F::f(v, w);
        let df_dt = DFu::dt(t, f, w);
        let df_dd = DFu::dd(t, v, w);
        
        df_dt + df_dd * dw_dt.d

    }


    fn dv(t:f64, v:f64, w:&W,) -> f64 {

        DFu::dv(t, v, w)

    }

}
// 


impl Residual for Cubic  {
    
    fn get_properties(&self)->&crate::parameters::Properties { &self.parameters.properties }

    fn molar_weight(&self)->&Array1<f64> {&self.parameters.properties.molar_weight}

    fn components(&self)->usize {self.parameters.a.len()}

    fn max_density(&self, x:&Array1<f64>)->f64 {

        let arr = Array1::from_vec(self.parameters.b.clone());
        1.0 / (arr.dot(x))  

    }

    fn r_helmholtz(&self,t:f64, rho:f64, x:&Array1<f64>) -> f64 {

        let v = 1. / rho;
        let w = self.parameters.options.mix.apply(t, v, x, &self.parameters);
        // let w = self.parameters.options.mix.apply(t, v, x, &self.parameters);
        let g = F::g(v, &w);
        let f = F::f(v, &w);
        let d = w.d;
        
        - g - f * d / t
        // let w = self.parameters.m
        // let mix = MixParameters::new(t, d, x, &self.parameters);
        // let inner = v / (mix.tv - mix.b);
        
        // inner.ln() - mix.q * mix.i

    }

    fn r_entropy(&self, t:f64, rho:f64, x:&Array1<f64>)->f64 {

        // let mix = MixParameters::new(t, d, x, &self.parameters);
        // let mix_dt = MixTemperatureDerivatives::new(t, x, &mix, &self.parameters);
        
        // let a = self.r_helmholtz(t, d, x);
        // let da_dt = -  mix.i * mix_dt.dq_dt;
        
        // -a - t * da_dt
        let v = 1. / rho;
        let w = self.parameters.options.mix.apply(t, v, x, &self.parameters);

        let dw_dt = self.parameters.options.mix.dw_dt(t, v, x, &w, &self.parameters);
        let df_dt = DFx::dt(t, v, &w, &dw_dt);
        let f = F::me(t, v, &w);

        - t * df_dt - f
        // // let mix = MixParameters::new(t, d, x, &self.parameters);
        // let mix_dt = MixTemperatureDerivatives::new(t, x, &mix, &self.parameters);
        
        // let a = self.r_helmholtz(t, d, x);
        // let da_dt = -  mix.i * mix_dt.dq_dt;
        
        // -a - t * da_dt
    }

    fn r_pressure(&self, t:f64, rho:f64, x:&Array1<f64>)->f64 {

        let v = 1. / rho;
        let w = self.parameters.options.mix.apply(t, v, x, &self.parameters);
        let df_dv = DFx::dv(t, v, &w);

        - df_dv
        // let mix = MixParameters::new(t, d, x, &self.parameters);
        // let rep = (mix.b - mix.c) / v / (mix.tv - mix.b);
        // let att = mix.a / mix.delta / R / t;
// 
        // rep - att

    }

    fn r_chemical_potential(&self,t:f64, rho:f64, x:&Array1<f64>)->Array1<f64> {
        
        let v = 1. / rho;
        let w = self.parameters.options.mix.apply(t, v, x, &self.parameters);
        let dw_dni = self.parameters.options.mix.dw_dni(t, v, x,& w, &self.parameters);

        DFx::dni(t, v, &w, &dw_dni)

        // let n = self.parameters.a.len();
        // let mix = MixParameters::new(t, d, x, &self.parameters);
        // let mix_dn = MixMoleDerivatives::new(t, x, &mix, &self.parameters);
        // let inner = v / (mix.tv - mix.b );
        // let ln = inner.ln();
        // let z_residual = self.r_pressure(t, d, x) * v;

        // Array1::from_shape_fn(n, |i| {

        //     mix_dn.dnb_dni[i] * z_residual / mix.b + ln - mix_dn.dnq_dni[i] * mix.i

        // })


    }
    

}


// pub struct MixParameters{
//     pub a:f64,
//     pub b:f64,  
//     pub c:f64,
//     pub q:f64,
//     pub i:f64,
//     pub delta:f64,
//     pub tv:f64, //translated molar volume
// }
// impl MixParameters{

//     pub fn new(t:f64, d:f64, x:&Array1<f64>, parameters:&CubicParameters) -> Self{
        
//         let n = parameters.a.len();
//         let epsilon = parameters.options.model.eps();  
//         let sigma = parameters.options.model.sig();  
//         let bin = &parameters.binary;
//         let alpha = &parameters.options.alpha.alpha(t, &parameters.tc);
//         let [ac, bc, vt] = [&parameters.a,&parameters.b,&parameters.c];
        
//         let [mut a,mut b, mut c] = [0.,0.,0.];

//         for i in 0..n {
            
//             b += x[i] * bc[i];
//             c += x[i] * vt[i];

//             for j in 0..n {

//                 let [ai, aj] = [ac[i] * alpha[i], ac[j] * alpha[j]];
                
//                 let sqrt = (ai * aj).sqrt();
//                 let kij = bin[(i,j)].kij + bin[(i,j)].lij * t;
//                 let aij = (1. - kij) * sqrt;
//                 a += x[i] * x[j] * aij;

//             }
//         }
        
//         let tv = 1. / d + c;
//         let q = a / b / R / t;
        
//         let r = (tv + sigma * b) / (tv + epsilon * b);
//         let i = 1. / (sigma - epsilon) * r.ln();
//         let delta = (tv + sigma * b) * (tv + epsilon * b);
        
//         // dbg!(a, b, c, q, i, delta, tv);
//         Self { a, b, c, q, i, delta, tv}
//     }


// }

// pub struct MixTemperatureDerivatives{
//     pub da_dt:f64,
//     pub dq_dt:f64,
// }

// impl MixTemperatureDerivatives {

//     pub fn new(t:f64, x:&Array1<f64>, mix: &MixParameters, parameters:&CubicParameters ) -> Self{

//         let n = x.len();
//         let alpha = parameters.options.alpha.alpha(t, &parameters.tc);
//         let dalpha_dt = parameters.options.alpha.dalpha_dt(t, &parameters.tc);
//         let ac = &parameters.a;
//         let bin = &parameters.binary;
//         let [a, q] = [mix.a, mix.q];

//         let mut da_dt = 0.;

//         for i in 0..n {
//             for j in 0..n{

//                 let dkij = bin[(i,j)].lij;
//                 let kij = bin[(i,j)].kij + dkij * t;
                 
//                 let [ai, aj] = [ac[i] * alpha[i] , ac[j] * alpha[j]];
//                 let [dai, daj] = [ac[i] * dalpha_dt[i] , ac[j] * dalpha_dt[j]];

//                 let sqrt = (ai * aj).sqrt();

//                 da_dt += 0.5 * (1. - kij) * ( dai * aj + daj * ai) / sqrt - sqrt * dkij;

//             }
//         }   

//         let dq_dt = q * (da_dt / a - 1./t);

//         Self { da_dt, dq_dt }
//     }
// }
// pub struct MixMoleDerivatives{
//     pub dna_dni:Vec<f64>,
//     pub dnb_dni:Vec<f64>,
//     pub dnc_dni:Vec<f64>,
//     pub dnq_dni:Vec<f64>,
// }

// impl MixMoleDerivatives {

//     pub fn new(t:f64, x:&Array1<f64>, mix: &MixParameters, parameters:&CubicParameters ) -> Self{

//         let n = x.len();
//         let [ac, bc, vt] = [&parameters.a,&parameters.b,&parameters.c];
//         let [a, b, q] = [mix.a, mix.b, mix.q];
//         let bin = &parameters.binary;
//         let alpha = &parameters.options.alpha.alpha(t, &parameters.tc);

//         let mut dna_dni = Vec::with_capacity(n);
//         let mut dnb_dni = Vec::with_capacity(n);
//         let mut dnc_dni = Vec::with_capacity(n);
//         let mut dnq_dni = Vec::with_capacity(n);
        
//         for i in 0..n {
            
//             let dnb = bc[i];
//             let dnc = vt[i];
//             let mut sum = 0.; 
//             for j in 0..n {
                
//                 let kij = bin[(i,j)].kij + bin[(i,j)].lij * t;

//                 let [ai, aj] = [ac[i] * alpha[i], ac[j] * alpha[j]];
//                 let sqrt = (ai * aj).sqrt();

//                 let aij = (1. - kij) * sqrt;

//                 sum += aij * x[j];

//             }

//             let dna = 2. * sum - a;
//             let dnq = q * (1. + dna / a - dnb / b);
//             dna_dni.push(dna);
//             dnb_dni.push(dnb);
//             dnc_dni.push(dnc);
//             dnq_dni.push(dnq);
            
//         }

//         Self { dna_dni, dnb_dni, dnc_dni, dnq_dni }
//     }
// }
