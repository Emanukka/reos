use core::f64;

use crate::models::cubic::mixing_rule::MixingRuleModel;
use crate::models::cubic::parameters::CubicParameters;
use crate::residual::Residual;
use ndarray::Array1;
use crate::models::R_GAS;


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


impl From<CubicParameters> for Cubic {
    fn from(value: CubicParameters) -> Self {


        Self {parameters: value}
    }
}
/// w = (B, D, d1, d2), where B = f(N), D = f(N, T) and d1,d2 are parameters of a generic cubic.
/// such that  w = u ∩ x
pub struct W {
    pub b:f64,
    pub c:f64,
    pub d:f64,
    pub d1:f64,
    pub d2:f64
}


pub struct DwDni{
    pub b:Vec<f64>,
    pub c:Vec<f64>,
    pub d:Vec<f64>
}

pub struct DwDt{
    pub d:f64,
}

struct F;

impl F {

    // fn me(t:f64, v:f64, w:&W) -> f64 {

    //     - Self::g(v, w) - w.d * Self::f(v, w) / t
    // }

    fn g(v:f64, w:&W) -> f64 {

        // let inner = (v - w.b) / v;
        let inner = (v + w.c - w.b) / v;
        inner.ln()
    }
    fn f(v:f64, w:&W) -> f64 {

        let d3 = w.d2 - w.d1;
        let r = (v + w.c + w.d2 * w.b ) / (v + w.c + w.d1 * w.b); 
        // let r = (v + c + d2 * b ) / (v + c + d1 * b); 
        
        r.ln() / R_GAS / w.b / d3
    }
}
/// Interface for derivatives of F in u = (n, T, V, B, C, D),
struct DFu ;

// Derivatives of aux. functions g(V, B) and f(V, B) 
impl DFu  {

    fn gv(v:f64, w:&W) -> f64{

        (w.b - w.c) / v / (v + w.c - w.b)
    }

    fn gb(v:f64, w:&W) -> f64{

        - 1. / (v + w.c - w.b)
    }

    fn gc(v:f64, w:&W) -> f64{

        - Self::gb(v, w)
    }

    fn fv(v:f64, w:&W) -> f64 {

        // let d3 = w.d2 - w.d1;
        - 1. / R_GAS / (v + w.c + w.d1 * w.b) / (v + w.c + w.d2 * w.b)
    }

    fn fb(v:f64, f:f64, fv:f64, w:&W) -> f64 {

        - (f + fv * (v + w.c) )/ w.b

    }

    fn fc(v:f64, w:&W) -> f64 {

        Self::fv(v, w)
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
    
    fn dc(t:f64,v:f64, w:&W) -> f64 {

        let gc = Self::gc(v, w);
        let fc = Self::fc(v, w);
        
        - gc - w.d / t * fc
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
        let df_dc = DFu::dc(t, v, w);
        
        let mut df_dni = Array1::from_elem(n, df_dn);
        
        for i in 0..n {

            df_dni[i] += {
                
                df_db * dw_dni.b[i] + 
                df_dd * dw_dni.d[i] +
                df_dc * dw_dni.c[i]
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

    fn components(&self)->usize {self.parameters.tc.len()}

    fn max_density(&self, x:&Array1<f64>)->f64 {

        let mut v_min = 0.0;
        
        for i in 0..x.len(){
            v_min += x[i] * (self.parameters.bij[(i, i)] - self.parameters.cij[(i, i)] )
        }

        0.9 / v_min

    }

    fn helmholtz(&self,t:f64, rho:f64, x:&Array1<f64>) -> f64 {

        let v = 1. / rho;
        let w = self.parameters.mix.apply(t, v, x, &self.parameters);
        // let w = self.parameters.mix.apply(t, v, x, &self.parameters);
        let g = F::g(v, &w);
        let f = F::f(v, &w);
        let d = w.d;
        
        - g - f * d / t

    }

    fn df_dt(&self, t:f64, rho:f64, x:&Array1<f64>)->f64 {


        let v = 1. / rho;
        let w = self.parameters.mix.apply(t, v, x, &self.parameters);

        let dw_dt = self.parameters.mix.dw_dt(t, v, x, &w, &self.parameters);
        let df_dt = DFx::dt(t, v, &w, &dw_dt);
        // let f = F::me(t, v, &w);

        // - t * df_dt - f
        df_dt

    }

    fn df_dv(&self, t:f64, rho:f64, x:&Array1<f64>)->f64 {

        let v = 1. / rho;
        let w = self.parameters.mix.apply(t, v, x, &self.parameters);
        let df_dv = DFx::dv(t, v, &w);

        df_dv


    }

    fn df_dn(&self,t:f64, rho:f64, x:&Array1<f64>)->Array1<f64> {
        
        let v = 1. / rho;
        let w = self.parameters.mix.apply(t, v, x, &self.parameters);
        let dw_dni = self.parameters.mix.dw_dni(t, v, x,& w, &self.parameters);

        DFx::dni(t, v, &w, &dw_dni)

    }
    

}

