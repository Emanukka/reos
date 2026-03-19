use crate::residual::Residual;

use super::recipes::water;

use approx::assert_relative_eq;
use ndarray::array;
// use crate::models::IDEAL_GAS_CONST as R;
// fn nhexane_yohann() -> Cubic {
       
//     let pr = CubicPureRecord::Twu91 { 
//         tc: 507.60, 
//         pc: 30.25 * 1e5, 
//         l: 0.2557, 
//         n: 2.1871, 
//         m: 0.8377, 
//         volt: Some(0.7925 / 1e6)
//     };
//     let pr = PureRecord::new(0.0, "n-hexane".into(), pr);
//     Cubic::from_parameters(CubicParameters::new(vec![pr], vec![], PR78.into()))
// }
// fn pressure_dehlouz(t:f64, v:f64, a0:f64, alpha:f64, b:f64, c:f64,)->f64{
    
//     let rep = R * t / (v + c - b);
//     let att = a0 * alpha / ( (v + c) * (v + c + b) + b * (v + c - b) );
//     rep - att
    
// }

pub const TOL:f64 = 1e-12;
pub const RHO:f64 = 1_000.;
pub const T:f64 = 298.15;

#[test]
fn water_pr78_helmholtz() {
    
    let cub = water(); 
    println!("{}", cub.parameters);
    
    assert_relative_eq!(cub.helmholtz(T, RHO, &array![1.]), -0.37074849150969297, epsilon = TOL)
}
    
#[test]
fn water_pr78_df_dt() {
    
    let cub = water(); 
    let x = &array![1.];

    assert_relative_eq!(cub.df_dt(T, RHO, x), 0.001913855202265413, epsilon = TOL)
}
#[test]
fn water_pr78_df_dn() {
    
    let cub = water(); 
    dbg!(&cub.parameters);
    let x = &array![1.];
    assert_relative_eq!(cub.df_dn(T, RHO, &x)[0], -0.73422735014086, epsilon = TOL)
}

#[test]
fn water_pr78_compressibility() {
    
    let cub = water(); 
    let x = &array![1.];
    assert_relative_eq!( - cub.df_dv(T, RHO, x) / RHO, -0.3634788586311659, epsilon = TOL)
}

#[cfg(test)]
mod volume_translation{



    use super::{Residual, super::Cubic, super::recipes};
    use approx::assert_relative_eq;
    use ndarray::{array, Array1};
    use super::*;

    pub enum Derivative {
        DT,
        DNi,
        DV,
    }

    fn df_du(cub: &Cubic, t:f64, d:f64, x: &ndarray::Array1<f64>, du: Derivative) -> ndarray::Array1<f64> {

        const EPS:f64 = 1e-5;

        match du {
            Derivative::DT => {

                let fwd = cub.helmholtz(t + EPS, d, x);
                let bwd = cub.helmholtz(t - EPS, d, x);

                ndarray::array![0.5 * (fwd - bwd) / EPS]
            }
            Derivative::DV => {

                let fwd = cub.helmholtz(t, d + EPS * d, x);
                let bwd = cub.helmholtz(t, d - EPS * d, x);

                let df_drho = 0.5 * (fwd - bwd) / (EPS * d);
                let df_dv = df_drho * (- d.powi(2));

                ndarray::array![df_dv]
            }

            Derivative::DNi => {
                
                let fwd = cub.helmholtz(t, d + EPS * d, x);
                let bwd = cub.helmholtz(t, d - EPS * d, x);
                let df_drho = 0.5 * (fwd - bwd) / (EPS * d);

                ndarray::array![cub.helmholtz(t, d, x) + d * df_drho]
                }
            }
        }

    #[test]
    fn valid_num_diff() {

        let cub = recipes::water(); 
        let x = &array![1.0];

        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DT)[0], cub.df_dt(T, RHO, x), max_relative = 1e-8);
        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DV)[0], cub.df_dv(T, RHO, x), max_relative = 1e-8);
        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DNi)[0], cub.df_dn(T, RHO, x)[0], max_relative = 1e-8)

    }
    // left  = 363.47883036746964
    // right = 363.47885863116585

    #[test]
    fn test_volume_translation() {

        let cub = recipes::water_vt(); 
        let x = &array![1.0];

        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DT)[0], cub.df_dt(T, RHO, x), max_relative = 1e-8);

        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DV)[0], cub.df_dv(T, RHO, x), max_relative = 1e-8);
        
        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DNi)[0], cub.df_dn(T, RHO, x)[0], max_relative = 1e-8)
        
    }


}