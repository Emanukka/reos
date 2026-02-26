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

const TOL:f64 = 1e-12;
const RHO:f64 = 1_000.;
const T:f64 = 298.15;

#[test]
fn water_pr78_helmholtz() {
    
    let cub = water(); 
    // println!("{}", cub.parameters);
    
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
    let x = &array![1.];
    assert_relative_eq!(cub.df_dn(T, RHO, &x)[0], -0.73422735014086, epsilon = TOL)
}

#[test]
fn water_pr78_compressibility() {
    
    let cub = water(); 
    let x = &array![1.];
    assert_relative_eq!( - cub.df_dv(T, RHO, x) / RHO, -0.3634788586311659, epsilon = TOL)
}

// #[test]
// fn test_volume_translation_tcpr78() {
    
//     let cub = super::recipes::nhexane_dehlouz();
//     let param = &cub.parameters;
//     let t = 198.15;
//     let d = 8499.433742;
//     let v = 1. / d;
//     let x = array![1.];
//     let c = param.vvolt[0];
//     let a0 = param.a[0];
//     let b = param.b[0];
//     let tr = t / param.tc[0];
//     let alpha = param.alpha.alpha(0, tr);
//     let p_dehlouz = pressure_dehlouz(t, v, a0, alpha, b, c);
    
//     let p = R * t * (d + cub.df_dv(t, d, &x));
//     assert_relative_eq!(p / 1e5, p_dehlouz / 1e5, epsilon = 1e-9);
    
//     assert_relative_eq!(p / 1e5, 2.0070345678300554, epsilon = 1e-9);
    
//     // let content = format!("R = {},\nerr_P = {:.6} %", R, (p - 1e5).abs()/1e5 * 100.);
//     // let arquivo = OpenOptions::new()
//     //     .create(true)   // cria se não existir
//     //     .append(true)   // escreve no final se existir
//     //     .open("errp.txt")
//     //     .expect("Erro ao abrir/criar o arquivo");
//     // use std::io::Write;
//     // let mut writer = BufWriter::new(arquivo);
//     // writeln!(writer,"{}\n",content).expect("Erro ao escrever");
//     // assert_relative_eq!(p_dehlouz, 1e5, epsilon = 1e0);
    
// }
// #[test]
// fn crit() {
//     let cub = super::recipes::nhexane_dehlouz();
//     // let param = cub.parameters;
//     let t =  507.60;
//     // let d = 8499.433742;
//     let p = 30.25 * 1e5;
//     let eos = std::sync::Arc::new(crate::state::E::from_residual(cub));
//     let state = crate::state::State::new_tp(eos, t, p, None).unwrap();
//     println!("{}", &state);
//     // let d = state.d;
//         // let v = 1. / d;
//     // let x = array![1.];
//     // let c = param.vvolt[0];
//     // let a0 = param.a0[0];
//     // let b = param.b[0];
//     // let tr = t / param.tc[0];
//     // let alpha = param.alpha.alpha(0, tr);
    
    
// }
