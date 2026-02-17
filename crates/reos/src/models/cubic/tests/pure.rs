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
#[test]
fn water_pr78_helmholtz() {
    
    let cub = water(); 
    let t = 298.15;
    let d= 47087.634250762916;
    let x = &array![1.];
    println!("{}", cub.parameters);
    
    // let val = cub.r_helmholtz(t, d, &x);
    assert_relative_eq!(cub.r_helmholtz(t, d, x), -9.68764867952888, epsilon = 1e-10)
}
    
#[test]
fn water_pr78_entropy() {
    
    let cub = water(); 
    let t = 298.15;
    let d= 47087.634250762916;
    let x = &array![1.];
    assert_relative_eq!(cub.r_entropy(t, d, x), -7.769961493718528, epsilon = 1e-10)
}
#[test]
fn water_pr78_chem_pot() {
    
    let cub = water(); 
    let t = 298.15;
    let d= 47087.634250762916;
    let x = &array![1.];
    assert_relative_eq!(cub.r_chemical_potential(t, d, &x)[0], -10.686791988676037, epsilon = 1e-10)
}

#[test]
fn water_pr78_pressure() {
    
    let cub = water(); 
    let t = 298.15;
    let d= 47087.634250762916;
    let x = &array![1.];
    assert_relative_eq!(cub.r_pressure(t, d, x), -47047.29470521823, epsilon = 1e-10)
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
    
//     let p = R * t * (d + cub.r_pressure(t, d, &x));
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
