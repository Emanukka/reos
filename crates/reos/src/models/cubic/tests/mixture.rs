use crate::residual::Residual;


use approx::assert_relative_eq;
pub use ndarray::array;

const TOL:f64 = 1e-12;

mod nbip {

    use super::*;
    use super::super::recipes::water_co2;
    use ndarray::Array1;

    const RHO:f64 = 41025.02831164824;
    const T:f64 = 298.15;
    const X:[f64;2] = [0.75, 0.25];

    #[test]
    fn water_co2_pr78_helmholtz() {
        
        let cub = water_co2(); 
        // println!("{}", cub.parameters);
        // dbg!(x);
        assert_relative_eq!(cub.helmholtz(T, RHO, &Array1::from_vec(X.to_vec())), -6.730496989050704, epsilon = TOL)

    }
        
    #[test]
    fn water_co2_pr78_df_dt() {
        
        let cub = water_co2(); 
        assert_relative_eq!(cub.df_dt(T, RHO, &Array1::from_vec(X.to_vec())), 0.04377347419896475, epsilon = TOL)
    }

    #[test]
    fn water_co2_pr78_df_dn() {
        
        let cub = water_co2(); 

        assert_relative_eq!(cub.df_dn(T, RHO, &Array1::from_vec(X.to_vec())), array![-9.98133876074576,-0.19974089691138], epsilon = TOL);

    }

    #[test]
    fn water_co2_pr78_compressibility() {
        
        let cub = water_co2(); 
        assert_relative_eq!(-cub.df_dv(T, RHO, &Array1::from_vec(X.to_vec())) / RHO, -0.8054423057364586, epsilon = TOL)
    }
}

// mod ybip {
    
//     use super::*;
//     use super::super::recipes::water_co2_bip;

//     #[test]
//     fn water_co2_pr78_helmholtz_bip() {
        
//         let cub = water_co2_bip(); 
//         let t = 298.15;
//         let d= 38.082099077791675;
//         let x = &array![0.6, 0.4];
//         println!("{}", cub.parameters);
        
//         // let val = cub.r_helmholtz(t, d, &x);
//         assert_relative_eq!(cub.r_helmholtz(t, d, x), 0.0593270749233682, epsilon = TOL)
//     }
        
//     #[test]
//     fn water_co2_pr78_entropy_bip() {
        
//         let cub = water_co2_bip(); 
//         let t = 298.15;
//         let d= 38.082099077791675;
//         let x = &array![0.6, 0.4];
//         assert_relative_eq!(cub.r_entropy(t, d, x), -0.0355454648447711, epsilon = TOL)
//     }
//     #[test]
//     fn water_co2_pr78_chem_pot_bip() {
        
//         let cub = water_co2_bip(); 
//         let t = 298.15;
//         let d= 38.082099077791675;
//         let x = &array![0.6, 0.4];

//         let lnphi = cub.r_chemical_potential(t, d, &x);
//         let reff = vec![0.091575667256, 0.159150223538];
//         for i in 0..lnphi.len() {
//             assert_relative_eq!(lnphi[i], reff[i] , epsilon = TOL)

//         }
//     }

//     #[test]
//     fn water_co2_pr78_pressure_bip() {
        
//         let cub = water_co2_bip(); 
//         let t = 298.15;
//         let d= 38.082099077791675;
//         let x = &array![0.6, 0.4];
//         assert_relative_eq!(cub.r_pressure(t, d, x), 2.257446467311823, epsilon = TOL)
//     }
// }




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
