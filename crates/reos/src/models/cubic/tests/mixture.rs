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
        println!("{}", cub.parameters);
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

mod ybip {
    
    use super::*;
    use super::super::recipes::water_co2_bip;

    #[test]
    fn water_co2_pr78_helmholtz_bip() {
        
        let cub = water_co2_bip(); 
        let t = 298.15;
        let d= 41025.02831164824;
        let x = &array![0.75, 0.25];
        println!("{}", cub.parameters);
        
        // let val = cub.r_helmholtz(t, d, &x);
        assert_relative_eq!(cub.helmholtz(t, d, x), -6.330936142001641, epsilon = TOL)
    }
        
    #[test]
    fn water_co2_pr78_entropy_bip() {
        
        let cub = water_co2_bip(); 
        let t = 298.15;
        let d= 41025.02831164824;
        let x = &array![0.75, 0.25];
        // assert_relative_eq!(cub.df_dt(t, d, x), 0.04165736620095273, epsilon = TOL)
        
    }
    #[test]
    fn water_co2_pr78_chem_pot_bip() {
        
        let cub = water_co2_bip(); 
        let t = 298.15;
        let d= 41025.02831164824;
        let x = &array![0.75, 0.25];

        assert_relative_eq!(cub.df_dn(t, d, &x), array![-9.52659247284169,  1.28896444175966], epsilon = TOL)

    }

    #[test]
    fn water_co2_pr78_compressibility_bip() {
        
        let cub = water_co2_bip(); 
        let t = 298.15;
        let d= 41025.02831164824;
        let x = &array![0.75, 0.25];
        assert_relative_eq!(- cub.df_dv(t, d, x) / d, -0.491767102189709, epsilon = TOL)
    }
}



mod volume_translation {

    const T: f64 = 298.15;
    const RHO: f64 = 41025.02831164824;
    use super::{Residual, super::Cubic, super::recipes};
    use approx::assert_relative_eq;
    use ndarray::{array, Array1};

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
                
                let nc = x.len();

                let [mut x_fwd, mut x_bwd] = [x.clone(), x.clone()];
                let mut df_dx = Array1::zeros(nc);
                
                let df_dv = df_du(cub, t, d, x, Derivative::DV)[0];
                
                for i in 0..nc {
                    x_fwd[i] += EPS / 10.;
                    x_bwd[i] -= EPS / 10.;
                    
                    let df_dx_fwd = cub.helmholtz(t, d, &x_fwd);
                    let df_dx_bwd = cub.helmholtz(t, d, &x_bwd);

                    df_dx[i] =  0.5 * (df_dx_fwd - df_dx_bwd) / (EPS / 10.);
                    x_fwd[i] -= EPS / 10.;
                    x_bwd[i] += EPS / 10.;

                }

                let f = cub.helmholtz(t, d, x);
                let x_df_dx = x.dot(&df_dx);
                
                Array1::from_shape_fn(nc, |i| {

                    f + (- df_dv / d ) + df_dx[i] - x_df_dx

                })

            }
            }
        }

    #[test]
    fn valid_num_diff() {

        let cub = recipes::water_co2(); 
        let x = &array![0.75, 0.25];

        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DT)[0], cub.df_dt(T, RHO, x), max_relative = 1e-8);
        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DV)[0], cub.df_dv(T, RHO, x), max_relative = 1e-8);
        // cub.df_dn(T, RHO, x)[0];
        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DNi), cub.df_dn(T, RHO, x), max_relative = 1e-7)

    }

    #[test]
    fn test_volume_translation() {

        let cub = recipes::water_co2_vt(); 
        let x = &array![0.75, 0.25];

        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DT)[0], cub.df_dt(T, RHO, x), max_relative = 1e-8);

        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DV)[0], cub.df_dv(T, RHO, x), max_relative = 1e-8);
        
        assert_relative_eq!(df_du(&cub, T, RHO, x, Derivative::DNi), cub.df_dn(T, RHO, x), max_relative = 1e-8)
        
    }


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
    
    
