

#[cfg(test)]
mod tests {

    use std::collections::HashMap;
    use std::fs::OpenOptions;
    use std::io::BufWriter;
    use std::sync::Arc;

    use approx::assert_relative_eq;

    use crate::models::IDEAL_GAS_CONST;
    use crate::scaling::dehlouz::{Scaling, ScalingRecord};
    use crate::{models::cubic::utilis,};
    use crate::state::{State, E};
    use crate::state::density_solver::DensityInitialization::{Liquid, Vapor};


    #[test]
    fn parameters_from_json(){

        let s = r#"
            {
                "n-hexane":[-0.174896, 0.023392, 0.569631, 0.187410, 0.512607, 0.585792]
            }
        "#;

        let map:HashMap<String, ScalingRecord> = serde_json::from_str(s).unwrap();
        let r = map.get("n-hexane").unwrap();
        
        assert_eq!(*r , [-0.174896, 0.023392, 0.569631, 0.187410, 0.512607, 0.585792]);

    }
    fn hexane_parameters() -> ScalingRecord {

        [-0.174896, 0.023392, 0.569631, 0.187410, 0.512607, 0.585792]
        
    }
    
    fn water_parameters() -> ScalingRecord {

        [0.819374, 0.159219, 0.602420, 1.278829, 0.682719, 0.681497]

    }


    #[test]
    fn liquid_hexane_dehlouz_case1() {

        let visc_exp = 0.001265;

        let parameters = hexane_parameters();

        let cub = utilis::nhexane_dehlouz();
        let eos = Arc::new(E::from_residual(cub));
        let t = 198.15;
        let p = 1e5;
        let scrit = -1.290042565;
        // let d = 8499.433742; 

        let state = Arc::new(State::new_tp(eos.clone(), t, p, Some(Liquid)).unwrap());
        let scrit = Arc::new(State::new_tp(eos.clone(), 507.60, 30.25 * 1e5, Some(Liquid)).unwrap()).entropy_isov() / IDEAL_GAS_CONST;

        // let p = state.p;
        // let d = state.d;
        
        // assert_relative_eq!(d, 8499.433742, epsilon = 1e-0); 
        
        let scaling = Scaling::new(&state, parameters);

        assert_relative_eq!(scaling.entropy, -9.019785182, epsilon = 1e-3);
        assert_relative_eq!(scaling.reference, 0.00005876, epsilon = 1e-3);

        // assert_relative_eq!(scaling.xscaling(scrit), -8.93659614 , epsilon = 1e-3);
        let visc = scaling.viscosity(scrit);
        // assert_relative_eq!(visc, 0.00125461 , epsilon = 1e-3);
        let err = format!("err_exp = {:.6} %", (visc - visc_exp).abs() / visc_exp * 100.);
        dbg!(visc);

        println!("{}",err);
        
        let content = format!("R = {},\nerr_ρ = {:.6} %,\nerr_s = {:.6} %,\nerr_η = {:.6} %,\n{}",
        IDEAL_GAS_CONST,
        (state.d - 8499.433742).abs() / 8499.433742 * 100.,
        (scaling.entropy - -9.019785182).abs() /9.019785182 * 100. ,
        (visc - 0.00125461).abs()/0.00125461 * 100.,
        err);

        let arquivo = OpenOptions::new()
            .create(true)   // cria se não existir
            .append(true)   // escreve no final se existir
            .open("scrit_meu.txt")
            .expect("Erro ao abrir/criar o arquivo");

        use std::io::Write;
        let mut writer = BufWriter::new(arquivo);
        writeln!(writer,"{}\n",content).expect("Erro ao escrever");
        // err = 0.8384488162755742 %
        // err do artigo = 0.820855864 %

    }
    
    #[test]
    fn liquid_hexane_dehlouz_case2() {

        let visc_exp = 0.0000095;

        let parameters = hexane_parameters();

        let cub = utilis::nhexane_dehlouz();
        let eos = E::from_residual(cub);

        let t = 440.;
        let p = 1e5;
        let scrit = -1.290042565;
        // let d = 8499.433742; 

        let state = Arc::new(State::new_tp(eos.into(), t, p, Some(Vapor)).unwrap());

        // let p = state.p;
        let d = state.d;
        
        // assert_relative_eq!(p, 1e5, epsilon = 1e-7); 
        assert_relative_eq!(d, 27.88530923, epsilon = 1e-0); 
        
        let scaling = Scaling::new(&state, parameters);

        assert_relative_eq!(scaling.entropy, -0.019341104, epsilon = 1e-3);
        assert_relative_eq!(scaling.reference, 0.00000193, epsilon = 1e-6);

        assert_relative_eq!(scaling.xscaling(scrit), 4.185205305 , epsilon = 1e-3);
        
        let visc = scaling.viscosity(scrit);
        assert_relative_eq!(visc, 0.00000948 , epsilon = 1e-3);

        println!("err = {} %", (visc - visc_exp).abs() / visc_exp * 100.)
        
        // err = 0.22714210444231467 %
        
        // err ref 0.227287219
    }


    #[test]
    fn water() {


        let parameters = hexane_parameters();

        let cub = utilis::nhexane_dehlouz();
        let eos = E::from_residual(cub);

        let t = 30. + 298.15;
        let p = 1e5;
        let scrit = -1.2845767568055306;
        let state = Arc::new(State::new_tp(eos.into(), t, p, Some(Liquid)).unwrap());
        let scaling = Scaling::new(&state, parameters);
        
        let visc = scaling.viscosity(scrit);

        dbg!(visc);

    }
    
} 