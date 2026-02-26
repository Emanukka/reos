use approx::assert_relative_eq;
use ndarray::array;

use crate::{residual::Residual};

use super::recipes;





const TOL:f64 = 1e-12;

const T:f64 = 298.15;
const RHO:f64 = 1_000.;

mod kontogeorgis {

    use super::*;
    
    pub fn water_methanol()->super::super::super::CPA{

        let water = recipes::water4c();
        let methanol = recipes::methanol3b();
        recipes::scpa(vec![water, methanol], vec![]).unwrap()

    }

    #[test]
    fn test_helmholtz() {
        
        let model = water_methanol();
        let x = array![0.6, 0.4];
        let val = model.helmholtz(T, RHO, &x);
        assert_relative_eq!(val, -1.4579443886240895 , epsilon = TOL)
    }
    #[test]
    fn test_df_dt() {
        
        let model = water_methanol();
        
        let x = array![0.6,0.4];
        let val = model.df_dt(T, RHO, &x);
        assert_relative_eq!(val, 0.018638018571771893, epsilon = TOL)
    }
    #[test]
    fn test_df_dn() {
        
        let model = water_methanol();
        let x = array![0.6,0.4];
        let val = model.df_dn(T, RHO, &x);
        assert_relative_eq!(val, array![-2.6026409814278524, -2.0602656114888185] , epsilon = 1e-9)
    }

    #[test]
    fn test_compressibility() {
        
        let model = water_methanol();
        let x = array![0.6,0.4];
        let val = - model.df_dv(T, RHO, &x) / RHO;
        
        
        assert_relative_eq!(val, -0.9277464448281495 , epsilon = TOL)
    }

}

mod carnahan_starling {
    
    use super::*;

    pub fn water_methanol()->super::super::super::CPA{

        let water = recipes::water4c();
        let methanol = recipes::methanol3b();
        recipes::cpa_cs(vec![water, methanol], vec![]).unwrap()

    }
    #[test]
    fn test_helmholtz() {
        
        let model = water_methanol();
        let x = array![0.6,0.4];
        let val = model.helmholtz(T, RHO, &x);
        assert_relative_eq!(val, -1.4605958993450774 , epsilon = TOL)
    }
    #[test]
    fn test_df_dt() {
        
        let model = water_methanol();
        
        let x = array![0.6,0.4];
        let val = model.df_dt(T, RHO, &x);
        assert_relative_eq!(val,  0.018660061178296765, epsilon = TOL)
    }
    #[test]
    fn test_df_dn() {
        
        let model = water_methanol();
        let x = array![0.6,0.4];
        let val = model.df_dn(T, RHO, &x);
        assert_relative_eq!(val, array![-2.6085353063412935, -2.0671551322928976] , epsilon = 1e-9)

    }

    #[test]
    fn test_compressibility() {
        
        let model = water_methanol();
        let x = array![0.6,0.4];
        let val = - model.df_dv(T, RHO, &x) / RHO;
        
        
        assert_relative_eq!(val, -0.9313873373768573 , epsilon = TOL)
    }
}
