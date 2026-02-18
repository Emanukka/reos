
use super::{recipes, super::{Associative,AssociativeParameters}};
use crate::parameters::Parameters;

use ndarray::array;
use approx::assert_relative_eq;


fn water() -> Associative {
    let pure_record = recipes::water();
    let parameters = AssociativeParameters::build_pure(pure_record, ()).unwrap();
    Associative::from_parameters(parameters)
}
fn methanol() -> Associative {
    let pure_record = recipes::methanol();
    let parameters = AssociativeParameters::build_pure(pure_record, ()).unwrap();
    Associative::from_parameters(parameters)
}


mod unbonded_sites{


    use super::*;

    #[test]
    fn test_3b(){
        // t = 298.15 K , d = 1000.0
        let asc = methanol();
        let k = array![[0.           , 0.76195179893041],
                                                    [0.76195179893041, 0.           ]];
                                                    
        let x = array![1.0];
        let unb = asc.unbonded_sites_fraction(&x, &k);
        let reff = array![0.73571946249752, 0.47143892500497];
        unb.iter().zip(reff.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });
        
        let m = &asc.sites_mole_frac(&x);
        let tan = asc.x_tan(m, &k).unwrap();
        
        unb.iter().zip(tan.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });
    }

    #[test]
    fn test_4c(){
        // t = 298.15 K , d = 1000.0
        let asc = water();
        let k = array![[0.           , 0.83518048731],
                                                    [0.83518048731, 0.           ]];
        let x = array![1.0];
        let unb = asc.unbonded_sites_fraction(&x, &k);
        let reff = array![0.530287110928 , 0.530287110928];
        unb.iter().zip(reff.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });
        
        let m = &asc.sites_mole_frac(&x);
        let tan = asc.x_tan(m, &k).unwrap();
        
        unb.iter().zip(tan.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });
    }

    #[test]
    fn test_1a_with_inert() {

        let asc = water();
        let k = array![[0.           , 0.83518048731],
                                                    [0.83518048731, 0.           ]];
        let x = array![1.0];
        let unb = asc.unbonded_sites_fraction(&x, &k);
        let reff = array![0.530287110928 , 0.530287110928];
        unb.iter().zip(reff.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });
        
        let m = &asc.sites_mole_frac(&x);
        let tan = asc.x_tan(m, &k).unwrap();
        
        unb.iter().zip(tan.iter()).for_each(|(x,y)| {
            assert_relative_eq!(x, y, epsilon = 1e-10)
        });
    }
}


#[test]
fn test_association_helmholtz() {
    let asc = water();
    let x = array![1.0];
    let unbonded = array![0.530287110928 , 0.530287110928];
    
    let a = asc.r_helmholtz(&x, &unbonded);
    assert_relative_eq!(a, -1.597921023379, epsilon = 1e-10)
}
#[test]
fn test_association_entropy() {
    let t = 298.15;
    let asc = water();
    let x = array![1.0];
    let k = array![[0.           , 0.83518048731],
                                                 [0.83518048731, 0.           ]];
    let s = asc.r_entropy(t, &x, &k);
    assert_relative_eq!(s, -4.713659269705, epsilon = 1e-9)
}
#[test]
fn test_association_pressure() {
    let d = 1000.0;
    let asc = water();
    let x = array![1.0];
    let unbonded = array![0.530287110928 , 0.530287110928];
    let h = asc.h(&x, &unbonded);
    let p = asc.r_pressure(h, d, 6.9352666490453005 * 1e-6);
    assert_relative_eq!(p, -945.9409464127781, epsilon = 1e-9)
}
#[test]
fn test_association_chem_pot() {
    let asc = water();
    let x = array![1.0];
    let unbonded = array![0.530287110928 , 0.530287110928];
    let h = asc.h(&x, &unbonded);
    
    let mu = asc.r_chemical_potential(h, &array![0.00693526664905], &unbonded);
    let reff = array![-2.54386196979185, -2.54386196979185];
    mu.iter().zip(reff.iter()).for_each(|(x,y)| {
        assert_relative_eq!(x, y, epsilon = 1e-10)
    });
}