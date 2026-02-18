use approx::assert_relative_eq;
use ndarray::array;


    #[test]
    fn test_assoc_consts_4c(){

        let t=298.15;
        let rho= 1000.0;

        let x= &array![1.0];
        // let comp1 = water4c();
        let water = super::pure::water();
        let assoc = water.assoc;
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;

        let k = assoc.assoc.association_constants(t, rho, x,&volf);
        let reff = array![[0.           , 0.83518048731],[0.83518048731, 0.           ]];

        assert_relative_eq!(k, reff, epsilon = 1e-10)        
        // k.iter().zip(reff.iter()).for_each(|(x,y)| {
        //     assert_relative_eq!(x, y, epsilon = 1e-10)
        // });

    }

    // #[test]
    // fn test_assoc_consts_3b(){

    //     let t=298.15;
    //     let rho= 1000.0;

    //     let x= &array![1.0];
    //     let comp1 = super::recipes::methanol3b();
    //     let assoc = super::recipes::scpa(vec![comp1],vec![],).unwrap().assoc;
    //     let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;

    //     let k = assoc.assoc.association_constants(t, rho, x,&volf);
    //     let reff = array![[0.           , 0.76195179893041],
    //                                                  [0.76195179893041, 0.           ]];  

    //     k.iter().zip(reff.iter()).for_each(|(x,y)| {
    //         assert_relative_eq!(x, y, epsilon = 1e-10)
    //     });

    // }


    // #[test]
    // fn test_assoc_consts_1a_inert(){
    //     let t=298.15;
    //     let rho= 7_500.;
    //     let x= &array![0.25,0.75];

    //     let comp1 = super::recipes::acetic1a();
    //     let comp2 = super::recipes::octane();
    //     let b = super::recipes::acoh_octane();
    //     let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;
        
    //     let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
    //     let k_rust = assoc.assoc.association_constants(t, rho, x,&volf);

    //     let ok = assert_relative_eq!(k_rust, array![ 1980.9483616 ], epsilon=1e-7);

    //     assert!(ok)
    // }
    // #[test]
    // fn unbonded_acetic_octane(){
    //     let t=298.15;
    //     let rho= 7_500.;
    //     let x= &array![0.25,0.75];

    //     let comp1 = super::recipes::acetic1a();
    //     let comp2 = super::recipes::octane();
    //     let b = super::recipes::acoh_octane();
    //     let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;
        
    //     let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
    //     let k = assoc.assoc.association_constants(t, rho, x,&volf);

    //     let m = &assoc.assoc.sites_mole_frac(x);
    //     let tan = assoc.assoc.x_tan(m, &k).unwrap();
    //     // let michelsen = assoc.assoc.x_michelsen(m, &k).unwrap();
    //     let analytic = assoc.assoc.unbonded_sites_fraction(x,&k);
        
    //     let ok = assert_relative_eq!(tan, analytic, epsilon=1e-10);


    //     assert!(ok)
    // }

    // #[test]
    // fn phi_acetic_octane(){
    //     let t=298.15;
    //     let rho= 7_500.;
    //     let x= &array![0.25,0.75];

    //     let comp1 = super::recipes::acetic1a();
    //     let comp2 = super::recipes::octane();
    //     let b = super::recipes::acoh_octane();
    //     // let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b]).unwrap().assoc;
        
    //     let cpa = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap();

    //     let s = S::new_trx(Arc::new(E::from_residual(cpa)), t, rho, x.clone());
    //     let phi_rust = s.lnphi().exp();
        
    //     let ok = assert_relative_eq!(phi_rust,array![0.0002995 , 0.00071522],epsilon=1e-6);

    //     assert!(ok)
    // }
    #[test]
    fn association_constants_water_co2() {

        let t=298.15;
        let rho= 55_000.0;
        let x= &array![0.5,0.5];

        let comp1 = super::recipes::water4c();
        let comp2 = super::recipes::co2();
        let b = super::recipes::water4c_co2();
        let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k_rust = assoc.assoc.association_constants(t, rho, x,&volf).flatten().to_owned();

        assert_relative_eq!(k_rust, array![ 0. , 25.04896564,  3.21025658, 25.04896564,  0., 0., 3.21025658,  0.,  0.], epsilon=1e-7);

    }



    #[test]
    fn association_constants_water_acetic() {

        let t=298.15;
        let rho= 21_500.0;
        let x= &array![0.25,0.75];

        let comp1 = super::recipes::water4c();
        let comp2 = super::recipes::acetic1a();
        let b = super::recipes::water4c_acetic1a();
        let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k_rust = assoc.assoc.association_constants(t, rho, x, &volf).flatten().to_owned();

        let reff =  array![0.00000000e+00, 1.84368158e+00, 3.00115845e+02,
                                              1.84368158e+00, 0.00000000e+00, 3.00115845e+02, 
                                              3.00115845e+02, 3.00115845e+02, 4.88530782e+04] ;

        assert_relative_eq!(k_rust,reff,epsilon=1e-4);
    }

    #[test]
    fn unbonded_water_acetic() {
        let t=298.15;
        let rho= 21_500.0;
        let x= &array![0.25,0.75];

        let comp1 = super::recipes::water4c();
        let comp2 = super::recipes::acetic1a();
        let b = super::recipes::water4c_acetic1a();
        let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = assoc.assoc.association_constants(t, rho, x, &volf);

        let unb = assoc.assoc.unbonded_sites_fraction(x,&k);

        let reff = array![0.1598388 , 0.1598388 , 0.00241471];

        assert_relative_eq!(unb,reff,epsilon=1e-6);

    }

    #[test]
    fn permutation_water_acetic() {

        let t=298.15;
        let rho= 21_500.0;
        let x= &array![0.75, 0.25];

        let comp1 = super::recipes::water4c();
        let comp2 = super::recipes::acetic1a();
        let b = super::recipes::water4c_acetic1a();
        let assoc = super::recipes::scpa(vec![comp2,comp1],vec![b],).unwrap().assoc;
        
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;

        let k = assoc.assoc.association_constants(t, rho, x, &volf);
        let unb = assoc.assoc.unbonded_sites_fraction(x,&k);

        let reff = array![0.00241471, 0.1598388 , 0.1598388 ];

        assert_relative_eq!(unb, reff, epsilon=1e-6);

    }
    // #[test]
    // fn phi_water_co2(){

    //     let p=500e5;
    //     let t=298.15;
    //     let x=Array1::from_vec(vec![0.5,0.5]);
    //     let water = water4c();
    //     let co2 = co2();
    //     let b = water4c_co2();
    //     let cpa = super::recipes::scpa(vec![water,co2],vec![b],).unwrap();
        
    //     let s = S::new_tpx(Arc::new(E::from_residual(cpa)), t, p, x.clone(), Some(DensityInitialization::Vapor)).unwrap();
    //     let phi_rust = s.lnphi().exp();
        
    //     let ok = assert_relative_eq!(phi_rust,array![2.14385745e-04, 5.65853284e-01],epsilon=1e-6);

    //     assert!(ok)


    // }

