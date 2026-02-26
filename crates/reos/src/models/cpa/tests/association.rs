use approx::assert_relative_eq;
use ndarray::{Array2, array};
use quantity::{NAV, MOL};
use super::*;

const TOL:f64 = 1e-12;

// Kontogeorgis radial distribution function, see https://doi.org/10.1021/ie051305v
mod kontogeorgis {

    use super::*;
    
    #[test]
    fn test_4c_unbonded_sites(){

        let t=298.15;
        let rho= 58_000.0;

        let x= &array![1.0];
        let water = pure::water();
        let assoc = water.assoc;
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.b;
        // let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = assoc.assoc.association_constants(t, rho, x,&volf);

        assert_relative_eq!(k,  rho * *(NAV * MOL) * array![[0., 2.2974432513007877e-27],
                                    [2.2974432513007877e-27, 0.]], epsilon = TOL);

        assert_relative_eq!(assoc.assoc.unbonded_sites_fraction(x, &k), array![0.07588168172556269, 0.07588168172556281], epsilon = TOL);
    }


    #[test]
    fn test_4c_3b_unbonded_sites(){

        let t = 298.15;
        let rho = 1_000.0;
        let x = &array![0.6, 0.4];

        let water = recipes::water4c();
        let methanol = recipes::methanol3b();

        let cpa = recipes::scpa(vec![water, methanol],vec![]).unwrap();

        let assoc = cpa.assoc;

        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.b;

        let k = assoc.assoc.association_constants(t, rho, x, &volf);
        let m = assoc.assoc.sites_mole_frac(x);
        let mm = Array2::from_shape_fn((m.len(),m.len()), |(i,j)| m[i] * m[j]);

        let k_clapeyron = mm * rho * *(NAV * MOL) * array![

            [0.0, 1.3933286407053757e-27, 0.0, 1.4407698427189598e-27],
            [1.3933286407053757e-27, 0.0, 1.4407698427189598e-27, 0.0],
            [0.0, 1.4407698427189598e-27, 0.0, 1.2583687227930472e-27],
            [1.4407698427189598e-27, 0.0, 1.2583687227930472e-27, 0.0]
            
        ];

        assert_relative_eq!(k,  k_clapeyron, epsilon = TOL);

        assert_relative_eq!(assoc.assoc.unbonded_sites_fraction(x, &k), array![0.5974059887556956, 0.4956687564644897 , 0.5992454354468074, 0.5037025677672323], epsilon = 1e-9);


    }

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

    #[test]
    fn test_4c_induced_association() {
        
        let t=298.15;
        let rho= 1_000.;
        let x= &array![0.7,0.3];
        let comp1 = super::recipes::water4c();
        let comp2 = super::recipes::co2();
        let b = super::recipes::water4c_co2();
        
        let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.b;
    // let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = assoc.assoc.association_constants(t, rho, x, &volf);
        let m = assoc.assoc.sites_mole_frac(x);
        let mm = Array2::from_shape_fn((m.len(),m.len()), |(i,j)| m[i] * m[j]);

        let k_clapeyron = mm * rho * *(NAV * MOL) * array![
            [0.0, 1.390825988133717e-27, 1.7812705976946612e-28],
            [1.390825988133717e-27, 0.0, 0.0],
            [1.7812705976946612e-28, 0.0, 0.0],

            
        ];

        assert_relative_eq!(k, k_clapeyron, epsilon =  TOL);
        assert_relative_eq!(assoc.assoc.unbonded_sites_fraction(x, &k), array![0.5786343099065254, 0.5957666698192122, 0.9200489870741287], epsilon = 1e-9);

    }



        #[test]
        fn test_4c_1a_elliot_cr() {

            let t=298.15;
            let rho= 1_000.0;
            let x= &array![0.75,0.25];

            let comp1 = super::recipes::water4c();
            let comp2 = super::recipes::acetic1a();
            let b = super::recipes::water4c_acetic1a();
            let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;

            // println!("{:?}",assoc.assoc.parameters.interactions);
            let volf = assoc.rdf.g(rho, x) * &assoc.rdf.b;
            let k_rust = assoc.assoc.association_constants(t, rho, x, &volf);
            // dbg!(assoc.rdf.g(rho, x));

            let kref = array![
                [  0.              ,   0.47210091350445,   8.53435728830102],
                [  0.47210091350445,   0.              ,   8.53435728830102],
                [  8.53435728830102,   8.53435728830102, 154.27899468296786]];

            println!("{}",k_rust);
            assert_relative_eq!(k_rust, kref, epsilon = TOL);
        }

        // #[test]
        // fn wa() {
        //     let t=298.15;
        //     let rho= 21_500.0;
        //     let x= &array![0.25,0.75];

        //     let comp1 = super::recipes::water4c();
        //     let comp2 = super::recipes::acetic1a();
        //     let b = super::recipes::water4c_acetic1a();
        //     let assoc = super::recipes::scpa(vec![comp1,comp2],vec![b],).unwrap().assoc;

        //     let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        //     let k = assoc.assoc.association_constants(t, rho, x, &volf);

        //     let unb = assoc.assoc.unbonded_sites_fraction(x,&k);

        //     let reff = array![0.1598388 , 0.1598388 , 0.00241471];

        //     assert_relative_eq!(unb,reff,epsilon=1e-6);

        // }


    }

// Carnahan-Starling radial distribution function, see https://doi.org/10.1021/ie051305v
mod carnahan_starling {

    use super::*;

    #[test]
    fn test_4c_unbonded_sites(){

        let t = 298.15;
        let rho = 58_000.0;
        let x= &array![1.0];

        let water = pure::water_cs();
        let assoc = water.assoc;
        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.b;
        // let volf = assoc.rdf.g(rho, x) * &assoc.rdf.bij;
        let k = assoc.assoc.association_constants(t, rho, x,&volf);

        assert_relative_eq!(k,  rho * *(NAV * MOL) * array![[0., 2.506550632893448e-27],
                                    [2.506550632893448e-27, 0.]], epsilon = TOL);

        assert_relative_eq!(assoc.assoc.unbonded_sites_fraction(x, &k), array![0.07276978646379922, 0.07276978646379931], epsilon = TOL);
    }

    #[test]
    fn test_4c_3b_unbonded_sites(){

        let t = 298.15;
        let rho = 1_000.0;
        let x= &array![0.6, 0.4];

        let water = recipes::water4c();
        let methanol = recipes::methanol3b();

        let cpa = recipes::cpa_cs(vec![water, methanol],vec![]).unwrap();
        let assoc = cpa.assoc;

        let volf = assoc.rdf.g(rho, x) * &assoc.rdf.b;

        let k = assoc.assoc.association_constants(t, rho, x, &volf);

        let m = assoc.assoc.sites_mole_frac(x);
        let mm = Array2::from_shape_fn((m.len(),m.len()), |(i,j)| m[i] * m[j]);

        let k_clapeyron = mm * rho * *(NAV * MOL) * array![

            [0.0, 1.3979300759188778e-27, 0.0, 1.445527951391357e-27],
            [1.3979300759188778e-27, 0.0, 1.445527951391357e-27, 0.0],
            [0.0, 1.445527951391357e-27, 0.0, 1.2625244560375012e-27],
            [1.445527951391357e-27, 0.0, 1.2625244560375012e-27, 0.0]
            
        ];

        assert_relative_eq!(k,  k_clapeyron, epsilon = TOL);

        assert_relative_eq!(
            assoc.assoc.unbonded_sites_fraction(x, &k), 
            array![0.596911541557299, 0.49505113889610836, 0.5987518789366115, 0.5030849658567943], epsilon = 1e-9);

    }
}