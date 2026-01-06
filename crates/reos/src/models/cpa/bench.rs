

// use nalgebra::{DMatrix, DVector, dvector};
// use ndarray::{Array1, array};
// use ndarray_linalg::Solve;

// use crate::{arr_eq, models::cpa::SCPA};

// use super::parameters::readyto::*;


// #[test]
// fn bench_unbonded_water_co2(){

//     let n = 100_000 ;
//     let t= 298.15;
//     let rho= 27_000.0;
//     let x= &Array1::from_vec(vec![0.5,0.5]);
    
//     let water = water4c();
//     let co2 = co2();
//     let b = water4c_co2();
//     let assoc = SCPA::from_records(vec![water,co2],vec![b]).assoc;
    
//     let volf = &(assoc.rdf.g(rho, x) * &assoc.rdf.bij);

//     // let unbonded = assoc.assoc.x_michelsen(m, kmat).unwrap();
//     let k = &assoc.assoc.association_constants(t, rho, x, volf);
//     let m = &assoc.assoc.m(x);
//     let mult = &assoc.assoc.parameters.multiplicity;
//     let s= m.len();

//     let mult = &assoc.assoc.parameters.multiplicity;

//     let k = &DMatrix::from_fn(s,s,|i,j| k[(i,j)]);
    
//     let m = &DVector::from_fn(s,|j,_| m[j] );
//     let mult = &DVector::from_fn(s,|j,_| mult[j] );
//     // dbg!(k);
    
//     for _ in 0..n {
        
//         // let unbonded = assoc.assoc.x_tan(m,k);
//         // let unbonded = assoc.assoc.x_tan(m,k); //1.08
//         let unbonded = assoc.assoc.x_michelsen(mult,m,k); //1.28

//     }
//     // let unbonded = assoc.assoc.x_michelsen2(m,k).unwrap();
    
//     // let ok = arr_eq!(unbonded, dvector! [0.03874146, 0.20484502, 0.66778819],1e-5);
//     // let ok = arr_eq!(unbonded, dvector! [0.03874146, 0.20484502, 0.66778819],1e-5);
//     // assert!(ok);
// }

// #[test]
// fn bench_unbonded_4c_3b(){

//     let n = 100_000;
//     let t= 273.15;
//     let rho= 6_000.0;
//     let x= &Array1::from_vec(vec![0.5,0.5]);
    
//     let water = water4c();
//     let co2 = methanol3b();
//     // let b = water4c_co2();
//     let assoc = SCPA::from_records(vec![water,co2],vec![]).assoc;
    
//     let volf = &(assoc.rdf.g(rho, x) * &assoc.rdf.bij);

//     // let unbonded = assoc.assoc.x_michelsen(m, kmat).unwrap();
//     let k = &assoc.assoc.association_constants(t, rho, x, volf);
//     let m = &assoc.assoc.m(x);
//     let s= m.len();
//     let mult = &assoc.assoc.parameters.multiplicity;

//     // let k = &DMatrix::from_fn(s,s,|i,j| k[(i,j)]);
    
//     // let m = &DVector::from_fn(s,|j,_| m[j] );
//     // let mult = &DVector::from_fn(s,|j,_| mult[j] );
//     // dbg!(k);
    
//     for _ in 0..n {
        
//         // let unbonded = assoc.assoc.x_tan(m,k);
//         let unbonded = assoc.assoc.x_tan(m,k); //1.71 
//         // let unbonded = assoc.assoc.x_michelsen(mult,m,k); //1.8

//     }
//     // let unbonded = assoc.assoc.x_michelsen(m,k).unwrap();
    
//     // let ok = arr_eq!(unbonded, array! [0.03874146, 0.20484502, 0.66778819],1e-5);
//     // assert!(ok);
// }

// #[test]
// fn bench_solver(){

//     let n = 700_000 * 60;

//     // let g = array![0.23601300218006926, -1.25802240456445, 0.23899129337722336, -1.1982072515705902];
//     // let h = array![[-21.65853364249086, -30.468733820029538, -0.0, -31.528877616985262],
//     //                                             [-30.468733820029538, -661.7102301229936, -31.528877616985262, -0.0],
//     //                                             [-0.0, -31.528877616985262, -20.913081512443068, -27.54595257908517],
//     //                                             [-31.528877616985262, -0.0, -27.54595257908517, -601.6033880260242]];
//     // 0.26
//     let g = -array![0.0006212612107894255, -0.018140418829538874, 0.0006427575348415537, -0.017272511482012476];
//     let h = array! [[-7.466664721913259, -30.468733820029538, -0., -31.528877616985262],
//                                                     [-30.468733820029538, -2157.133832778745, -31.528877616985262, -0.],
//                                                     [-0., -31.528877616985262, -7.226324282621713, -27.54595257908517],
//                                                     [-31.528877616985262, -0., -27.54595257908517, -1958.8241978339759]];

//     // 0.25 -> cte
    
//     for _ in 0..n {

//         // let a= h.clone();
//         let _=  h.solve(&g).unwrap();

//     }
// }


