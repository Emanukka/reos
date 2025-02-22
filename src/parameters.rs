pub mod srk{
    

    use ndarray::{Array1,Array2};
    use std::collections::HashMap;
    use std::fmt;

    pub struct RawSRK{
        pub va0:Vec<f64>,
        pub vb:Vec<f64>,
        pub vkappa:Vec<f64>,
        pub vtc:Vec<f64>,
        pub kij: HashMap<(usize,usize),f64>
    }
    pub struct SRKParameters{
        pub ncomp:  usize,
        pub va0:    Array1<f64>,
        pub vb:     Array1<f64>,
        pub vkappa: Array1<f64>,
        pub vtc:    Array1<f64>,
        pub kij:    Array2<f64>,
    }
    
    
    #[allow(dead_code)]
    impl SRKParameters {
    
        pub fn fill_kij(data:HashMap<(usize,usize),f64>,ncomp:usize)->Array2<f64>{
    
            let mut mkij:Array2<f64> = Array2::zeros((ncomp,ncomp));
            let keys = data.keys();
    
            for k in keys{
    
                let (i,j )= *k;
                let val = *data.get(k).unwrap();
    
                mkij[(i,j)] = val;
                mkij[(j,i)] = val;
                
                // println!("{i},{j} , {val}")
            }
            // println!("{mkij}");
            mkij
        }
    
        pub fn from_raw(ncomp:usize,data:RawSRK)-> Self {
    
            Self{
                ncomp: ncomp,
                va0:    Array1::from_vec(data.va0),
                vb:     Array1::from_vec(data.vb),
                vkappa: Array1::from_vec(data.vkappa),
                vtc:    Array1::from_vec(data.vtc),
                
                kij: Self::fill_kij(data.kij, ncomp)   
    
            }
        }
    
        // Fill kij, w knowledg that kij is Square AND Symmetric
    
        // 
    
    }
    
    
    impl fmt::Display for SRKParameters {
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
    
            write!(f, 
"ncomp ={},
             
kappa ={},

a0  = {},

vb  = {},

vtc = {},

kij = 
{}",
    self.ncomp, self.va0,self.vb,self.vkappa,self.vtc,self.kij)
        }
    }
    
}


pub mod asc{

    use crate::core::core_equations::{NEG,POS};
    use ndarray::{Array1,Array2};
    use std::collections::HashMap;
    use std::{fmt, vec};


    // 1A , 2B , 3B , 4C , SOLV
    // implementar 
    #[allow(dead_code)]
    #[derive(PartialEq)]
    pub enum SchemeType{
        OneSite,
        TwoSites,
        ThreeSites,
        FourSites,
        Solvation,
        Inert,
        TypeDoesntExist
    }

    // pub struct TypeDoesntExist(String);


        impl SchemeType {

        // Dado: 1A , 2B , 3B , 4C , SOLV, 
        // util qndo for possivel extrair dado do json
        
        pub fn from_str(data:&str)-> Self {
            
            match data{
                "1A"=> Self::OneSite,
                "2B"=> Self::TwoSites,
                "3B"=> Self::ThreeSites,
                "4C"=> Self::FourSites,
                "solvation"=> Self::Solvation,
                "inert"=> Self::Inert,
                &_ => panic!("{data} doesn't exist."),
                
            }
        }
        
    }

    // impl SchemeType {

    //     // Dado: 1A , 2B , 3B , 4C , SOLV, 
    //     // util qndo for possivel extrair dado do json
        
    //     pub fn from_str(data:&str)-> Result<Self,TypeDoesntExist> {
            
    //         match data{
    //             "1A"=> Ok(Self::OneSite),
    //             "2B"=> Ok(Self::TwoSites),
    //             "3B"=> Ok(Self::ThreeSites),
    //             "4C"=> Ok(Self::FourSites),
    //             "solvation"=> Ok(Self::Solvation),
    //             &_ => Err(TypeDoesntExist(data.to_string())),
                
    //         }
    //     }
        
    // }
    
    // impl std::fmt::Display for TypeDoesntExist {
    //     fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
    //             write!(f,"{} doens't exist. Ensure that the type passed is: 1A,2B,3B,4C or solvation.",self.0)
    //         }
    //     }


    #[allow(dead_code)]
    pub struct ASCParameters{
        pub ncomp:    usize,
        pub veps:     Array1<f64>,
        pub vb:       Array1<f64>,
        pub vbeta:    Array1<f64>,
        pub s_matrix: Array2<f64>,
        pub associative_components: AssociativeComponents,
        pub non_assoc_comps: NonAssociativeComponents,
        pub solvation: HashMap<(usize,usize),(f64,f64)>,
    }
    #[derive(Debug)]
    pub struct AssociativeComponents{
        pub index: Vec<usize>,
        pub ncomp: usize
    }
    #[derive(Debug)]
    pub struct NonAssociativeComponents{
        pub index: Vec<usize>,
        pub ncomp: usize
    }
    impl ASCParameters {

            pub fn mk_s_matrix(ncomp:usize,scheme: &Vec<SchemeType>)->Array2<f64>{

                // modificar pra diminuir no futuro
                let mut s_mat:Array2<f64> = Array2::zeros((2,ncomp));


                for (i,s) in scheme.iter().enumerate(){

                    match s {
                        SchemeType::OneSite => s_mat[(POS,i)] = 1.0,

                        SchemeType::TwoSites => {
                            s_mat[(POS,i)] = 1.0;
                            s_mat[(NEG,i)] = 1.0;
                        }
                        SchemeType::ThreeSites=> {
                            s_mat[(NEG,i)] = 2.0;
                            s_mat[(POS,i)] = 1.0;
                        }
                        SchemeType::FourSites=> {
                            s_mat[(NEG,i)] = 2.0;
                            s_mat[(POS,i)] = 2.0;
                        }
                        SchemeType::Solvation=> {
                            s_mat[(NEG,i)] = 1.0;
                        }
                        SchemeType::Inert=> {}

                        SchemeType::TypeDoesntExist=>
                        {
                            panic!("Found a type that doens't exist in scheme.")
                        }
                    }
                }
                s_mat
            }

            pub fn define_association_comps(scheme: &Vec<SchemeType>)->(AssociativeComponents,NonAssociativeComponents){

                let mut assoc: Vec<usize> = Vec::new();
                let mut inert: Vec<usize> = Vec::new();
                for (i,s) in scheme.iter().enumerate(){

                    if s == &SchemeType::Inert{
                        inert.push(i);
                    }

                    else {
                        assoc.push(i);
                    }
                }

                (AssociativeComponents{ncomp:assoc.len(),index:assoc}, NonAssociativeComponents{ncomp:inert.len(),index:inert})

            }



        }

        impl fmt::Display for ASCParameters {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        
                write!(f, 
                    
" ncomp = {},
        
vb = {},
        
veps = {},
        
vbeta = {}, 
        
S =
{},

Solvation Map = 
{:#?}",
                    
                self.ncomp, self.vb,self.veps,self.vbeta,self.s_matrix,self.solvation)
            }
        }

        
    }

