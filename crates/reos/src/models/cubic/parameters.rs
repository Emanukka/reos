use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize};
use crate::{models::cubic::{alpha::{Alpha, AlphaParameter}, models::{ CubicModel, CubicModels}},
parameters::{ BinaryMap, Parameters, records::{BinaryParameter, Properties}}};


type Pure = CubicPureRecord;
type Binary = CubicBinaryRecord;
    

// #[derive(Clone,Debug)]
// pub struct CubicParameters {
//     pub ncomp: usize,
//     pub a: Array1<f64>,
//     pub b: Array1<f64>,
//     pub tc: Array1<f64>,
//     pub alpha: Alpha,
//     pub vvolt: Array1<f64>,
//     pub properties: Properties,
//     pub binary:BinaryMap<CubicBinaryRecord>
// }
#[derive(Clone,Debug)]
pub struct CubicParameters {
    pub ncomp: usize,
    pub a: Array1<f64>,
    pub b: Array1<f64>,
    pub tc: Array1<f64>,
    pub aij: Array2<f64>,
    pub bij: Array2<f64>,
    pub epsilon: f64,
    pub binary:BinaryMap<CubicBinaryRecord>,
    pub sigma: f64,
    pub alpha: Alpha,
    pub vvolt: Array1<f64>,
    pub properties: Properties,
}

impl Parameters<Pure, Binary, CubicModels> for CubicParameters{

    fn from_raw(pure:Vec<Pure>, binary: crate::parameters::BinaryMap<Binary>, properties: Option<Properties>, opt: CubicModels) -> Self {
        
    // fn from_raw(pure:Vec<Pure>, binary: Vec<BinaryParameter<Binary>>, properties: Option<Properties>, opt: CubicModels) -> Self {
        let model = opt;


        let n = pure.len();
        let mut ma0    = Array1::<f64>::zeros(n);
        let mut mb     = Array1::<f64>::zeros(n);
        let mut mtc    = Array1::<f64>::zeros(n);
        let mut mkappa = Array1::<f64>::zeros(n);
        let mut maij   = Array2::<f64>::zeros((n,n));
        let mut mbij   = Array2::<f64>::zeros((n,n));
        let mut vvolt = Array1::<f64>::zeros(n);
        let mut alpha_parameters = Vec::<AlphaParameter>::with_capacity(n);

        for (i,record) in pure.iter().enumerate(){

            match record {
                
                &CubicPureRecord::Set1 { a0, b, c1, tc, volt } =>
                {

                    let volt = volt.unwrap_or_default();
                    ma0[i] = a0;
                    mb[i] = b;
                    mkappa[i] = c1;
                    mtc[i] = tc;

                    vvolt[i] = volt;

                    alpha_parameters.push(AlphaParameter::Soave(c1));
                }
                
                &CubicPureRecord::Set2 { tc, pc, w , volt} =>
                {
                    let volt = volt.unwrap_or_default();
                    ma0[i] = model.acrit(tc, pc);
                    mb[i] = model.bcrit(tc, pc);
                    mtc[i] = tc;
                    vvolt[i] = volt;



                    let c1 = model.kappa_from_w(w);
                    alpha_parameters.push(AlphaParameter::Soave(c1));

                }

                &CubicPureRecord::Twu91 { tc, pc, l, n, m, volt } => 
                {   
                    let volt = volt.unwrap_or_default();
                    ma0[i] = model.acrit(tc, pc);
                    mb[i] = model.bcrit(tc, pc);
                    mtc[i] = tc;
                    vvolt[i] = volt;

                    alpha_parameters.push(AlphaParameter::Twu91 { l, n, m });

                    
                }
            }
        }

        binary.iter().for_each(|(key,b)|{
            
            let comp1 = key.0;
            let comp2 = key.1;

            maij[(comp1,comp2)] = b.kij; 
            maij[(comp2,comp1)] = b.kij;
            
            mbij[(comp1,comp2)] = b.lij; 
            mbij[(comp2,comp1)] = b.lij;
            // match b{
            //     &CubicBinaryRecord::TemperatureDependent { aij, bij } => {

            //         maij[(comp1,comp2)] = aij; 
            //         maij[(comp2,comp1)] = aij;
                    
            //         mbij[(comp1,comp2)] = bij; 
            //         mbij[(comp2,comp1)] = bij;

            //     }
            //     &CubicBinaryRecord::TemperatureIndependent { kij } => {
            //         maij[(comp1,comp2)] = kij; 
            //         maij[(comp2,comp1)] = kij;
            //     }
            // }
        });

        let alpha = Alpha::new(alpha_parameters);

        // Self{
        //     alpha,
        //     ncomp:n,
        //     a: ma0,
        //     b: mb,
        //     tc: mtc,
        //     vvolt: vvolt,
        //     properties: properties.unwrap_or_default(),
        //     binary
        // }
        Self{
            binary,
            alpha,
            ncomp:n,
            a: ma0,
            b: mb,
            tc: mtc,
            vvolt: vvolt,
            properties: properties.unwrap_or_default(),
            aij:maij,
            bij:mbij,
            sigma:model.sig(),
            epsilon:model.eps()
        }


    }
}


// #[derive(Clone,Debug,Serialize, Deserialize)]
// #[serde(untagged)]
// pub enum CubicBinaryRecord{
//     TemperatureDependent{aij:f64, bij:f64},
//     TemperatureIndependent{kij:f64}
// }
#[derive(Clone,Debug,Serialize, Deserialize)]
pub struct CubicBinaryRecord{
    pub kij:f64,
    #[serde(default)]
    pub lij:f64
}

impl Default for CubicBinaryRecord {
    fn default() -> Self {
        Self { kij: 0., lij: 0. }
    }
} 
#[derive(Clone,Debug,Serialize, Deserialize)]
#[serde(untagged)]
pub enum CubicPureRecord{

    Set1{
        a0:f64,
        b:f64,
        c1:f64,
        tc:f64,
        #[serde(default)]
        volt:Option<f64>
    },

    Set2{
        tc:f64,
        pc:f64,
        w:f64,
        #[serde(default)]
        volt:Option<f64>
    },

    Twu91 {
        tc:f64,
        pc:f64,
        l:f64,
        n:f64,
        m:f64,
        #[serde(default)]
        volt:Option<f64>
    },

}




impl CubicPureRecord {
    
    pub fn new_set1(a0:f64, b:f64, c1:f64, tc:f64, volt:Option<f64>)->Self{

        Self::Set1{
            a0,
            b,
            c1,
            tc,
            volt 
        }
    }
    pub fn new_set2(tc:f64,pc:f64,w:f64, volt:Option<f64>)->Self{

        Self::Set2{
            tc,
            pc,
            w,
            volt
        
        }
    }

    pub fn twu91(tc:f64, pc:f64, l:f64, n:f64, m:f64, volt:Option<f64>) -> Self {

        Self::Twu91 { tc, pc, l, n, m, volt }
    }
}

impl std::fmt::Display for CubicParameters {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        unimplemented!()
        // if self.ncomp == 1 {
        //     write!(f, "CubicParameters(\n\ta0={},\n\tb={},\n\talpha={},\n\ttc={})",
        //         self.a.to_string(),
        //         self.b.to_string(),
        //         self.alpha.to_string(),
        //         self.tc.to_string(),
        //     )

        // } else {

        //     let saij = self.aij
        //     .rows()
        //     .into_iter()
        //     .map(|row| row.to_string())
        //     .collect::<Vec<_>>()
        //     .join("\n\t  ");
        //     let sbij = self.bij
        //     .rows()
        //     .into_iter()
        //     .map(|row| row.to_string())
        //     .collect::<Vec<_>>()
        //     .join("\n\t  ");

        //     write!(f, "CubicParameters(\n\ta0={},\n\tb={},\n\talpha={},\n\ttc={},\n\taij=\n\t  [{}],\n\tbij=\n\t  [{}])",
        //         self.a.to_string(),
        //         self.b.to_string(),
        //         self.alpha.to_string(),
        //         self.tc.to_string(),
        //         saij,
        //         sbij
        //     )
        // }

    }
}
// impl PureRecord for CubicPureRecord {
    
// // }
// impl<T:CubicModel> Parameters<CubicPureRecord,CubicBinaryRecord,String> for CubicParameters<T> {
    
//     fn from_recordss(pure_records:Vec<PureRecord<CubicPureRecord>>,binary_records: Vec<crate::parameters::records::BinaryRecord<CubicBinaryRecord,String>>)->Self {

//         todo!()
//     }
    
// }
impl std::fmt::Display for CubicPureRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        match self {
            
            Self::Set1 { a0, b, c1, tc , volt} => {
                write!(f, "CubicPureRecord(a0={}, b={}, c1={}, tc={}, volt={})", a0, b, c1, tc, volt.unwrap_or_default())
            }
            Self::Set2 { tc, pc, w, volt } => {
                write!(f, "CubicPureRecord(tc={}, pc={}, w={}, volt={})", tc, pc, w, volt.unwrap_or_default())
            }
            Self::Twu91 { tc, pc, l, n, m, volt } => {
                write!(f, "CubicPureRecord(tc={}, pc={}, l={}, n={}, m={}, c={})", tc, pc, l, n, m, volt.unwrap_or_default())

            }
        }
    }
}

impl std::fmt::Display for CubicBinaryRecord {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        unimplemented!()
        // match self {
            
        //     Self::TemperatureDependent { aij, bij } => {
        //         write!(f, "CubicBinaryRecord(aij={}, bij={})", aij, bij   )
        //     }
        //     Self::TemperatureIndependent { kij } => {
        //         write!(f, "CubicBinaryRecord(kij={})", kij)
        //     }
        // }
    }
}


// impl<T:CubicModel > CubicParameters<T>{


//     pub fn set_kij(&mut self,i:usize,j:usize,kij:f64){

//         self.aij[(i,j)] = kij;
//         self.aij[(j,i)] = kij;
//     }
//     pub fn set_kij_temperature_dependent(&mut self,i:usize,j:usize,aij:f64,bij:f64){

//         self.aij[(i,j)] = aij;
//         self.aij[(j,i)] = aij; 
 
//         self.bij[(i,j)] = bij;
//         self.bij[(j,i)] = bij;
//     }

// }

// impl<C:CubicModel> Parameters<CubicPureRecord> for CubicParameters<C> {
//     fn from_records(records:Vec<CubicPureRecord>) -> Self {
        
//         Self::new(C::model(), records)
//     }
// }

// impl<C:CubicModel > CubicParameters<C> {
//     pub fn from_records(pure_records:Vec<PureRecord<CubicPureRecord>>,binary_records:Vec<BinaryRecord<CubicBinaryRecord>>) -> Self {
        
//         Self::new(pure_records, binary_records)
//     }

// }


// impl <T:CubicModel> Deserialize for CubicParameters<T> {
//     fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
//         where
//             D: serde::Deserializer<'de> {
        
//     }
// }

#[cfg(test)]
mod tests{
    use ndarray::array;
    use serde_json::from_str;

    // use crate::{models::cubic::{Cubic, models::CubicModels, parameters::{CubicBinaryRecord, CubicPureRecord}}, parameters::{Parameters, records::{BinaryRecord, PureRecord}}, residual::{Residual, ResidualDerivedProperties}};

    use crate::residual::{Residual, ResidualDerivedProperties};

    use super::*;
    use super::super::Cubic;
    use super::super::models::{PR76,PR78};
    use crate::parameters::{PureRecord,BinaryRecord};
    #[test]
    fn test_pure_records_from_json(){
        let data1 = r#"
            {   
                "name": "resin",
                "tc":507.317,
                "pc":3290000,
                "w":0.275000
            }
        "#;
        let data2 = r#"
            {   
                "name": "asphaltene",
                "tc": 860.372, 
                "pc": 1246000, 
                "w":  1.02200 
            }
        "#;
        let pr1: PureRecord<CubicPureRecord> = from_str(data1).unwrap();
        let pr2: PureRecord<CubicPureRecord> = from_str(data2).unwrap();
        // let json = serde_json::to_string_pretty(&pr1).unwrap();
        
        // pr78
        // let params = super::CubicParameters::<PR78>::new(vec![c1,c2],vec![]);
        let model = PR78.into();
        let params = CubicParameters::new(vec![pr1,pr2],vec![], model);

        // let json = serde_json::to_string_pretty(&params).unwrap();

        
        // println!("{}", json);

        let cub = Cubic::from_parameters(params);
        let t = 298.15;
        let rho = 1000.0;
        let x = array![0.5,0.5];
        let test = cub.all_derived_properties(t, rho, &x);
        
        let data = r#"
            {
              "a": -6.739499683331773,
              "dadv": -5538.109023777433,
              "dadni": [
                -4.772941005259919,
                -19.78227640895849
              ]
            }
            "#;

        let answer:ResidualDerivedProperties = from_str(data).unwrap();

        answer.comparison(test,None);
        
    }

    #[test]
    fn test1(){
        let tc = array![ 507.317  , 860.372];
        let pc = array![ 3290000. , 1246000.];
        let w  = array![ 0.27500  , 1.022  ];

        let c1 = CubicPureRecord::new_set2(tc[0], pc[0], w[0], None);
        let c2 = CubicPureRecord::new_set2(tc[1], pc[1], w[1], None);

        let pr1 = PureRecord::new(0.0, "".into(), c1);
        let pr2 = PureRecord::new(0.0, "".into(), c2);        

        let model = PR78.into();
        let params = CubicParameters::new(vec![pr1,pr2],vec![], model);
        // println!("{}",params);
        // let json = serde_json::to_string_pretty(&params).unwrap();

        
        // println!("{}", json);

        let cub = Cubic::from_parameters(params);
        let t = 298.15;
        let rho = 1000.0;
        let x = array![0.5,0.5];
        let test = cub.all_derived_properties(t, rho, &x);
        // let json = serde_json::to_string_pretty(&r_properties).unwrap();
        // println!("{}",json);
        let data = r#"
            {
              "a": -6.739499683331773,
              "dadv": -5538.109023777433,
              "dadni": [
                -4.772941005259919,
                -19.78227640895849
              ]
            }
            "#;

        let answer:ResidualDerivedProperties = from_str(data).unwrap();

        answer.comparison(test,None);
        
    }
    #[test]
    fn test2(){

        // let tc = array![ 507.317  , 860.372 ];
        // let pc = array![ 3290000. , 1246000.];
        // let w  = array![ 0.27500  , 1.02200 ];

        let data1 = r#"
            {   
                "tc":507.317,
                "pc":3290000,
                "w":0.275000
            }
        "#;
        let data2 = r#"
            {   
                "tc": 860.372, 
                "pc": 1246000, 
                "w":  1.02200 
            }
        "#;
        let c1: CubicPureRecord = from_str(data1).unwrap();
        let c2: CubicPureRecord = from_str(data2).unwrap();
        let pr1 = PureRecord::new(0.0, "".into(), c1);
        let pr2 = PureRecord::new(0.0, "".into(), c2);        

        let model = PR78.into();
        let params = CubicParameters::new(vec![pr1,pr2],vec![], model);
        // println!("{}",params);
        let s = params.to_string();

        
        println!("{}", s);

        let cub = Cubic::from_parameters(params);
        let t = 298.15;
        let rho = 1000.0;
        let x = array![0.5,0.5];
        let test = cub.all_derived_properties(t, rho, &x);
        // let json = serde_json::to_string_pretty(&r_properties).unwrap();
        // println!("{}",json);
        let data = r#"
            {
              "a": -6.739499683331773,
              "dadv": -5538.109023777433,
              "dadni": [
                -4.772941005259919,
                -19.78227640895849
              ]
            }
            "#;

        let answer:ResidualDerivedProperties = from_str(data).unwrap();

        answer.comparison(test,None);
        
    }

    #[test]
    fn test3(){
        let data1 = r#"
            {   
                "tc":507.317,
                "pc":3290000,
                "w":0.275000
            }
        "#;
        let data2 = r#"
            {   
                "tc": 860.372, 
                "pc": 1246000, 
                "w":  1.02200 
            }
        "#;

        let c1: CubicPureRecord = from_str(data1).unwrap();
        let c2: CubicPureRecord = from_str(data2).unwrap();
        let pr1 = PureRecord::new(0.0, "".into(), c1);
        let pr2 = PureRecord::new(0.0, "".into(), c2);        

        let model = PR76.into();
        let params = super::CubicParameters::new(vec![pr1,pr2],vec![], model);
        let s = params.to_string();

        
        println!("{}", s);

        let cub = Cubic::from_parameters(params);
        let t = 298.15;
        let rho = 1000.0;
        let x = array![0.5,0.5];
        let test = cub.all_derived_properties(t, rho, &x);
        // let json = serde_json::to_string_pretty(&r_properties).unwrap();
        // println!("{}",json);
        let data = r#"
            {
              "a": -6.552525570828891,
              "dadv": -5381.453914153913,
              "dadni": [
                -4.708714861213096,
                -19.15924410875251
              ]
            }
            "#;

        let answer:ResidualDerivedProperties = from_str(data).unwrap();

        answer.comparison(test,None);

    }
    #[test]
    fn test4(){

        let data1 = r#"
            {   
                "a0":2.47272329e0,
                "b":9.97464159e-5,
                "c1":7.78348800e-1,
                "tc":5.073170e2
            }
        "#;
        let data2 = r#"
            {   
                "a0":1.87787673e1,
                "b":4.46665087e-4,
                "c1":1.66890260e0,
                "tc":8.603720e2
            }
        "#;

        let c1: CubicPureRecord = from_str(data1).unwrap();
        let c2: CubicPureRecord = from_str(data2).unwrap();
        let pr1 = PureRecord::new(0.0, "".into(), c1);
        let pr2 = PureRecord::new(0.0, "".into(), c2);        
        let model = PR76.into();
        let params = CubicParameters::new(vec![pr1,pr2],vec![], model);
        let s = params.to_string();

        
        println!("{}", s);

        let cub = Cubic::from_parameters(params);
        let t = 298.15;
        let rho = 1000.0;
        let x = array![0.5,0.5];
        let test = cub.all_derived_properties(t, rho, &x);
        // let json = serde_json::to_string_pretty(&r_properties).unwrap();
        // println!("{}",json);
        let data = r#"
            {
              "a": -6.552525570828891,
              "dadv": -5381.453914153913,
              "dadni": [
                -4.708714861213096,
                -19.15924410875251
              ]
            }
            "#;

        let answer:ResidualDerivedProperties = from_str(data).unwrap();

        answer.comparison(test,Some(1e-5));

    }

    #[test]
    fn test_binary_from_json(){

        let data1 = r#"
        {   
                "name": "first",
                "a0":2.47272329e0,
                "b":9.97464159e-5,
                "c1":7.78348800e-1,
                "tc":5.073170e2
            }
        "#;
        let data2 = r#"
            {   
                "name": "second",
                "a0":1.87787673e1,
                "b":4.46665087e-4,
                "c1":1.66890260e0,
                "tc":8.603720e2
            }
        "#;

        let data3 = r#"
            {   
                "name": "third",
                "a0":0.87787673e1,
                "b": 0.46665087e-4,
                "c1":0.66890260e0,
                "tc":3.603720e2
            }
        "#;

        let bin = r#"
            {   
                "id1": "first",
                "id2": "second",
                "kij": 0.12

            }
        "#;
        let pr1: PureRecord<CubicPureRecord> = from_str(data1).unwrap();
        let pr2: PureRecord<CubicPureRecord> = from_str(data2).unwrap();
        let pr3: PureRecord<CubicPureRecord> = from_str(data3).unwrap();
        let bin: BinaryRecord<CubicBinaryRecord> = from_str(bin).unwrap();
        let model = PR76.into();
        let params = CubicParameters::new(vec![pr1,pr2,pr3],vec![bin], model);
        assert_eq!(params.aij[(0,1)], 0.12);
        assert_eq!(params.aij[(0,1)], 0.12);
        let s = params.to_string();

        
        println!("{}", s);

    }

}

    // // left  = -6.552525570828891
    // right = -6.552527221312991