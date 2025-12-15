use ndarray::{Array1, Array2};
use serde::{Deserialize, Serialize, ser::SerializeStruct};
use crate::{models::cubic::{CubicModel}, parameters::{Parameters, PureRecord}};


#[derive(Debug)]
pub struct CubicParameters<T:CubicModel>{
    pub ncomp:usize,
    pub a0:   Array2<f64>,
    pub b:    Array2<f64>,
    pub tc:   Array2<f64>,
    pub kappa:Array2<f64>,
    pub aij:  Array2<f64>, //kij = Aij + Bij*T
    pub bij:  Array2<f64>,
    // pub w: Option<Array<f64>>
    pub model: T

}

#[derive(Clone,Debug,Serialize, Deserialize)]
#[serde(untagged)]
pub enum CubicPureRecord{
    Set1{a0:f64,b:f64,kappa:f64,tc:f64},
    Set2{tc:f64,pc:f64,w:f64},

}



impl CubicPureRecord {
    
    pub fn new_set1(a0:f64,b:f64,kappa:f64,tc:f64)->Self{
        Self::Set1{
            a0,
            b,
            kappa,
            tc
        }
    }
    pub fn new_set2(tc:f64,pc:f64,w:f64)->Self{
        Self::Set2{
            tc,
            pc,
            w,
        }
    }
}

impl PureRecord for CubicPureRecord {
    
}

impl<T:CubicModel> CubicParameters<T>{
    pub fn new(model:T,records:Vec<CubicPureRecord>)->Self {
        let ncomp =records.len();

        let mut ma0:Array2<f64>    = Array2::<f64>::zeros((ncomp,1));
        let mut mb:Array2<f64>     = Array2::<f64>::zeros((ncomp,1));
        let mut mtc:Array2<f64>    = Array2::<f64>::zeros((ncomp,1));
        let mut mkappa:Array2<f64> = Array2::<f64>::zeros((ncomp,1));
        let aij:Array2<f64>       = Array2::<f64>::zeros((ncomp,ncomp));
        let bij:Array2<f64>       = Array2::<f64>::zeros((ncomp,ncomp));

        for (i,record) in records.iter().enumerate(){

            match record {
                
                &CubicPureRecord::Set1 { a0, b, kappa, tc } =>
                {
                    ma0[(i,0)] = a0;
                    mb[(i,0)] = b;
                    mkappa[(i,0)] = kappa;
                    mtc[(i,0)] = tc;
                }
                
                &CubicPureRecord::Set2 { tc, pc, w } =>
                {
                    ma0[(i,0)] = model.acrit(tc, pc);
                    mb[(i,0)] = model.bcrit(tc, pc);
                    mkappa[(i,0)] = model.kappa_from_w(w);
                    mtc[(i,0)] = tc;
                }
            }
        }    
        Self{
            ncomp,
            a0: ma0,
            b: mb,
            kappa: mkappa,
            tc: mtc,
            aij,
            bij,
            model
        } 
    }

    pub fn set_kij(&mut self,i:usize,j:usize,kij:f64){

        self.aij[(i,j)] = kij;
        self.aij[(j,i)] = kij;
    }
    pub fn set_kij_temperature_dependent(&mut self,i:usize,j:usize,aij:f64,bij:f64){

        self.aij[(i,j)] = aij;
        self.aij[(j,i)] = aij; 
 
        self.bij[(i,j)] = bij;
        self.bij[(j,i)] = bij;
    }

}

impl<C:CubicModel> Parameters<CubicPureRecord> for CubicParameters<C> {
    fn from_records(records:Vec<CubicPureRecord>) -> Self {
        
        Self::new(C::model(), records)
    }
}

impl<T:CubicModel> Serialize for CubicParameters<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: serde::Serializer {
        
        let mut state = serializer.serialize_struct("CubicParameters", 5)?;
        state.serialize_field("cubic_model", &self.model.which())?;

        state.serialize_field("a0", &self.a0
        .flatten()
        .iter()
        .map(|v| format!("{:.9e}", v))
        .collect::<Vec<_>>())?;
        
        state.serialize_field("b", &self.b
        .flatten()
        .iter()
        .map(|v| format!("{:.9e}", v))
        .collect::<Vec<_>>())?;

        state.serialize_field("kappa", &self.kappa
        .flatten()
        .iter()
        .map(|v| format!("{:.9e}", v))
        .collect::<Vec<_>>())?;
        
        state.serialize_field("tc", &self.tc
        .flatten()
        .iter()
        .map(|v| format!("{:.9e}", v))
        .collect::<Vec<_>>())?;

        state.end()        
    }
}

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

    use crate::{models::cubic::{Cubic, PR76, PR78, parameters::CubicPureRecord}, parameters::Parameters, residual::{Residual, ResidualDerivedProperties}};


    #[test]
    fn test1(){
        let tc = array![ 507.317  , 860.372];
        let pc = array![ 3290000. , 1246000.];
        let w  = array![ 0.27500  , 1.022  ];

        let c1 = CubicPureRecord::new_set2(tc[0], pc[0], w[0]);
        let c2 = CubicPureRecord::new_set2(tc[1], pc[1], w[1]);
        // let params = super::CubicParameters::from_tpw(PR78, tc, pc, w);
        let params = super::CubicParameters::<PR78>::from_records(vec![c1,c2]);
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
        // let params = super::CubicParameters::from_tpw(PR78, tc, pc, w);
        let params = super::CubicParameters::<PR78>::from_records(vec![c1,c2]);
        // println!("{}",params);
        let json = serde_json::to_string_pretty(&params).unwrap();

        
        println!("{}", json);

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
        let params = super::CubicParameters::<PR76>::from_records(vec![c1,c2]);

        let json = serde_json::to_string_pretty(&params).unwrap();

        
        println!("{}", json);

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
                "kappa":7.78348800e-1,
                "tc":5.073170e2
            }
        "#;
        let data2 = r#"
            {   
                "a0":1.87787673e1,
                "b":4.46665087e-4,
                "kappa":1.66890260e0,
                "tc":8.603720e2
            }
        "#;

        let c1: CubicPureRecord = from_str(data1).unwrap();
        let c2: CubicPureRecord = from_str(data2).unwrap();
        let params = super::CubicParameters::<PR76>::from_records(vec![c1,c2]);

        let json = serde_json::to_string_pretty(&params).unwrap();

        
        println!("{}", json);

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

}

    // // left  = -6.552525570828891
    // right = -6.552527221312991