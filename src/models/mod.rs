use std::ops::Mul;

pub mod associative;
pub mod cubic;
pub mod cpa;


// Association
pub const A: usize = 0;
pub const B: usize = 1;
pub const C: usize = 2;
pub const NS:  usize = 3;
pub const IDEAL_GAS_CONST: f64 = 8.31446261815324;

// pub const SINT: [(usize,usize);4]=[(A,B),(A,C),(B,C),(C,C)];

#[derive(Debug, Clone, PartialEq)]
pub enum Site{
    A(f64),
    B(f64),
    C(f64),
    Null
}

impl Site {
    
    pub fn n(&self)->f64{
        match self {
            Site::A(m)=>*m,
            Site::B(m)=>*m,
            Site::C(m)=>*m,
            Site::Null=>0.0
        }
    }
}

impl Default for Site {
    
    fn default() -> Self {
        Site::Null
    }
}
pub const SITES:[usize;3]=[A,B,C];

// impl Mul {
    
// }
impl<'a, 'b> Mul<&'b Site> for &'a Site {
    type Output = f64;

    fn mul(self, rhs: &'b Site) -> f64 {

        if (self.n()==0.0) || (rhs.n()==0.0){
            0.0
        }else{
            match (self,rhs) {
            (Site::A(_),Site::A(_))=>{0.0}
            (Site::B(_),Site::B(_))=>{0.0}
            (_,_)=>1.0
            }
        }
    }
}

// impl<'a> From<&'a Site> for usize {
//     fn from(s: &'a Site) -> Self {
//         match s {

//             Site::A=>0,
//             Site::B=>1,
//             Site::C=>2,

//         }
//     }
// }

