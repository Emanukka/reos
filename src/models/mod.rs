use core::fmt;

pub mod associative;
pub mod cubic;
pub mod cpa;


// Association
pub const A: usize = 0;
pub const B: usize = 1;
pub const C: usize = 2;
pub const NS:  usize = 3;
pub const IDEAL_GAS_CONST: f64 = 8.31446261815324;


#[derive(Debug, Clone, PartialEq)]
pub enum Site{
    A(usize),
    B(usize),
    C(usize),
}

impl Site {
    
    ///Get component of site alpha
    pub fn i(&self)->usize{
        match self {
            Self::A(i)=>*i,
            Self::B(i)=>*i,
            Self::C(i)=>*i,
        }
    }
    ///Get type index site alpha
    pub fn t(&self)->usize{
        match self {
            Self::A(_)=>A,
            Self::B(_)=>B,
            Self::C(_)=>C,
        }
    }
}

pub const SITES:[usize;3]=[A,B,C];

impl fmt::Display for Site {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Site::A(i) => write!(f, "A({})", i),
            Site::B(i) => write!(f, "B({})", i),
            Site::C(i) => write!(f, "C({})", i),
        }
    }
}
