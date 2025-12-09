use core::fmt;

use ndarray::Array1;


// Association
pub const A: usize = 0;
pub const B: usize = 1;
pub const C: usize = 2;
pub const NS:  usize = 3;


#[derive(Debug, Clone, PartialEq)]
pub enum Site{
    A(usize),
    B(usize),
    C(usize),
}

impl Site {
    
    ///Get component of site alpha
    pub fn c(&self)->usize{
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
            Site::A(i) => write!(f, "A(owner={})", i),
            Site::B(i) => write!(f, "B(owner={})", i),
            Site::C(i) => write!(f, "C(owner={})", i),
        }
    }
}
impl From<Site> for (usize,usize) {

    fn from(value: Site) -> Self {
        
        (value.t(),value.c())
    }
}

#[cfg(test)]
pub mod tests{
    use ndarray::array;

    use super::Site;


    #[test]
    fn print_sites(){

        let s1=Site::A(0);
        let s2=Site::B(0);

        println!("Site1 = {}",s1);
        println!("Site2 = {}",s2);


        let sites=array![s1,s2];

        println!("Sites = {}",sites)
        
    }
}