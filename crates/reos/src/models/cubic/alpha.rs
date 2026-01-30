use std::clone;

use enum_dispatch::enum_dispatch;


pub fn alpha_soave(t:f64, tc:f64, kappa:f64) -> f64 {
    
    let tr = t /tc;
    (1.0 + kappa * (1.0 - tr.sqrt())).powi(2)

}

pub fn dalpha_dt_soave(t:f64, tc:f64, kappa:f64) -> f64 {
    // let a0 = p.a0[i];
    let tr = t / tc;
    let alpha = alpha_soave(t, tc, kappa);
    - kappa * (alpha * tr).sqrt() / t
}

#[derive(Clone, Debug)]
pub enum AlphaParameter {
    Soave(f64),
    Twu91{l:f64, n:f64, m:f64}
}


#[derive(Clone, Debug)]
pub struct Alpha {
    pub parameters: Vec<AlphaParameter>,
}

impl std::fmt::Display for Alpha {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}",
            self.parameters
        )
    }
}

impl std::fmt::Display for AlphaParameter {



    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        match self {
            Self::Soave(kappa) => write!(f, "Soave(kappa={})",kappa),
            Self::Twu91 { l, n, m } => write!(f, "Twu91(l={},n={},m={})",l,n,m)
        }

    }

}

impl Alpha {
    
    pub fn new(parameters:Vec<AlphaParameter>) -> Self {
        Self { parameters }
    }
    pub fn dalpha_dt(&self, i:usize, t:f64, tc:f64) -> f64 {
        
        let tr = t / tc;

        match self.parameters[i] {
        
            AlphaParameter::Soave(kappa) => {
                
                let alpha = self.alpha(i, tr);

                - kappa * (alpha * tr).sqrt() / t
                
            }

            AlphaParameter::Twu91 { l, n, m } => {

                // - (n * tr.powf((m - 1.) * n) * (l * m * tr.powf(m * n) - m + 1.) * (l * (1. - tr.powf(m * n))).exp()) / t

                let alpha = self.alpha(i, tr);
                alpha * n * (m - 1. - l * m * tr.powf(n * m)) / t

            }
        }

    }
    pub fn alpha(&self,i:usize, tr:f64) -> f64 {
        

        match self.parameters[i] {
        
            AlphaParameter::Soave(kappa) => {
                
                (1.0 + kappa * (1.0 - tr.sqrt())).powi(2)
                
            }

            AlphaParameter::Twu91 { l, n, m } => {

                tr.powf(n * (m - 1.)) * (l * (1. - tr.powf(n * m))).exp()

            }
        }
    }
}

mod tests {

}
