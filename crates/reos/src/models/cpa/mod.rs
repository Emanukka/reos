pub mod parameters;
pub mod rdf;
pub mod association;
pub mod tests;
use crate::models::IDEAL_GAS_CONST;
use crate::models::associative::parameters::CombiningRule;
use crate::models::cpa::parameters::CPAParameters;
use crate::models::cpa::rdf::{CS, Kontogeorgis, RdfModel};
use crate::models::cubic::{CubicModel, SRK};


use crate::models::{cubic::Cubic,cpa::association::Associative};
use crate::residual::Residual;
use crate::state::eos::EosResult;
use crate::state::State;
use core::{f64, str};
use ndarray::{Array1, Array2};


pub struct CR1;
impl CombiningRule for CR1 {
    
    fn to_string()->String {
        "CR1".to_string()
    }
    fn which()->Self where  Self: Sized {
        CR1
    }
    fn association_strength_jl(
            t:f64,
            g_contact:f64,
            bij:&Array2<f64>,
            epsilon:&Array2<f64>,
            beta:&Array2<f64>,
        )->Array2<f64> {
        
        
    let x = epsilon / (IDEAL_GAS_CONST * t);
    // dbg!(&x);

    let exp_term = x.mapv(|v| v.exp());
    // dbg!(&exp_term);

    let expm1_term = &exp_term - 1.0;
    // dbg!(&expm1_term);

    let result = bij * beta * g_contact * expm1_term;
    // dbg!(&result);

    result

    }
}
pub struct ECR;
impl CombiningRule for ECR {
    
    fn to_string()->String {
        "ECR".to_string()
    }
    fn which()->Self where  Self: Sized {
        ECR
    }
    fn association_strength_jl(
            t:f64,
            g_contact:f64,
            bij:&Array2<f64>,
            epsilon:&Array2<f64>,
            beta:&Array2<f64>,
        )->Array2<f64> {
        
        let ns = epsilon.ncols(); 
        let cr1 = CR1::association_strength_jl(t, g_contact, bij, epsilon, beta);
        let binding = cr1.diag();
        let diag=binding.to_shape((ns,1)).unwrap();
        let dd=diag.dot(&diag.t());

        dd.sqrt()

    }
}



// #[derive(Clone)]
pub struct CPA <C:CubicModel,T:CombiningRule,R:RdfModel>{
    pub cubic: Cubic<C>,
    pub assoc: Associative<T,R>,
}

impl<C:CubicModel,T:CombiningRule,R:RdfModel> CPA<C,T,R> {

    pub fn from_parameters(parameters:CPAParameters<C>)->Self{
        let cubic=Cubic::from_parameters(parameters.cubic);
        let l_t = &parameters.assoc.lambda_t;
        let g_t = &parameters.assoc.gamma_t;
        let b_components = &cubic.parameters.b;
        let b_sites = l_t.dot(&g_t.dot(b_components));
        let bij = 0.5 * (&b_sites + &b_sites.t());
        let rdf = R::new(b_components.flatten().to_owned(),bij);
        let assoc: Associative<T,R> = Associative::from_parameters(parameters.assoc,rdf);
        
        // dbg!(&assoc.parameters.p);
        // dbg!(&assoc.parameters.epsilon);
        // dbg!(&assoc.parameters.beta);
        // dbg!(&assoc.parameters.lambda);
        // dbg!(parameters.assoc.lambda_t);
        // dbg!(parameters.assoc.lambda_;);





        Self
        {
            cubic,
            assoc,
        }
    }
    // pub fn from_models(c:Cubic<C>,a:Associative<R>)->Self{
    //     Self{cubic:c,assoc:a}
    // }


}

impl<C:CubicModel,T:CombiningRule,R:RdfModel> Residual for CPA<C,T,R> {
    
    fn components(&self)->usize {
        self.cubic.components()
    }
    fn r_pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        Ok(
        self.cubic.r_pressure(t, rho, x)?
        +self.assoc.r_pressure(t, rho, x)?
        )

    }
    fn r_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        Ok(
            self.cubic.r_chemical_potential(t, rho, x)?
            +self.assoc.r_chemical_potential(t, rho, x)?
        )
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.cubic.bmix(x)
    }

    fn r_helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        Ok(
        self.cubic.r_helmholtz(t, rho, x)?+self.assoc.r_helmholtz(t, rho, x)?)
        
    }
}


impl<C:CubicModel,T:CombiningRule,R:RdfModel> State<CPA<C,T,R>> {
    
    pub fn non_bonded_sites(&self)->Array1<f64> {

        let (t,rho,x)=(self.t,self.rho,&self.x);
        self.eos.residual.assoc.x_tan(t, rho, x).unwrap()
    }

}


pub type SCPAsrkCR1 = CPA<SRK,CR1,Kontogeorgis>;
// pub type CPAsrkCR1 = CPA<SRK,CR1,CS>;
pub type SCPAsrkECR = CPA<SRK,ECR,Kontogeorgis>;
// pub type SCPAsrkECR = CPA<SRK,ECR,Kontogeorgis>;


