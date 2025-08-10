// mod's
pub mod eos;
pub mod density_solver;




use ndarray::{Array1, Array2};
use std::fmt::{Display, Write};
use std::sync::Arc;

use crate::models::cpa::CPA;
use crate::models::IDEAL_GAS_CONST;
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult, EquationOfState};
use crate::state::density_solver::{DensityInitialization,density};
pub type E<R>=EquationOfState<R>;
pub type S<R>=State<R>;

pub type StateResult<R>=Result<S<R>,EosError>;


#[derive(Clone)]
pub struct State<R:Residual>{

    pub eos:Arc<EquationOfState<R>>,
    // Kelvin
    pub t:f64,
    pub p:f64,
    pub rho:f64,
    pub x:Array1<f64>,
    // pub cache: Option<Cache>
}
pub fn fmt_array2<W: Write>(
    f: &mut W,
    name: &str,
    mat: &Array2<f64>,
    decimals: usize,
) -> std::fmt::Result {
    writeln!(f, "{} = [", name)?;
    for row in mat.rows() {
        write!(f, "    [")?;
        for (j, val) in row.iter().enumerate() {
            if j > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{:.*}", decimals, val)?;
        }
        writeln!(f, "],")?;
    }
    writeln!(f, "]")?;
    Ok(())
}
impl Display for State<CPA> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        
        let assoc=&self.eos.residual.assoc;
        let t=self.t;
        let p=self.p;
        let rho=self.rho;
        let x=&self.x;
        writeln!(f,"State Variables:")?;
        writeln!(f,"Temperature    = {t}   K")?;
        writeln!(f,"Pressure       = {p}   Pa")?;
        writeln!(f,"Density        = {rho} mol/m3")?;
        writeln!(f,"Molar Fraction = {x}")?;

        //Association
        let gmix=assoc.g_func(rho, x);
        let dlng_drho=assoc.dlngdrho(rho, x);
        let dlng_dni=assoc.dlngdni(rho, x);

        let frac_non_assoc_sites=assoc.X_tan(t, rho, x).unwrap_or_default();
        let delta=assoc.delta(t, gmix);
        let mu_assoc=assoc.residual_chemical_potential(t, rho, x).unwrap_or_default();

        let m=assoc.get_m(x);
        let m_as_matrix=m.to_shape((m.len(),1)).unwrap_or_default();

        let m_mul_by_m_transp=m_as_matrix.dot(&m_as_matrix.t());
        let h=assoc.h_func(x, &frac_non_assoc_sites);
        
        let k=&delta*m_mul_by_m_transp*&assoc.parameters.pmat;

        writeln!(f,"Association Parameters:")?;

        let assocp_as_string=format!("{}",assoc.parameters);

        writeln!(f,"{}",assocp_as_string)?;

        writeln!(f,"Association Residual Properties:")?;


        writeln!(f, "h/2        = {}", h*0.5)?;
        writeln!(f, "dln(g)/dρ  = {}", dlng_drho)?;
        writeln!(f, "dln(g)/dni = {}", dlng_dni)?;
        writeln!(f, "m          = {}", m)?;
        writeln!(f, "X          = {}", frac_non_assoc_sites)?;
        writeln!(f, "μ_assoc    = {}", mu_assoc)?;
        fmt_array2(f, "Δ matrix", &delta, 10)?;
        fmt_array2(f, "K matrix", &k, 10)?;

        Ok(())
        
    }
}
impl<R:Residual> State<R> {

    pub fn new_trx(eos:&Arc<E<R>>,t:f64,rho:f64,x:Array1<f64>)->StateResult<R>{

        let p =eos.pressure(t, rho, &x);

        match p {
            Ok(pres)=>Ok(
                Self{
                    eos:Arc::clone(&eos),
                    t,
                    p:pres,
                    rho,
                    x,
                }
            ),
            Err(e)=>Err(e)
        }

    }
    pub fn new_tpx(eos:&Arc<E<R>>,t:f64,p:f64,x:Array1<f64>,phase:DensityInitialization)->StateResult<R>{

        match phase {
            DensityInitialization::Liquid=>{

                let guess=0.99;
                // let guess=0.9;
                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }

            }

            DensityInitialization::Vapor=>{

                let bm=eos.residual.bmix(&x);
                let guess = bm/(bm + (IDEAL_GAS_CONST*t)/p);

                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }
            }

            DensityInitialization::Guess(old_density)=>{
                let bm=eos.residual.bmix(&x);
                let guess=old_density*bm;

                density(eos, t, p, x, guess)
                // match eos.density_solver(t, p, &x, guess){

                //     Ok(rho)=>{State::new_trx(eos, t, rho, x)}
                //     Err(e)=>{Err(e)}
                // }
            }
            
        }
        
    }




}


impl<R:Residual> State<R> {
    
    pub fn lnphi(&self)->EosResult<Array1<f64>>{
        self.eos.lnphi(self.t, self.rho, &self.x)
    }
    pub fn pressure(&self)->EosResult<f64>{
        self.eos.pressure(self.t, self.rho, &self.x)
    }
    pub fn bmix(&self)->f64{
        self.eos.residual.bmix(&self.x)
    }
    pub fn composition(&self)->Array1<f64>{
        self.x.clone()
    }
    pub fn temperature(&self)->f64{
        self.t
    }
}
