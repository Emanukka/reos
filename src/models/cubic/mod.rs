use crate::models::cubic::parameters::CubicParameters;
use crate::state::eos::EosResult;
use crate::residual::Residual;
use ndarray::Array1;
use ndarray::Array2;
use crate::models::IDEAL_GAS_CONST as R;

pub mod parameters;
#[derive(Clone)]
pub struct Cubic{
    pub parameters:CubicParameters,
    pub model:CubicModel,
}

impl Cubic {
    pub fn new(p:CubicParameters,model:CubicModel)->Self{
        Self { parameters: p, model }
    }

    pub fn set_binary(&mut self,i:usize,j:usize,kij_a:Option<f64>,kij_b:f64){
        self.parameters.pbinary.insert((i,j), (kij_a.unwrap_or(0.0),kij_b));
    } }

#[derive(Clone)]
pub enum CubicModel{
    SRK,
    PR
}
// Models (only SRK and PR)
impl Cubic{
    fn eps(&self)->f64{
        match self.model {
            CubicModel::SRK=>{0.0}
            CubicModel::PR=>{1.0-2.0_f64.sqrt()}
        }
    }
    fn sig(&self)->f64{
        match self.model {
            CubicModel::SRK=>{1.0}
            CubicModel::PR=>{1.0+2.0_f64.sqrt()}
        }    
    }
}
impl Cubic {
    fn calc_sqrt_aij_matrix(
        &self,t:f64
    )->Array2<f64>{
        
        let p=&self.parameters;
        let vtr =  t/&p.vtc;

        let alpha =  (1.0+&p.vkappa*(1.0 - vtr.sqrt())).pow2();

        let va = &p.va0*alpha;
        
        // println!("a({} K)={}",t,&va);
        let n = p.va0.len();
        let mut aij = Array2::<f64>::zeros((n,n));
        for i in 0..n {
            for j in i..n{

                let ij=&(i,j);
                let ji=&(j,i);

                let binary_parameters= p.pbinary.get(ij).or_else(||p.pbinary.get(ji));

                // if binary not contains (a,b), then kij=0.0

                let mut kij=0.0;
                match binary_parameters{
                Some(tup)=>{
                    let (a,b)=(tup.0,tup.1);
                    kij=a*t+b;
                }
                None=>{}
                }
                let prod = va[i]*va[j];
                let sqrt_aij = prod.sqrt();
                aij[(i,j)] = sqrt_aij*(1.0-kij); 
                aij[(j,i)] = aij[(i,j)]; 
            }
        }

        aij
    }
    fn calc_amix(&self, t: f64, vx:&Array1<f64>)->f64{

        // a = x*Aij*x
        let aij_mat = self.calc_sqrt_aij_matrix(t);
        
        let aij_mul_x: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>> = aij_mat.dot(vx);

        let a = vx.dot(&aij_mul_x);
        a
    }


    
}

impl Residual for Cubic {

    fn components(&self)->usize {
        self.parameters.ncomp
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.parameters.vb.dot(x)
    }

    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {
        let bm=self.bmix(x);
        let am=self.calc_amix(t, x);
        let vm=1.0/rho;
        let delta= (vm+self.sig()*bm)*(vm+self.eps()*bm);

        let repulsive= (R*t*bm)/(vm*(vm-bm));
        let attractive = am/delta;

        Ok(
            repulsive-attractive
        )

    }

    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>> {
        let bm=self.bmix(x);
        let am = self.calc_amix(t, x);
        let dbdni=&self.parameters.vb;
        let vm=1./rho;
        let ln=(vm/(vm-bm)).ln();
        let q=am/(bm*R*t);

        let sig=self.sig();
        let eps=self.eps();
        let maij= self.calc_sqrt_aij_matrix(t);
        let aij_dot_vx:Array1<f64> = maij.dot(x);

        let dadni=2.0*aij_dot_vx - am;
        
        let dqni=q*(1.0 + dadni/am - dbdni/bm);
        let i: f64 = (1.0/(sig-eps))*f64::ln((vm + sig*bm )/(vm + eps*bm));
        let z_residual=(self.pressure(t,rho,x)?)/(rho*R*t);

        Ok(
        (dbdni/bm)*z_residual +ln - dqni*i
        )
        

    }
    
}

