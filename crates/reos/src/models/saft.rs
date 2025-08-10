
use std::f64::consts::PI;

use nalgebra::coordinates::X;
use ndarray::arr1;
use ndarray::arr2;
use ndarray::array;
use ndarray::Array1;
use ndarray::Array2;
use ndarray::ArrayView1;
use ndarray::ArrayView2;
use crate::eos::analytic_diff::associative::Associative;

use super::associative::ASCParameters;
use super::cubic::CubicParameters;
use super::density_solver::bissec;
use super::density_solver::Bissection;
use super::equation_of_state::EosError;
// use super::cubic::SRKeos;
// use super::associative::ASCeos;
pub use super::equation_of_state::EquationOfState;
// 1:Gross, J., & Sadowski, G. (2001). 
// Perturbed-Chain SAFT: An Equation of State 
// Based on a Perturbation Theory for Chain Molecules
// Joachim. Industrial and Engineering Chemistry Research, 40, 1244–1260.
//  https://doi.org/10.1021/ie0003887
    
// Gross, J., & Sadowski, G. (2019). Reply to Comment on 
// “perturbed-Chain SAFT: An Equation of State Based on 
// a Perturbation Theory for Chain Molecules.” Industrial 
// and Engineering Chemistry Research, 58(14), 5744–
// 5745. https://doi.org/10.1021/acs.iecr.9b01515

// Michelsen, M. L., & Hendriks, E. M. (2001). 
// Physical properties from association models. 
// Fluid Phase Equilibria, 180(1–2), 165–174. 
// https://doi.org/10.1016/S0378-3812(01)00344-2

// X association algorithm:
// Tan, S. P., Adidharma, H., & Radosz, M. (2004).
// Generalized Procedure for Estimating the Fractions of 
// Nonbonded Associating Molecules and Their Derivatives in 
// Thermodynamic Perturbation Theory. Industrial and Engineering 
// Chemistry Research, 43(1), 203–208. https://doi.org/10.1021/ie034041q

// Delta funtion:
// Tan, S. P., Adidharma, H., & Radosz, M. (2004). Generalized Procedure 
// for Estimating the Fractions of Nonbonded Associating Molecules and 
// Their Derivatives in Thermodynamic Perturbation Theory. Industrial 
// and Engineering Chemistry Research, 43(1), 203–208. 
// https://doi.org/10.1021/ie034041q

const AP:[[f64;7]; 3]= [
[9.10563145e-01, 
 6.36128145e-01, 
 2.68613479e+00, 
-2.65473625e+01, 
 9.77592088e+01, 
-1.59591541e+02, 
 9.12977741e+01],

 [-3.08401692e-01,
 1.86053116e-01,
-2.50300473e+00,
 2.14197936e+01,
-6.52558853e+01,
 8.33186805e+01,
-3.37469229e+01],

 [-9.06148351e-02,
 4.52784281e-01,
 5.96270073e-01,
-1.72418291e+00,
-4.13021125e+00,
 1.37766319e+01,
-8.67284704e+00],
 ];

const BP:[[f64;7]; 3]= [
[7.24094694e-01,        
 2.23827919e+00,        
-4.00258495e+00,        
-2.10035768e+01,        
 2.68556414e+01,        
 2.06551338e+02,        
-50.8003365888685e0 * 7.],

 [-5.75549808e-01,       
 6.99509552e-01,       
 3.89256734e+00,       
-1.72154716e+01,       
 1.92672264e+02,       
-1.61826462e+02,       
-23.6010990650801e0 *7.],

 [ 9.76883116e-02,
-2.55757498e-01,
-9.15585615e+00,
 2.06420760e+01,
-3.88044301e+01,
 9.36267741e+01,
-4.23812936930675e0 * 7.],
 ];



// const AB:[f64; 21]= [
//  7.24094694e-01,          -5.75549808e-01,          9.76883116e-02,
//  2.23827919e+00,           6.99509552e-01,         -2.55757498e-01,
// -4.00258495e+00,           3.89256734e+00,         -9.15585615e+00,
// -2.10035768e+01,          -1.72154716e+01,          2.06420760e+01,
//  2.68556414e+01,           1.92672264e+02,         -3.88044301e+01,
//  2.06551338e+02,          -1.61826462e+02,          9.36267741e+01,
// -50.8003365888685e0 * 7., -23.6010990650801e0 *7., -4.23812936930675e0 * 7.];

        

pub fn navo()->f64{ 6.02214076 * f64::powf(10., 23.)}
// eq. 21( molecules/volume)
pub fn molecules_total_density(rho:f64)->f64{rho*navo()/1.0e30}
pub struct PCSaft{
    assoc:Associative,
    hs:HardShere,
    disp:Dispersion,

}

pub struct HardShere{

}

// (dens, T, x,self.ap,self.bp,self.ncomp,self.sigma,self.epsilon_k,self.m ,self.kbi)
pub struct Dispersion{

    pub m:Array1<f64>,
    /// (dim:Angstromn)
    pub s:Array1<f64>,
    /// epsilon dividid by Boltzmann Const (dim:Kelvin)
    pub e:Array1<f64>,
    /// Combination eps
    pub eij:Array2<f64>,
    /// Combination sigma
    pub sij:Array2<f64>,
    /// kij matrix
    pub kbi:Array2<f64>,

}

impl Dispersion {

    // def pc_a_e_b(x,ap,bp,m):

    // mmed = pc_mmed(x,m)

    // a = ap[:, 0] + (mmed-1)*ap[:, 1]/mmed + \
    //     (1-1/mmed)*(1-2/mmed)*ap[:, 2]
    // b= bp[:, 0] + (mmed-1)*bp[:, 1]/mmed + \
    //     (1-1/mmed)*(1-2/mmed)*bp[:, 2]

    // return a, b

    
    pub fn new(m:Array1<f64>,s:Array1<f64>,e:Array1<f64>,kbi:Array2<f64>)->Self{


        assert!(m.len()==s.len());
        assert!(m.len()==e.len());
        assert!((kbi.len()==e.len()*e.len()));
        let n=m.len();

        let sij= 
        Array2::from_shape_fn((n,n),
        |(i,j)|0.5*(s[i]+s[j]) );

        let eij= 
        Array2::from_shape_fn((n,n),
        |(i,j)| (e[i]*e[j]).sqrt()*(1.-kbi[(i,j)]) );

        Dispersion{
            m,s,e,eij,sij,kbi
        }
    }
    pub fn mmed(&self,x:&Array1<f64>)->f64 {self.m.dot(x) }

    pub fn a_and_b(&self,x:&Array1<f64>)->(Array1<f64>,Array1<f64>){

        let mmed=self.mmed(x);
        
        // Ok (AP e AB são arrays fixos na memória, de modo que basta 
        // utilizar referências aos seus dados, evitando alocação desnecessária
        // de vetores ou matriz de 21 elementos)

        // let ap = ArrayView2::from_shape((7, 3), &AP).unwrap();
        // let bp = ArrayView2::from_shape((7, 3), &AB).unwrap();

        let coef1=(mmed-1.)/mmed;
        let coef2=(1.-1./mmed)*(1.-2./mmed);

        let ap0=&ArrayView1::from_shape(7, &AP[0]).unwrap();
        let ap1=&ArrayView1::from_shape(7, &AP[1]).unwrap();
        let ap2=&ArrayView1::from_shape(7, &AP[2]).unwrap();

        let bp0=&ArrayView1::from_shape(7, &BP[0]).unwrap();
        let bp1=&ArrayView1::from_shape(7, &BP[1]).unwrap();
        let bp2=&ArrayView1::from_shape(7, &BP[2]).unwrap();
        let a = ap0+coef1*ap1+coef2*ap2;
        let b = bp0+coef1*bp1+coef2*bp2;
        (a,b)
    } 


    // # A.16 and A.17 ok!
    pub fn int_1(&self,rho:f64,t:f64,x:&Array1<f64>)->f64{
        let eta3=self.csi(rho, t, x)[3];
        let (a,b)=self.a_and_b(x);

        (Array1::from_shape_fn
        (a.len(), |i| eta3.powi(i as i32))*a).sum()
    }
    pub fn int_2(&self,rho:f64,t:f64,x:&Array1<f64>)->f64{
        let eta3=self.csi(rho, t, x)[3];
        let (a,b)=self.a_and_b(x);

        let i2=
        (Array1::from_shape_fn
        (b.len(), |i| eta3.powi(i as i32))*b).sum();

        i2
    }
    pub fn deta_i1_and_deta_i2_deta(&self,rho:f64,t:f64,x:&Array1<f64>)->(f64,f64){
        let eta3=self.csi(rho, t, x)[3];
        let (a,b)=self.a_and_b(x);
        let deta_int1_deta=
        (Array1::from_shape_fn
        (a.len(), |i| (i as f64 + 1.0)*eta3.powi(i as i32))*a).sum();
       
        let deta_int2_deta=
        (Array1::from_shape_fn
        (b.len(), |i|  (i as f64 + 1.0)*eta3.powi(i as i32))*b).sum();

        (deta_int1_deta,deta_int2_deta)
    }
    
    pub fn c1(&self,rho:f64,t:f64,x:&Array1<f64>)->f64{
        let mmed=self.mmed(x);
        let eta=self.csi(rho, t, x)[3];

        (1. + mmed*(8.*eta-2.*eta.powi(2))/(1.-eta).powi(4) + (1.-mmed)*(20.*eta -
        27.*eta.powi(2) + 12.*eta.powi(3) - 2.*eta.powi(4))/((1.-eta)*(2.-eta)).powi(2)).powi(-1)

    }
    pub fn dc1dx(&self,rho:f64,t:f64,x:&Array1<f64>)->Array1<f64>{
        let c1= self.c1(rho, t, x);
        let c2=self.c2(rho, t, x);
        let dcsidx=self.dcsidx(rho, t, x);
        let eta=self.csi(rho, t, x)[3];
        
        let p1=&dcsidx.row(3)*c2;
        let p2=&self.m*(8.*eta-2.*eta.powi(2))/((1.-eta).powi(4));
        let p3=-&self.m*
        (20.-eta-27.*eta.powi(2)+12.*eta.powi(3)-2.*eta.powi(4))/((1.-eta)*(2.-eta)).powi(2);

        p1-(p2+p3)*(c2.powi(2))
    }
    pub fn c2(&self,rho:f64,t:f64,x:&Array1<f64>)->f64{
        let c1=self.c1(rho, t, x);
        let mmed=self.mmed(x);
        let eta=self.csi(rho, t, x)[3]; 

        -c1.powi(2)*
        (
        mmed*(-4.*eta.powi(2)+20.*eta+8.)/(1.-eta).powi(5)+
        (1.-mmed)*(2.*eta.powi(3)+12.*eta.powi(2)-48.*eta+40.)/((1.-eta)*(2.-eta)).powi(3)
        )
         
    }

    // # EQ A.12 and A.13 ok!
    pub fn mes_and_me2s(&self,t:f64,x:&Array1<f64>)->(f64,f64){

        let n = x.len();
        let m=&self.m;
        let sij=&self.sij;
        let eij=&self.eij;
        // let s=&self.s;

        // dbg!(sij);
        // dbg!(eij);
        let mut si_1=0.0;
        let mut si_2=0.0;
        for i in 0..n{
            let mut sj_1=0.0;
            let mut sj_2=0.0;
            for j in 0..n{

                let eijt=eij[(i,j)]/t;
                let xxmmsij3=x[i]*x[j]*m[i]*m[j]*sij[(i,j)].powi(3);
                let prod1=xxmmsij3*eijt;
                sj_1+=prod1;
                sj_2+=prod1*eijt;
            }
            si_1+=sj_1;
            si_2+=sj_2;
        }

        (si_1,si_2)
    }
    // eq. A.21;
    pub fn csi(&self,rho:f64,t:f64,x:&Array1<f64>)->Array1<f64>{

        let mdens=  molecules_total_density(rho);   
        let vdt= self.td_hard_segment_diameter(t);
        let powd=Array2::from_shape_fn((x.len(),4), |(i,j)| vdt[i].powi(j as i32));
        let csi=(x*&self.m).dot(&powd)*PI*mdens/6.0;
        csi        
        
    }


    pub fn td_hard_segment_diameter(&self,t:f64)->Array1<f64>{
        &self.s*(1.0-0.12* ((-3.*&self.e)/t).exp())
    }

    // A.10
    pub fn helmholtz(&self,t:f64,rho:f64,x:&Array1<f64>)->f64{
        let mmed=self.mmed(x);
        let i2= self.int_2(rho, t, x);
        let i1= self.int_1(rho, t, x);
        let c1=self.c1(rho, t, x);
        let (mes,me2s)=self.mes_and_me2s(t, x);
        let mdens=molecules_total_density(rho);

        -2.*PI*mdens*i1*mes - PI*mdens*mmed*c1*i2*me2s


    }

    pub fn dcsidx(&self,rho:f64,t:f64,x:&Array1<f64>)->Array2<f64>{
        let mdens=  molecules_total_density(rho);   
        let vdt= self.td_hard_segment_diameter(t);
        let powd=Array2::from_shape_fn((x.len(),4), |(i,j)| vdt[i].powi(j as i32));

        let m=&self.m;
        
        let mut mat=Array2::zeros((4,x.len()));

        for i in 0..4{
            for k in 0..m.len(){
                mat[(i,k)]=m[k]*powd[(k,i)];
            }
        }
        // 4Xn
        mat*mdens*PI/6.        
    }
    pub fn dadx_and_dbdx(&self,x:&Array1<f64>)->(Array2<f64>,Array2<f64>){
        let mmed=self.mmed(x);
        let ap1=&ArrayView1::from_shape(7, &AP[1]).unwrap();
        let ap2=&ArrayView1::from_shape(7, &AP[2]).unwrap();

        let bp1=&ArrayView1::from_shape(7, &BP[1]).unwrap();
        let bp2=&ArrayView1::from_shape(7, &BP[2]).unwrap();
        let m=&self.m;
        let mut ma=Array2::zeros((7,m.len()));
        let mut mb=Array2::zeros((7,m.len()));

        for i in 0..7{
            for k in 0..m.len(){

                ma[(i,k)]=ap1[i]*m[k]/(mmed.powi(2)) + ap2[i]*m[k]*(3.-4./mmed)/(mmed.powi(2));
                mb[(i,k)]=bp1[i]*m[k]/(mmed.powi(2)) + bp2[i]*m[k]*(3.-4./mmed)/(mmed.powi(2));

            }
        } 
        (ma,mb)

    }
    pub fn dint1dx_and_dint2dx(&self,rho:f64,t:f64,x:&Array1<f64>)->(Array1<f64>,Array1<f64>){

        let (a,b)=self.a_and_b(x);
        let eta3= self.csi(rho, t, x)[3];
        let dcsidx=self.dcsidx(rho, t, x); //ok


        let dcsidx_3=dcsidx.row(3);
        let n = 7;
        let ncomp=x.len();
        let mut d1dx=Array1::zeros(ncomp);
        let mut d2dx=Array1::zeros(ncomp);

        let (dadx,dbdx)=self.dadx_and_dbdx(x);
        //ok    

        // dbg!(dcsidx_3);
        // println!("dbg={dadx}");
        // dbg!(&dadx);
        // println!("bx={dbdx}");
        for k in 0..ncomp{
            let mut sum_1=0.0;
            let mut sum_2=0.0;
            
            for i in 0..n{
                
                sum_1+= a[i]*(i as f64)*dcsidx_3[i]*eta3.powi( i as i32-1) + dadx[(i,k)]*eta3.powi(i as i32);
                sum_2+= b[i]*(i as f64)*dcsidx_3[i]*eta3.powi( i as i32-1) + dbdx[(i,k)]*eta3.powi(i as i32)

            }
            d1dx[k]=sum_1;
            d2dx[k]=sum_2

        }
        dbg!(&d1dx);
        (d1dx,d2dx)

    }
    pub fn dmesdx_and_dme2sdx(&self,t:f64,x:&Array1<f64>)->(Array1<f64>,Array1<f64>){
        let n = x.len();
        let m=&self.m;
        let sij=&self.sij;
        let eij=&self.eij;
        // let s=&self.s;
        let mut dmesdx=Array1::zeros(n);
        let mut dme2sdx=Array1::zeros(n);

        for k in 0..n{
            let mut sj_1=0.0;
            let mut sj_2=0.0;
            for j in 0..n{

                let eijt=eij[(k,j)]/t;
                let xjmjsij3=x[j]*m[j]*sij[(k,j)].powi(3);
                let prod1=xjmjsij3*eijt;
                sj_1+=prod1;
                sj_2+=prod1*eijt;
            }
            dmesdx[k]=2.*m[k]*sj_1;
            dme2sdx[k]=2.*m[k]*sj_2;
        }

        (dmesdx,dme2sdx)
    }

    pub fn z(&self,rho:f64,t:f64,x:&Array1<f64>)->f64{

        let mmed=self.mmed(x);//ok
        let i2= self.int_2(rho, t, x);//ok
        let c1=self.c1(rho, t, x); //ok
        let c2=self.c2(rho, t, x); //ok
        let eta3=self.csi(rho, t, x)[3]; //ok
        let (mes,me2s)=self.mes_and_me2s(t, x); //ok
        let (di1,di2)=self.deta_i1_and_deta_i2_deta(rho, t, x); //ok

        // dbg!(di1,di2); ok

        // dbg!(mes,me2s); ok

        // dbg!(c1);
        // dbg!(c2);
        // dbg!(i2);

        let mdens=molecules_total_density(rho);

        -2.*PI*mdens*di1*mes - PI*mdens*mmed*(c1*di2 + c2*eta3*i2)*me2s
    }

    pub fn dadx(&self,t:f64,rho:f64,x:&Array1<f64>)->Array1<f64>{
        
        let mmed=self.mmed(x);
        let i2= self.int_2(rho, t, x);
        let i1= self.int_1(rho, t, x);
        let c1=self.c1(rho, t, x);
        let mdens=molecules_total_density(rho);

        let (di1dx,di2dx)= self.dint1dx_and_dint2dx(rho, t, x);
        println!("d1dx={di1dx}");
        println!("d2dx={di2dx}");
        let (mes,me2s)=self.mes_and_me2s(t, x);
        let (dmes,dme2s)=self.dmesdx_and_dme2sdx(t, x);
        let dc1dx=self.dc1dx(rho, t, x);

        // dbg!(&di1dx,&di2dx);

        let p1=-2.*PI*mdens*(di1dx*mes+i1*dmes);
        let p2=(&self.m*c1*i2+mmed*dc1dx*i2+mmed*c1*di2dx)*me2s;
        let p3=mmed*c1*i2*dme2s;

        p1-PI*mdens*(p2+p3)


    }


}
// impl EquationOfState for PCSaft {
//     fn ncomp(&self)->usize {
//         todo!()
//     }
//     fn bmix(&self,vx:&Array1<f64>)->f64 {
//         todo!()
//     }
//     fn density_solver(&self,t:f64,p:f64,vx:&Array1<f64>,guess:f64)->Result<f64,EosError> {
//         todo!()
//     }
//     fn dpdrho(&self,t:f64,rho:f64,vx:&Array1<f64>)->Result<(f64,f64,),EosError> {
//         todo!()
//     }

//     fn pressure(&self,t:f64,rho:f64,vx:&Array1<f64>)->Result<f64,EosError> {
//         todo!()
//     }
//     fn lnphi(&self,t:f64,rho:f64,vx:&Array1<f64>)->Result<Array1<f64>,EosError> {
//         todo!()
//     }
// }


mod tests{
    use std::sync::Arc;

    use approx::relative_eq;
    use ndarray::{Array1, Array2};


    use super::{ASCParameters, Dispersion, HardShere};


    fn get_disp()->(Dispersion,usize){
        let m=Array1::from_vec(
        vec![ 1.6686,1.2053, 2.0729,1.0   , 1.6069,2.0002,2.3316,2.2616,2.6896,2.5620,3.0576,3.4831, 3.8176 ,4.2079, 4.6627, 1.0  ]
        );
        let s=Array1::from_vec(
        vec![3.0349 ,3.313 , 2.7852,3.7039, 3.5206,3.6184,3.7086,3.7574,3.7729,3.8296,3.7983,3.8049, 3.8373 ,3.8448, 3.8384, 3.04 ]
        );
        let e=Array1::from_vec(
        vec![229.   ,90.96 , 169.21,150.03, 191.42,208.11,222.88,216.53,231.2 ,230.75,236.77,238.4 , 242.78 ,244.51, 243.87, 204.7]
        );
        let n=m.len();
        let mut kbi=Array2::zeros((n,n));

        kbi[(2,3)]=0.046;
        kbi[(3,2)]=0.046;

        kbi[(0,3)]=0.055;
        kbi[(3,0)]=0.055;

        kbi[(2,0)]=0.062;
        kbi[(0,2)]=0.062;

        ;
        let disp=Dispersion::new(m,s, e, kbi);
        (disp,n)
    }

    fn req(x:f64,a:f64)->bool{

        println!("Calc={x}");
        println!("ref={a}");
        relative_eq!
        (
        x, a, epsilon = 1e-12, max_relative = 1e-8)

    }
    #[test]
    fn test_1(){

        // ap = np.zeros([7, 3])
        // bp = np.zeros([7, 3])

        let zer1=Array1::zeros(2);
        let zer2=Array2::zeros((2,2));

        let disp=Dispersion::new(zer1.clone(),zer1.clone() , zer1.clone(), zer2.clone());
        let n=disp.a_and_b(&zer1.clone());
    }

    #[test]
    fn test_2(){

        let (d,n)=get_disp();

        let x=&Array1::from_elem(n, 1./n as f64  );

        let csi=d.csi(1e3, 300.,x );
        dbg!(csi);
    }   

    #[test]
    fn z_disp(){

        let (d,n)=get_disp();
        let x=&Array1::from_elem(n, 1./n as f64  );
        let rho=1e3;
        let t=300.;

        let z_disp=d.z(rho, t, x);
        let a=-0.8516168440922453;
        assert!(req(z_disp, a));


    }   

    #[test]
    fn a_disp(){
        let (d,n)=get_disp();
        let x=&Array1::from_elem(n, 1./n as f64  );
        let rho=1e3;
        let t=300.;

        let hel=d.helmholtz(t, rho, x);
        assert!(req(hel, -0.8859051456075235))
    }
    #[test]
    fn dadx_disp(){
        let (d,n)=get_disp();
        let x=&Array1::from_elem(n, 1./n as f64  );
        let rho=1e3;
        let t=300.;

        let dadx=d.dadx(t, rho, x);
        // dbg!(&dadx);

        assert!(req(dadx.product(), 175.66517851791482))
    }



}