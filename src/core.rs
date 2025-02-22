// Module that defines math.functions for EOSs

pub mod core_equations{

    use crate::parameters::{ asc::ASCParameters, srk::SRKParameters};
    use ndarray::{s, Array1, Array2, Array4, Axis, Zip};
    
    
    const EPS: f64 = 0.0;
    const SIG: f64 = 1.0;

    #[allow(unused)]
    pub const NEG: usize = 0;
    #[allow(unused)]
    pub const POS: usize = 1;

    pub const IDEAL_GAS_CONST: f64 = 8.31446261815324;
    
    pub struct SRKEquations{}

    pub struct ASCEquations{}
    

    impl SRKEquations {

        pub fn calc_bmix(vx:&Array1<f64>,vb:&Array1<f64>) -> f64 {vx.dot(vb)}
    
        pub fn calc_a_soave(t:f64,parameters:&SRKParameters)->Array1<f64>{
    
            let p = parameters;
            let va0 = &p.va0;
            let vtc = &p.vtc;
            let vkappa = &p.vkappa;
            let vtr =  t/vtc;
    
            let alpha_srk =  1.0+vkappa*(1.0 - vtr.sqrt());
            let alpha_srk = &alpha_srk*&alpha_srk;
    
            va0*alpha_srk
    
        }
        
    
        pub fn calc_sqrt_aij_matrix(t:f64,parameters:&SRKParameters)->Array2<f64>{
    
            let va = SRKEquations::calc_a_soave(t,parameters);
            let kij = &parameters.kij;
            let n = parameters.ncomp;
            let mut aij = Array2::<f64>::zeros((n,n));
    
            for i in 0..n {
                for j in i..n{
                    let prod = va[i]*va[j];
                    let sqrt_aij = prod.sqrt();
                    aij[(i,j)] = sqrt_aij*(1.0-kij[(i,j)]); 
                    aij[(j,i)] = aij[(i,j)]; 
                }
            }
            aij
        }
        pub fn calc_amix(t:f64,vx:&Array1<f64>,parameters:&SRKParameters)->f64{

            // a = x*Aij*x
            let aij_mat = SRKEquations::calc_sqrt_aij_matrix(t,parameters);
            let aij_mul_x: ndarray::ArrayBase<ndarray::OwnedRepr<f64>, ndarray::Dim<[usize; 1]>> = aij_mat.dot(vx);
            let a = vx.dot(&aij_mul_x);
            a
        }
    
        pub fn calc_dadn(t:f64,vx:&Array1<f64>,parameters:&SRKParameters)->Array1<f64>{
    
            let amix = Self::calc_amix(t,vx,parameters);
            let maij= Self::calc_sqrt_aij_matrix(t,parameters);
            let aij_dot_vx:Array1<f64> = maij.dot(vx);

            2.0*aij_dot_vx - amix
        }
    
        pub fn calc_press(t:f64,rho:f64,vx:&Array1<f64>,parameters:&SRKParameters)->f64{
    
            let am = Self::calc_amix(t,vx,parameters);
            let bm = Self::calc_bmix(vx,&parameters.vb);
            let r = IDEAL_GAS_CONST;

            (r * t) / (1.0/rho - bm) - am / ((1.0/rho + SIG * bm) * (1.0/rho + EPS * bm))
        }
    
        pub fn calc_phi(t:f64,rho:f64,vx:&Array1<f64>,parameters:&SRKParameters)->Array1<f64>{

            let vmolar: f64 = 1.0/rho;
            let t: f64 = t;
            let vb: &Array1<f64> = &parameters.vb;
            let vdbdn: &Array1<f64> = vb;
            let bmix: f64 = Self::calc_bmix(vx,&parameters.vb);
            let amix: f64 = Self::calc_amix(t,vx,parameters);
            let r: f64 = IDEAL_GAS_CONST;
            let vdadn: Array1<f64> = Self::calc_dadn(t,vx,parameters);
            let q: f64 = amix/(bmix*r*t);
            let pressure: f64 = Self::calc_press(t,rho,vx,parameters);
    
            let vdqdn: Array1<f64> = q*(1.0 + vdadn/amix - vdbdn/bmix);
            
            let i: f64 = (1.0/(SIG-EPS))*f64::ln((vmolar + SIG*bmix )/(vmolar + EPS*bmix));
    
            let z: f64 = (pressure*vmolar)/(r*t);
            
            let beta: f64 = z*rho*bmix;
            
            let vlnphi: Array1<f64> = (vdbdn/bmix)*(z-1.0) - f64::ln(z-beta) - vdqdn*i;
    
            vlnphi.exp()
        }
    
}


    // esse tipo de enum será útil na criação de metodos numericos
    #[derive(Debug)]
    pub enum NotConverging{
        NaNValues,
        MaxIterations(i32)
    }

    impl std::fmt::Display for NotConverging {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            match self {
                Self::NaNValues=> write!(f,"Method converged to NaN values"),
                Self::MaxIterations(it)=> write!(f,"Method didn't converged in {it} iterations."),
            }
        }
    }
    pub struct Converging {
            pub result:Array2<f64>,
            pub it:i32,
            pub tol: f64
    }

    impl std::fmt::Display for Converging {
        fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f,"Result = {} , it = {} , tol = {} " ,self.result,self.it,self.tol)
        }
    }
    impl ASCEquations {

        const R:f64 = IDEAL_GAS_CONST;
        
        pub fn calc_bmix(vx:&Array1<f64>,vb:&Array1<f64>) -> f64 {
            vx.dot(vb)
        }

        // g of mixture
        pub fn calc_radial_dist_mix(rho:f64,vx: &Array1<f64>,vb:&Array1<f64>)->f64{

            let bm = Self::calc_bmix(vx, vb);
            1.0 / (1.0 - 1.9 * (rho * bm / 4.0))
        }

        pub fn calc_dlngdrho(rho:f64,vx: &Array1<f64>,parameters:&ASCParameters)->f64{

            let bm = Self::calc_bmix(vx, &parameters.vb);
            let gm = Self::calc_radial_dist_mix(rho,vx,&parameters.vb);

            (1.0 / gm) * (-1.0) * (gm*gm) * (-1.9 * bm / 4.0)

        }

        pub fn calc_dlngdni(rho:f64,vx: &Array1<f64>,vb:&Array1<f64>)->Array1<f64>{

            let gmix = Self::calc_radial_dist_mix(rho, vx,vb);
            gmix * 1.9 *(rho*vb/4.0)
        }

        pub fn calc_vdelta(t:f64,rho:f64,vx:&Array1<f64>,parameters:&ASCParameters)->Array1<f64>
        {   
            let veps = &parameters.veps;
            let vbeta = &parameters.vbeta;
            let vb = &parameters.vb;
            let gmix = Self::calc_radial_dist_mix(rho, vx,&parameters.vb);

            gmix * ((veps / (Self::R * t)).exp() - 1.0) * vb * vbeta
        }

        pub fn calc_delta_cross(t:f64,rho:f64,vx: &Array1<f64>,vb:&Array1<f64>,eps:f64,beta:f64)->f64{

            let gmix = Self::calc_radial_dist_mix(rho, vx,vb);
            let bm = Self::calc_bmix(vx, vb);
            gmix * ((eps / (Self::R * t)).exp() - 1.0) *  bm* beta
        }

        pub fn calc_delta_cross_ecr(vdelta:&Array1<f64>,i:usize,j:usize)->f64
        {
            (vdelta[i]*vdelta[j]).sqrt()
        }


        pub fn calc_delta_mat(t: f64,rho:f64,vx: &Array1<f64>,parameters:&ASCParameters)->Array4<f64>{

            let ncomp = parameters.ncomp;

            // por enquanto clona
            let idx = &parameters.associative_components.index;
            // let n_assoc_comps = parameters.associative_components.ncomp;
            let vdelta = Self::calc_vdelta(t,rho,vx,parameters);

            let solvation = &parameters.solvation;
            // S - Comp - S - Comp
            let mut mdelta = Array4::zeros((2,ncomp,2,ncomp));

            // 
            for &comp1 in idx{
            // for comp1 in 0..ncomp{
                // tentar dps in comp1+1..ncomp
                for &comp2 in idx{
                // for comp2 in comp1..ncomp{

                    let k = &(comp1,comp2);
                    let kinv = &(comp2,comp1);

                    if comp1==comp2{

                        mdelta[(NEG,comp1,POS,comp2)] = vdelta[comp1];
                        mdelta[(POS,comp1,NEG,comp2)] = vdelta[comp1];

                    }
                    else {

                        if solvation.contains_key(k) | solvation.contains_key(kinv){

                            // cuidado com k,(mudar dps)

                            let (eps,beta) = match solvation.get(k){
                                Some(val) => *val,
                                None => *solvation.get(kinv).unwrap()
                            };

                            let delta_cross: f64 = Self::calc_delta_cross(t, rho, vx,&parameters.vb, eps, beta);

                            mdelta[(NEG,comp1,POS,comp2)] = delta_cross; 
                            mdelta[(POS,comp1,NEG,comp2)] = delta_cross;
                            
                            mdelta[(NEG,comp2,POS,comp1)] = delta_cross; 
                            mdelta[(POS,comp2,NEG,comp1)] = delta_cross;
                        }

                        else {

                            let delta_cross: f64 = Self::calc_delta_cross_ecr(&vdelta, comp1, comp2);
                            mdelta[(NEG,comp1,POS,comp2)] = delta_cross; 
                            mdelta[(POS,comp1,NEG,comp2)] = delta_cross;
                            
                            mdelta[(NEG,comp2,POS,comp1)] = delta_cross; 
                            mdelta[(POS,comp2,NEG,comp1)] = delta_cross;
                        }
                    }
                }
            }
            mdelta

        }


        pub fn calc_non_assoc_sites_mat(rho:f64,vx:&Array1<f64>,x_assoc:Option<Array2<f64>>,parameters:&ASCParameters,delta:&Array4<f64>)-> Result<Converging, NotConverging>{

            let ncomp = parameters.ncomp;
            let s_matrix = &parameters.s_matrix;
            const TOL:f64 = 1e-9;
            const MAX:i32 = 1000;
            let assoc_comps = parameters.associative_components.index.clone();
            let mut x_assoc = match x_assoc {
                
                Some(val)=> val,
                None=> { 
                    // Array2::from_elem((2,ncomp), 0.5)
                    let mut x_assoc:Array2<f64> = Array2::from_elem((2,ncomp), 1.0);
            
                    for &i in &assoc_comps{
        
                        x_assoc[(0,i)] = 0.5;
                        x_assoc[(1,i)] = 0.5;
        
                    }
                    x_assoc

                }
            };

            // let mut x_assoc: Array2<f64> = Array2::from_elem((2,ncomp), 0.5);
            let mut mat_error:Array2<f64> = Array2::zeros((2,ncomp));
            let mut x_old: Array2<f64>;


            let mut res = 1.0;
            let mut it:i32 = 0; 
            let mut error1=1.;
            let mut error0=1.;


            // println!("Dentro de calc_xassoc");
            // println!("{:#?}",&x_assoc);

            while (res>TOL) & (it<MAX) {

                // dbg!(res);
                it+=1;

                x_old= x_assoc.clone();

                
                for &i in &assoc_comps{

                    for j in  0..2{

                        let mut sum1: f64 = 0.0;

                        for &k in &assoc_comps{


                            let mut sum2 = 0.0;

                            for l in 0..2{
                                sum2+= s_matrix[(l,k)]*x_old[(l,k)]*delta[(l,k,j,i)]
                                // sum2+= s_matrix[(l,k)]*x_assoc[(j,i)]*delta[(l,k,j,i)]
                            }

                            sum1+= vx[k]*sum2
                        }

                        x_assoc[(j,i)] = 1.0/(1.0 + rho*sum1);
                        mat_error[(j,i)] = ((x_assoc[(j,i)]-x_old[(j,i)])/x_assoc[(j,i)]).abs()
                        // error1 = ((x_assoc[(j,i)]-x_old[(j,i)])/x_assoc[(j,i)]).abs();

                        // if error0<error1{
                        //     error0 = error1
                        // }

                    }

                }
                // reduce(f64::max(self, other))
                res = *mat_error.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
                // res = error1
            
            }

        
        // println!("it={}",it);

        if it == MAX{
            
            panic!("O cálculo da matriz de sítios não-associados não convergiu. Verifique as condições (T,P,x).")

        }

        // else if x_assoc.is_any_nan(){

        //     return Err(NotConverging::NaNValues);
        // }
        else {
            return Ok(Converging{ result: (x_assoc), it: (it), tol: (TOL)});
        }
        
        }


        pub fn assoc_molecules_frac(vx: &Array1<f64>,parameters:&ASCParameters,xassoc:&Array2<f64>)->f64{
            let s_matrix = &parameters.s_matrix;

            // let mut sum1 = 0.0;
            let ncomp = parameters.ncomp;
        
            let minus_one_mat = 1.0 - xassoc;

            // element wise product 
            let hprod = minus_one_mat*s_matrix;

            // 1xN
            let sum_2_mat = hprod.sum_axis(ndarray::Axis(0));

            // for i in 0..ncomp{

            //     // sum2 significa a quantidade total de sitios associados (tanto - quanto +)
            //     // para uma dada substância
            //     // dps, essa quantidade é ponderada pela fração desse componente na mistura
            //     // resultando, assim, no valor efetivo de sítios associados da subs. na mistura

            //     // de tal modo que sum1 será a quantidade efetiva de associação da mistura
            //     let mut sum2 = 0.0;

            //     for j in 0..2{

            //         sum2+= (1.0-xassoc[(j,i)])*s_matrix[(j,i)]

            //     }
            //     sum1+= vx[i]*sum2
            // }
            // sum1

            vx.dot(&sum_2_mat)

        }   
        pub fn calc_press(t: f64, rho:f64,vx: &Array1<f64>,parameters:&ASCParameters,xassoc:&Array2<f64>)->f64{

            let f_asc_molec = Self::assoc_molecules_frac( vx, parameters, xassoc,);
            let dlngdrho = Self::calc_dlngdrho(rho, vx, parameters);
            -IDEAL_GAS_CONST*t*((1.0/(2.*(1./rho)))*((f_asc_molec)*(1.+(1./(1./rho))*dlngdrho)))
        }
        pub fn calc_phi(rho:f64,vx: &Array1<f64>,xassoc:&Array2<f64>,parameters:&ASCParameters,)->Array1<f64>{

            let ncomp = parameters.ncomp;
            let s_matrix = &parameters.s_matrix;

            let f_asc_molec = Self::assoc_molecules_frac( vx, parameters, xassoc);

            let dlngdni = Self::calc_dlngdni(rho, vx, &parameters.vb);

            let ln_xassoc = xassoc.ln();
            let lnx_mult_by_smat:Array2<f64> = ln_xassoc*s_matrix;

            let slogx = lnx_mult_by_smat.sum_axis(ndarray::Axis(0));

            let lnphi = slogx - (f_asc_molec/2.0)*dlngdni;
            // let mut lnphi = Array1::<f64>::zeros(ncomp);
            // for i in 0..ncomp{
            //     let mut slogx = 0.0;
            //     for j in 0..2{

            //         slogx += xassoc[(j,i)].ln()*s_matrix[(j,i)]  
            //     }

            //     lnphi[i] = slogx - (f_asc_molec/2.0)*dlngdni[i]

            //     // slogx = 
            // }
            lnphi.exp()

        }
}
}
