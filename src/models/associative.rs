use ndarray::{Array1, Array2};
use crate::models::SITES;
use crate::parameters::association::{ASCParameters, AssociationRule};
use crate::residual::Residual;
use crate::state::eos::{EosError, EosResult};
use super::{IDEAL_GAS_CONST,NS};




// impl Parameters for ASCParameters {

//     fn new(data:JsonStruct,comps:Vec<&str>)->Self {
        
//         let (components_map,binary_map)=
//         Self::get_components_and_binary_parameters(data, &comps);  
//         let ncomp = comps.len();
//         let mut vb = Array1::<f64>::zeros(ncomp);

//         #[allow(non_snake_case)]
//         let mut S=Array2::<f64>::zeros((NS,ncomp));

//         let mut mepsilon_cross=Array2::<f64>::zeros((ncomp,ncomp));
//         let mut mbeta_cross=Array2::<f64>::zeros((ncomp,ncomp));
//         let mut nself:Vec<usize>=Vec::new();
//         let mut nsolvation:Vec<usize>=Vec::new();


//         for (i,comp_i) in components_map.iter().enumerate(){

//             let option_assoc = comp_i.assoc.clone();
//             vb[i]=comp_i.cubic.b;

//             if let Some(assoc_term)=option_assoc{

//                 match assoc_term {

//                     AssociativeTerm::Inert=>{continue;}

//                     AssociativeTerm::Solvate{na,nb,nc}=>{
//                         S[(A,i)]=na;
//                         S[(B,i)]=nb;
//                         S[(C,i)]=nc;
//                         nsolvation.push(i);
//                     }
//                     AssociativeTerm::Associative { epsilon, beta, na,nb,nc   }=>{

//                         S[(A,i)]=na;
//                         S[(B,i)]=nb;
//                         S[(C,i)]=nc;

//                         nself.push(i);
                        
//                         mepsilon_cross[(i,i)]= epsilon;
//                         mbeta_cross[(i,i)]=beta;
//                     }
//                 }
//             }
//         }
//         let n_asc=nself.len();
//         let n_s=nsolvation.len();
//         let nassoc=([nsolvation.clone(),nself.clone()]).concat();
        
//         // dbg!(n_asc);
//         // dbg!(n_s);
//         // dbg!(&asc_comps);

        
//         // o numero de vezes que if let Some(assoc) funciona tem que ser
//         // exatamente c
//         // se it>c, nao foi especificado corretamente o tipo de interação 
//         // na parte binaria (tem menos informação sobre a interação binaria)

//         // se it<c, tem componente que nao é associativo sendo dito associativo

//         let mut pbinary:HashMap<(usize,usize),(f64,f64)>=HashMap::new();
//         let mut interaction: Vec<AssociativeInteraction> = Vec::new();
        
//         let c= ((n_asc*(n_asc-1)))/2 + n_s;
//         // dbg!(c);

//         let mut it =0;

//         let interactions=binary_map.keys();
        
//         for (i,j) in interactions{

//             let (i,j)=(*i,*j);

//             if let Some(assoc)=&binary_map.get(&(i,j)).unwrap().assoc{
//                 it+=1;
//                 match assoc{
//                     AssociativeBinary::Normal { rule, eps_rule, a, b }=>{
//                         if a.is_some()||b.is_some(){
//                             assert!(a.is_some(),"Interaction {i}{j}: you forget a parameter");
//                             assert!(b.is_some(),"Interaction {i}{j}: you forget b parameter"); 
//                             pbinary.insert((i,j), (a.unwrap(),b.unwrap())); 
//                         }
//                         let eps_rule=eps_rule.unwrap_or_default();
//                         let int_new = AssociativeInteraction::new(i, j, *rule,eps_rule);
//                         interaction.push(int_new);

//                         assert_ne!(mepsilon_cross[(j,j)],0.0,"Parametro associativo eps de {j} deu 0. Verifique os parâmetros.");
//                         assert_ne!(mepsilon_cross[(i,i)],0.0,"Parametro associativo eps de {i} deu 0. Verifique os parâmetros.");
//                         assert_ne!(mbeta_cross[(j,j)],0.0,"Parametro associativo beta de {j} deu 0. Verifique os parâmetros.");
//                         assert_ne!(mbeta_cross[(i,i)],0.0,"Parametro associativo beta de {i} deu 0. Verifique os parâmetros.");

//                         // let epscross=0.5*(mepsilon_cross[(i,i)]+mepsilon_cross[(j,j)]);

//                         let epscross = eps_rule.cross(mepsilon_cross[(i,i)], mepsilon_cross[(j,j)]);
//                         // println!("epscross={epscross}");
//                         mepsilon_cross[(i,j)]=epscross;
//                         mepsilon_cross[(j,i)]=epscross;

//                         let bcross=(mbeta_cross[(i,i)]*mbeta_cross[(j,j)]).sqrt();
//                         mbeta_cross[(i,j)]= bcross;
//                         mbeta_cross[(j,i)]= bcross;
//                     }

//                     AssociativeBinary::Solvation { rule,epsilon_cross, beta_cross }=>{
                    
//                     let int_new = AssociativeInteraction::new(i, j, *rule,EpsilonRule::Art);
//                     interaction.push(int_new);
                    
//                     let (assoc_idx, solv_idx, eps_assoc) = if nsolvation.contains(&j) {
//                         assert_ne!(mepsilon_cross[(i, i)], 0.0);
//                         assert_eq!(mepsilon_cross[(j, j)], 0.0);
//                         (i, j, mepsilon_cross[(i, i)])
//                     } else {
//                         assert_ne!(mepsilon_cross[(j, j)], 0.0);
//                         assert_eq!(mepsilon_cross[(i, i)], 0.0);
//                         (j, i, mepsilon_cross[(j, j)])
//                     };
//                     let eps_cross_val = match rule {
//                         AssociationRule::EXP => epsilon_cross.unwrap(),
//                         AssociationRule::MCR1 => eps_assoc * 0.5,
//                         _=>{panic!("rule of solvation interaction must be 'exp' or 'mcr1'!")}
//                     };
                    
//                     mepsilon_cross[(assoc_idx, solv_idx)] = eps_cross_val;
//                     mepsilon_cross[(solv_idx, assoc_idx)] = eps_cross_val;
//                     mbeta_cross[(assoc_idx, solv_idx)] = *beta_cross;
//                     mbeta_cross[(solv_idx, assoc_idx)] = *beta_cross;
//                     }
//                 }   
//             }
//         }


//         // balanço pra verificar se n to fazendo associação de forma equivocada já de inicio
//         if c>it{ panic!("Bsible Interactions = {c} > Valid Interactions verified (from json) = {it} (underspecified) ") }
//         if c<it{ panic!("Bsible Interactions = {c} < Valid Interactions verified (from json) = {it} (overspecified)") }
//         let kronecker=S.map(
//             |&x|
//             if x!=0.0{
//                 x/x
//             }else{
//                 0.0
//             }
            
//         );
//         let associative_parameters = ASCParameters{
//             ncomp,
//             vb,
//             nself,
//             map_interactions:HashMap::new(),
//             nassoc,
//             nsolv: Vec::new(),
//             site_multiplicity: S,
//             kronecker,
//             eps_cross_mat:mepsilon_cross,
//             beta_cross_mat:mbeta_cross,
//             interaction: interaction,
//             pbinary
//         };

//         associative_parameters
//     }
// }

pub struct Associative{
    pub parameters:ASCParameters
}
// pub struct Associative{
//     pub parameters:ASCParameters,
//     pub X:RefCell< X>,
// }

// struct X{
//     pub data:Array2<f64>,
//     pub rho:f64,

// }

// impl X {
    
//     fn default(rho:f64,ncomp:usize)->Self{
//         Self{
//             data:Array2::<f64>::from_elem((NS,ncomp),0.2),
//             rho
//         }
//     }
//     pub fn change_coordinate(&mut self,rho:f64,vb:&Array1<f64>,x:&Array1<f64>){
        
//         let rhomax=1.0/vb.dot(x);
//         let ncomp=x.len();
//         let s0=self.rho/rhomax;
//         let s1=rho/rhomax;
//         let eps=0.4;

//         // println!("old xassoc{}",&self.data);
//         let dif = (s1-s0).abs();
//         // dbg!(rho);
//         // dbg!(dif);
//         // println!("{}",x);
//         // dbg!(dif);
//         self.rho=rho;
//         // println!("h={dif}");
//         if dif<=eps{
//             //Nao modificar X
//             // apenas no final do Xassoc que será atribuido o novo valor 
//             // associado a rho novo

//         }else {
            
//             //retorna pro chute inicial
//             // println!("Mudança de coordenada");
//             // println!("Mudança de coordenada");
//             self.data=Array2::<f64>::from_elem((NS,ncomp),0.2)
//         }

//     }
// }

impl Associative {
    
    // pub fn new(p:ASCParameters)->Self{
    //     let ncomp=p.ncomp;
    //     Self{parameters:p,X:
    //         RefCell::new(X::default(0.0, ncomp))}
    // }
    pub fn new(p:ASCParameters)->Self{
        let ncomp=p.ncomp;
        Self{parameters:p}
    }

    pub fn set_binary(
        &mut self,
        i:usize,
        j:usize,
        rule:AssociationRule,
        eps:Option<f64>,
        beta:Option<f64>){self.parameters.set_binary_interaction(i, j, rule, eps, beta);}
            



}
impl Associative {
    pub fn g_func(&self,rho:f64,vx:&Array1<f64>)->f64{

        let bm = self.parameters.vb.dot(vx);
        1.0 / (1.0 - 1.9 * (rho * bm / 4.0))
    }

    pub fn dlngdrho(&self,rho:f64,vx:&Array1<f64>)->f64{

        let bm = self.parameters.vb.dot(vx);
        let gm =self.g_func(rho,vx);

        let result = (1.0 / gm) * (-1.0) * (gm*gm) * (-1.9 * bm / 4.0);
        // println!("dlngdrho={result}");
        result

    }

    // !!!
    pub fn dlngdni(&self,rho:f64,vx:&Array1<f64>)->Array1<f64>{

        let gmix = self.g_func(rho,vx);
        gmix * 1.9 *(rho*&self.parameters.vb/4.0)
    }

        // !!!
    pub fn h_func(&self,vx:&Array1<f64>,xassoc:&Array2<f64>)->f64{

        let s_matrix = &self.parameters.site_multiplicity;
       #[allow(non_snake_case)]
        let mut sxX=0.0;

        for i in 0..self.parameters.ncomp{
            #[allow(non_snake_case)]
            let mut sX=0.0;

            for j in 0..NS{

                sX+=(1.0-xassoc[(j,i)])*s_matrix[(j,i)].n() 

            }
            sxX+=vx[i]*sX;

        }

        sxX

    }   
}

impl Associative {
    #[allow(non_snake_case)]
    
    pub fn CR1(&self,t:f64,rho:f64,vx:&Array1<f64>,i:usize,k:usize)->f64{

        let p = &self.parameters;
        // let ij=&(i,j);
        // let ji=&(j,i);

        // let binary_parameters= p.pbinary.get(ij).or_else(||p.pbinary.get(ji));

        let eps_cross= p.eps_cross_mat[(i,k)];
        let beta_cross= p.beta_cross_mat[(i,k)];

        // let lij=binary_parameters.map_or(0.0, |(a,b)| a*t+b);

        // let eps_cross_modf=eps_cross*(1.0-lij);        
        // let eps_cross_modf=eps_cross*(1.0-lij);        
        let gmix = self.g_func(rho,vx);

        ((p.vb[i]+p.vb[k])*0.5)*beta_cross * gmix *(
            (( 
                ( eps_cross) / (IDEAL_GAS_CONST * t) ) ).exp() -1.0 )
            
    }

    #[allow(non_snake_case)]
    pub fn ECR(&self,t:f64,rho:f64,vx:&Array1<f64>,i:usize,k:usize)->f64
    {
        (self.CR1(t, rho, vx, i, i)*self.CR1(t, rho, vx, k,k)).sqrt()
    }

    pub fn delta(&self,t:f64,rho:f64,vx:&Array1<f64>)->Array2<f64>{
        
        let nassoc=&self.parameters.nassoc;
        let n=self.parameters.ncomp;
        let mut delta=Array2::zeros((n,n));
        for &i in nassoc{
        for &k in nassoc{
            if i>k{continue;}

            
            else{

                match self.parameters.binary[(i,k)].rule{
                AssociationRule::CR1=>{
                    delta[(i,k)]=self.CR1(t, rho, vx, i, k)
                }
                AssociationRule::ECR=>{
                    delta[(i,k)]=self.ECR(t, rho, vx, i, k)
                }
                //mCR1 e exp
                _=>{
                    delta[(i,k)]=self.CR1(t, rho, vx, i, k)
                    }
                }

                delta[(k,i)]=delta[(i,k)]
            }
        }
    }
        delta
    }

    #[allow(non_snake_case)]
    pub fn X_tan(&self,t:f64,rho:f64,vx:&Array1<f64>)-> Result<Array2<f64>,EosError>{

        let ncomp = self.parameters.ncomp;
        let S = &self.parameters.site_multiplicity;
        let delta=self.delta(t, rho, vx);
        // let kronecker=&self.parameters.kronecker;

        const TOL:f64 = 1e-11; // 11 estava OK
        // const MAX:i32 = 1000;
        const MAX:i32 = 10000;

        let assoc_comps = self.parameters.nassoc.clone();
        // let nassoc=assoc_comps.len();

        let mut x_assoc =  Array2::from_elem((NS,ncomp),0.2);

        let mut mat_error:Array2<f64> = Array2::zeros((NS,ncomp));
        let mut x_old: Array2<f64>;
        let mut res = 1.0;
        let mut it:i32 = 0; 

        while (res>TOL) & (it<MAX) {

            it+=1;
            x_old= x_assoc.clone();
            
            
            for &i in &assoc_comps{

                for s1 in SITES {

                    let mut sum1: f64 = 0.0;

                    for &k in &assoc_comps{

                        let mut sum2 = 0.0;

                        for s2 in SITES{

                            sum2+= x_old[(s2,k)]*S[(s2,k)].n()*(&S[(s1,i)]*&S[(s2,k)])
                        }

                        sum1+= delta[(i,k)]*vx[k]*sum2
                    }
                    // println!("order={}",1.0/(1.0 + rho*sum1));
                    x_assoc[(s1,i)] = 1.0/(1.0 + rho*sum1)  ;
                    mat_error[(s1,i)] = ( (x_assoc[(s1,i)]-x_old[(s1,i)])/x_assoc[(s1,i)] ).abs();
                }
            }
            // println!("{}  {}  {} {}",x_assoc[(0,0)], x_assoc[(1,0)], x_assoc[(2,1)], it);
            res = *mat_error.iter().max_by(|a,b| a.total_cmp(b)).unwrap();
            //res = mat_error.sum();
            // println!("{}",res);
        }

    // dbg!(&x_assoc);
    // dbg!(res);
    if it == MAX{
        return Err(EosError::NotConverged("X_tan".to_string()));
        
    }

    else {
        
        // x_grande.data=x_assoc.clone();

        // println!("X={}",&x_assoc);
        return Ok(x_assoc);
    }
    
}   


}


impl Residual for Associative {

    fn components(&self)->usize {
        self.parameters.ncomp
    }
    fn bmix(&self,x:&Array1<f64>)->f64 {
        self.parameters.vb.dot(x)
    }
    fn pressure(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<f64> {

        // let delta=&self.calc_delta_mat(t, rho, x);
        let xassoc=self.X_tan(t,rho, x)?;

        let h = self.h_func(x,&xassoc);
        let dlngdrho = self.dlngdrho(rho,x);

        Ok(
        -IDEAL_GAS_CONST*t*((1.0/(2.*(1./rho)))*((h)*(1.+(1./(1./rho))*dlngdrho)))

        )
    }
    fn residual_chemical_potential(&self,t:f64,rho:f64,x:&Array1<f64>)->EosResult<Array1<f64>>{

        // let delta=&self.calc_delta_mat(t, rho, x);
        let xassoc=self.X_tan(t,rho, x)?;
        let h = self.h_func(x,&xassoc);
        let nassoc=&self.parameters.nassoc;
        let mut mu = -0.5*h*self.dlngdni(rho,x);

        let ln_xassoc = xassoc.ln();
        let sites=&self.parameters.site_multiplicity;
        // let lnx_mult_by_smat:Array2<f64> = ln_xassoc*&self.parameters.site_multiplicity;
        // let lnx_mult_by_smat:Array2<f64> = ln_xassoc;
        // let mut mu=

        //varre tbm components inertes , pois o residual deles é não-nulo
        for &i in nassoc{
            let mut slogx=0.0;
            for j in SITES{
                slogx+=ln_xassoc[(j,i)]*sites[(j,i)].n()
            }
            mu[i]+=slogx 

        } 
        // let slogx = lnx_mult_by_smat.sum_axis(ndarray::Axis(0));
        // println!("mu={}",&mu);
        Ok(
        mu
        )

        // Ok(
        // slogx - (h/2.0)*dlngdni
        // )
    }
    
}

impl std::fmt::Display for ASCParameters {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, " === Associative Parameters ===")?;
        writeln!(f, "  Number of Components (NC): {}", self.ncomp)?;

        writeln!(f, "  Hard-sphere Volumes (vb):")?;
        writeln!(f, "    {:?}", self.vb.to_vec())?;

        writeln!(f, "  Epsilon Cross Matrix (εᵢⱼ):")?;
        for row in self.eps_cross_mat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }

        writeln!(f, "  Beta Cross Matrix (βᵢⱼ):")?;
        for row in self.beta_cross_mat.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }

        writeln!(f, "  Association Site Matrix (Sⱼᵢ):")?;
        for row in self.site_multiplicity.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }

        writeln!(f, "  Components (Self-Assoc. + Solvates): {:?}", self.nassoc)?;

        writeln!(f, "  Binary Parameters (rule and lᵢⱼ=aT+b):")?;
        for row in self.binary.rows() {
            writeln!(f, "    {:?}", row.to_vec())?;
        }
        // writeln!(f, "  εᵢⱼ Rule): {:?}", self.)?;


        // writeln!(f, "  Solvated Components: {:?}", self.non_assoc_comps)?;

        // writeln!(f, "  Association Interactions: {:?}", self.interaction)?;
        Ok(())
    }
}


mod tests{


    #[test]
    fn init(){

        // let comps = vec!["water","methanol"];

        // let ncomp = comps.len();

        // let set_1= "konteo_1_mod.json";

        // let p=50e5;
        // let t=300.0;
        // let vx=Array1::from_vec(vec![0.2,0.8]);
        // // let cubic = SRKParameters::new(set_1, comps.clone());
        // // let assoc= ASCParameters::new(set_1, comps.clone());
        // let cubic = SRKParameters::from_file(set_1, comps.clone());
        // let assoc= ASCParameters::from_file(set_1, comps.clone());


        // // dbg!(assoc);
        // // let y_dry_gas=Array1::from_vec(vec![1.0,1e-16]);

        // let cpa = CPAeos::new(cubic.clone(),assoc.clone());
        // // let rho = 1.0/0.02;

        // let s = State::
        // new_tpx(Arc::new(cpa), t, p, vx, DensityInitialization::Vapor).unwrap();

        // // let s= 
        // // State::new_trx(Arc::new(cpa), t, rho, vx).unwrap();

        // let phi=s.lnphi().unwrap().exp();
        // println!("volume={}",1./s.rho);
        // println!("phi={phi}");


    }

    // #[test]
    // fn test_dxdrho(){

    //     let comps = vec!["co2","water"];
    //     // let ncomp = comps.len();

    //     let set_1= "konteo_1.json";

    //     // let p=500e5;
    //     // let t=298.15;
    //     let vx=Array1::from_vec(vec![0.5,0.5]);
    //     let cubic = SRKParameters::new(set_1, comps.clone());
    //     let assoc= ASCParameters::new(set_1, comps.clone());

    //     // let y_dry_gas=Array1::from_vec(vec![1.0,1e-16]);
    //     let cpa = CPAeos::new(cubic.clone(),assoc.clone());
    //     // let assoc_eos= ASCeos::new(assoc.clone());
    //     let eos = Arc::new(cpa);

    //     let t=298.15;
    //     let vp = Array1::linspace(1.0, 700., 100);
    //     let mut vrho=Array1::zeros(100);

    //     let mut dFdr = Array1::zeros(100);
    //     let mut rho_dlng=Array1::zeros(100);

    //     for (i,p) in vp.iter().enumerate(){

    //         let p = 1e5*p;
    //         vrho[i]= State::new_tpx(eos.clone(), t, p, vx.clone(), crate::eos::density_solver::DensityInitialization::Vapor).unwrap().rho;

    //         let rho_mais =vrho[i] + 1e-5;
    //         let rho_menos = vrho[i] - 1e-5;

    //         let delta_mais = &eos.asc.calc_delta_mat(t, rho_mais, 
    //             &vx);
    //         let delta_menos = &eos.asc.calc_delta_mat(t, rho_menos, &vx);


    //         let xassoc_mais=eos.asc.calc_non_assoc_sites_mat(  rho_mais,&vx,None, delta_mais);
    //         let xassoc_menos=eos.asc.calc_non_assoc_sites_mat( rho_menos,&vx,None, delta_menos);

    //         let df = xassoc_mais.unwrap()-xassoc_menos.unwrap();
    //         let dfdrho = (df/(2.*1e-5)) * &eos.asc.parameters.s_matrix;
            
    //         rho_dlng[i]=vrho[i]*eos.asc.calc_dlngdrho(vrho[i], &vx);

    //         dFdr[i]=-dfdrho.sum()*vrho[i]*(1.+rho_dlng[i]);

    //     }

    //     println!("vr={}",vrho);
    //     println!("vp={}",vp);
    //     println!("dFdr={}",dFdr);
    //     println!("rdlngdrho={}",rho_dlng)

    // }



}
