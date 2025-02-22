// Module that defines the interface for EOS's equations (which are impl. in core.rs)
pub mod eos{

    pub use ndarray::Array1;
    use ndarray::Array2;
    use crate::parameters::srk::SRKParameters;
    use crate::parameters::asc::{ASCParameters, SchemeType};
    use crate::core::core_equations::{ASCEquations, Converging,NotConverging, SRKEquations};
    use std::collections::HashMap;
    use std::fs::File;
    use std::io::Read;
    use serde_json::Value;

    use crate::core::core_equations::IDEAL_GAS_CONST as R;

    #[allow(dead_code)]
    
    // Cubic Term //
    pub struct SRKeos
    {
        pub parameters: SRKParameters,
    }

    #[allow(dead_code)]
    impl SRKeos {
        
        pub fn bmix(&self,vx:&Array1<f64>)->f64{SRKEquations::calc_bmix(vx,&self.parameters.vb)}

        pub fn asoave(&self,t:f64)->Array1<f64>{SRKEquations::calc_a_soave(t,&self.parameters)}

        pub fn amix(&self,t:f64,vx:&Array1<f64>)->f64{SRKEquations::calc_amix(t,vx,&self.parameters)}  

        pub fn phi(&self,t:f64,rho:f64,vx:&Array1<f64>)->Array1<f64> {SRKEquations::calc_phi(t,rho,vx,&self.parameters)}

        pub fn pressure(&self,t:f64,rho:f64,vx:&Array1<f64>)->f64 {SRKEquations::calc_press(t,rho,vx,&self.parameters)}

    }


    // Association Term //
    #[allow(dead_code)]

    pub struct ASCeos
    {
        pub parameters: ASCParameters,
    }
    impl ASCeos {

        #[warn(dead_code)]
        pub fn pressure(&self,t:f64,rho:f64,vx:&Array1<f64>)->f64{
            //     Passoc=-_R*T*((1/(2*(1/rho)))*((fAscMolec)*(1+(1/(1/rho))*dg))) 
            
            let delta = &ASCEquations::calc_delta_mat(t, rho, vx, &self.parameters);
            let result = &ASCEquations::calc_non_assoc_sites_mat(rho, vx, None,&self.parameters, delta)
            .unwrap();

            let xassoc = &result.result;
            let it = result.it;
            // println!("Iterações Xassoc = {}",it);
            ASCEquations::calc_press(t, rho, vx, &self.parameters, xassoc)
            
        }

        pub fn calc_pressure_and_xassoc(&self,t:f64,rho:f64,vx:&Array1<f64>,x_assoc:Option<Array2<f64>>)->(f64,Array2<f64>){
            //     Passoc=-_R*T*((1/(2*(1/rho)))*((fAscMolec)*(1+(1/(1/rho))*dg))) 
            
            let delta = &ASCEquations::calc_delta_mat(t, rho, vx, &self.parameters);
            let result = ASCEquations::calc_non_assoc_sites_mat(rho, vx, x_assoc,&self.parameters, delta)
            .unwrap();

            let x_assoc = result.result;
            let it = result.it;
            // println!("Iterações Xassoc = {}",it);
            (ASCEquations::calc_press(t, rho, vx, &self.parameters, &x_assoc),x_assoc)
            
        }
        #[warn(dead_code)]
        pub fn x_assoc(&self,t:f64,rho:f64,vx:&Array1<f64>)-> Result<Converging,NotConverging>{

            let delta = ASCEquations::calc_delta_mat(t,rho,vx,&self.parameters);

            // let default:Array2<f64> = Array2::zeros((2,self.parameters.ncomp));
            let res: Result<Converging,NotConverging> = ASCEquations::calc_non_assoc_sites_mat(rho,vx,None,&self.parameters,&delta);
            res
        }
        #[warn(dead_code)]
        pub fn assoc_molecules_frac(&self,t:f64,rho:f64,vx:&Array1<f64>)->f64{

            let delta = ASCEquations::calc_delta_mat(t,rho,vx,&self.parameters);
            let res = ASCEquations::calc_non_assoc_sites_mat(rho, vx, None,&self.parameters, &delta)
            .unwrap()
            ;
            let xassoc = res.result;
            // let it = res.it;
            // dbg!(&xassoc,it);
            ASCEquations::assoc_molecules_frac(vx,&self.parameters,&xassoc)
        }       
        #[warn(dead_code)]
        pub fn phi(&self,t:f64,rho:f64,vx:&Array1<f64>)->Array1<f64>{
            
            let delta = ASCEquations::calc_delta_mat(t, rho, vx, &self.parameters);
            let xassoc = ASCEquations::calc_non_assoc_sites_mat(rho, vx, None,&self.parameters, &delta)
            .unwrap()
            .result;

            ASCEquations::calc_phi(rho, vx, &xassoc, &self.parameters)
        }
        #[warn(dead_code)]
        pub fn calc_pres_and_phi(&self,t:f64,rho:f64,vx:&Array1<f64>)->(f64,Array1<f64>){

            let delta = ASCEquations::calc_delta_mat(t, rho, vx, &self.parameters);
            let xassoc = ASCEquations::calc_non_assoc_sites_mat(rho, vx, None,&self.parameters, &delta)
            .unwrap()
            .result;

            (ASCEquations::calc_press(t, rho, vx, &self.parameters, &xassoc),ASCEquations::calc_phi(rho, vx, &xassoc, &self.parameters))
            
        }
    }


    // CPA //
    #[allow(dead_code)]

    pub struct CPAeos
    {
        pub srk: SRKeos,
        pub asc: ASCeos
    }

    // #[allow(dead_code)]
    impl CPAeos {
        // 
        pub fn new(comps: Vec<&str>)->Self{

            let mut path = std::env::current_dir().unwrap();

            path.push("p.json");

            let mut file: File = File::open(path).unwrap();

            let mut content = String::new();

            file.read_to_string(&mut content).unwrap(); 

            let s = content.as_str();

            let data: Value = serde_json::from_str(s).unwrap(); 

            let all_comps= data["Components"].as_object().unwrap();
            let cross_parameters= data["CrossParameters"].as_object().unwrap();

            let ncomp = comps.len();

            let mut va0 = Array1::<f64>::zeros(ncomp);
            let mut vb = Array1::<f64>::zeros(ncomp);
            let mut vtc = Array1::<f64>::zeros(ncomp);
            let mut vkappa = Array1::<f64>::zeros(ncomp);
            let mut veps = Array1::<f64>::zeros(ncomp);
            let mut vbeta = Array1::<f64>::zeros(ncomp);

            let mut map_solv:HashMap<(usize,usize), (f64,f64)> = HashMap::new();
            let mut map_kij:HashMap<(usize,usize), f64> = HashMap::new();

            let mut scheme:Vec<SchemeType> = Vec::with_capacity(ncomp);

            for i in 0..ncomp{

                let string_name_uppercase = comps[i].to_uppercase();

                let i_name = string_name_uppercase.as_str();

                let component = match all_comps.get(i_name){
                    Some(val)=> val.as_object().unwrap(),
                    None=> panic!("{i_name} is not presented in parameter's file."),
                };

                let comp_p = component.get("Parameters")
                .unwrap().as_object().unwrap();

                va0[i] = comp_p.get("va0").unwrap().as_f64().unwrap();
                vb[i] = comp_p.get("vb").unwrap().as_f64().unwrap();
                vtc[i] = comp_p.get("vtc").unwrap().as_f64().unwrap();
                vkappa[i] = comp_p.get("vc1").unwrap().as_f64().unwrap();
                veps[i] = comp_p.get("epsilon").unwrap().as_f64().unwrap();
                vbeta[i] = comp_p.get("beta").unwrap().as_f64().unwrap();

                let scheme_str= comp_p.get("scheme").unwrap().as_str().unwrap();

                scheme.push(SchemeType::from_str(scheme_str));

                if ncomp>1{

                    for j in (i+1)..ncomp{

                        let string_name_uppercase = comps[j].to_uppercase();

                        let j_name = string_name_uppercase.as_str();

                        let i_and_j = i_name.to_owned() + "_" + j_name;
                        let j_and_i = j_name.to_owned() + "_" + i_name;

                        if cross_parameters.contains_key(&i_and_j)|cross_parameters.contains_key(&j_and_i){

                            let cross_values = match cross_parameters.get(&i_and_j){
                                Some(val)=> val.as_object().unwrap(),
                                None => cross_parameters.get(&j_and_i).unwrap().as_object().unwrap()
                            };

                            let do_solvation = cross_values["solvation"].as_bool().unwrap();
                            let has_kij = cross_values.contains_key("kij");

                            if do_solvation{

                                let eps:f64 = cross_values["epsilon_cross"].as_f64().unwrap();
                                let beta:f64 = cross_values["beta_cross"].as_f64().unwrap();

                                map_solv.insert((i,j),(eps,beta));

                            }

                            if has_kij{

                                let kij_val =  cross_values["kij"].as_f64().unwrap();

                                map_kij.insert((i,j), kij_val);
                            }
                        }

                        else {
                            continue;
                        }
                    }
                }

            }

            let p_srk = SRKParameters{
                ncomp,
                va0,
                vb:vb.clone(),
                vkappa,
                vtc,
                kij: SRKParameters::fill_kij(map_kij, ncomp)
            };
            let tup =ASCParameters::define_association_comps(&scheme);
            let p_asc = ASCParameters{

                ncomp,
                veps,
                vbeta,
                vb,
                associative_components: tup.0,
                non_assoc_comps: tup.1,
                s_matrix: ASCParameters::mk_s_matrix(ncomp, &scheme),
                solvation: map_solv
            };

            Self{
                asc: ASCeos { parameters: p_asc },
                srk: SRKeos { parameters: p_srk}
            }

            }
        
        // 

        // não é utilizado
        pub fn pressure(&self,t:f64,rho:f64,vx:&Array1<f64>)->f64 {self.srk.pressure(t, rho, vx) + self.asc.pressure(t,rho,vx)}

        // versão menos custosa computacionalmente 
        pub fn calc_p_and_new_xassoc_from_xassoc(&self,t:f64,rho:f64,vx:&Array1<f64>,x_assoc:Array2<f64>)->(f64,Array2<f64>){
            let (p_asc,x_assoc)= self.asc.calc_pressure_and_xassoc(t,rho,vx,Some(x_assoc));

            (self.srk.pressure(t, rho, vx) + p_asc, x_assoc)
        }

        pub fn calc_p_from_xassoc(&self,t:f64,rho:f64,vx:&Array1<f64>,x_assoc:&Array2<f64>)->f64{

            let p_asc= ASCEquations::calc_press(t, rho, vx, &self.asc.parameters, x_assoc);
            self.srk.pressure(t, rho, vx) + p_asc
        }


        // evitar calculo repetido de xassoc e delta (sao op. custosas)
        pub fn calc_pres_cpa_and_phi_asc(&self,t:f64,rho:f64,vx:&Array1<f64>)->(f64,Array1<f64>){

            let delta = ASCEquations::calc_delta_mat(t, rho, vx, &self.asc.parameters);
            let xassoc = ASCEquations::calc_non_assoc_sites_mat(rho, vx, None,&self.asc.parameters, &delta)
            .unwrap()
            .result;

            let phi_asc = ASCEquations::calc_phi(rho, vx, &xassoc, &self.asc.parameters);
            // let p_cpa = self.pressure(t,rho,vx);
            let p_asc = ASCEquations::calc_press(t,rho,vx,&self.asc.parameters,&xassoc);
            let p_srk = SRKEquations::calc_press(t,rho,vx,&self.srk.parameters);

            (p_asc+p_srk,phi_asc)
        }
        pub fn phi(&self,t:f64,rho:f64,vx:&Array1<f64>)->Array1<f64>{

            let (p_cpa,phi_asc) = self.calc_pres_cpa_and_phi_asc(t, rho, vx);
            
            let (p_srk,phi_srk) = (self.srk.pressure(t, rho, vx), self.srk.phi(t, rho, vx));

            // println!("phi asc = {} , phi srk = {}",&phi_asc,&phi_srk);
            phi_srk*phi_asc*(p_srk/p_cpa)

        

        }

        // Definir metodo pra calcular volume pra cpa
        // Nao sei se defino tudo aqui dentro da cpa
        // ja que muita coisa seria apenas função auxiliar 
        // pra que o solver funcionasse
        // pub fn volume(&self,t:f64,rho_old:Option<f64>,p:f64,vx:&Array1<f64>,phase:&str)->Result<VolumeResult,VolumeSolverError>{
        pub fn volume(&self,t:f64,rho_old:Option<f64>,p:f64,vx:&Array1<f64>,phase:&str)->f64{

            let bm = self.srk.bmix(vx);
            let rhomax = 1./bm;
            let eps = 1e-5;

            // função objetivo do método de busca de raiz
            let f = |s:f64,x_assoc:Array2<f64>| { 
                let rho = s*rhomax;
                let (p_cpa,new_x_assoc) = self.calc_p_and_new_xassoc_from_xassoc(t, rho, vx, x_assoc);

                ( (1.0-s)*(p_cpa - p) , new_x_assoc ) 
                // (1.0-s)*(self.pressure(t, rho, vx) - p) 

            };
            let dfds = |s:f64,x_assoc:&Array2<f64>|{

                let s_mais = s+eps;
                let s_menos = s-eps;

                let rho_mais = s_mais*rhomax;
                let rho_menos = s_menos*rhomax;

                // parto da suposição de o x_assoc(s) nao é muito difrente de x_assoc(s+eps)

                let p_cpa_mais = self.calc_p_from_xassoc(t, rho_mais, vx, &x_assoc);
                let p_cpa_menos = self.calc_p_from_xassoc(t, rho_menos, vx, &x_assoc);


                let fwd = (1.0 - (s+eps) )*(p_cpa_mais - p) ;
                let bwd = (1.0 - (s-eps) )*(p_cpa_menos - p) ;

                (fwd - bwd)/(2.*eps)

            };


            let phase = Self::which_phase_is(t, p, bm, phase);

            let mut s1 = match rho_old {

                //rho= s*rhomax
                Some(val)=> val/rhomax,
                None=> match phase{
                    Phase::Liquid{initial_guess} => initial_guess,
                    Phase::Vapor{initial_guess} => {
                        initial_guess
                    },
                }
            };


            // var do metodo iterativo
            let mut f0: f64 = 1.;

            let mut s_min =0.;
            let mut s_max =1.;

            let tol = 1e-6;

            // let it_max = 100;

            let mut it =0;
            let ncomp = self.asc.parameters.ncomp;

            let assoc_comps = self.asc.parameters.associative_components.index.clone();


            let mut x_assoc:Array2<f64> = Array2::from_elem((2,ncomp), 1.0);
            
            for i in assoc_comps{

                x_assoc[(0,i)] = 0.5;
                x_assoc[(1,i)] = 0.5;

            }
            let mut df_at_s0:f64;
            let mut res: f64 = 1.;

            let mut x_assoc:Array2<f64> = Array2::from_elem((2,ncomp), 0.5);

            while (res.abs()>tol) & (it<100) {


                // let xassoc =...

                it+=1;

                let s0: f64 = s1;

                // println!("itaração dentro do solver={}",it);
                // let mut x_assoc:Array2<f64> = Array2::from_elem((2,ncomp), 0.5);

                // dbg!(s0);
                // dbg!(f0.abs());
                // dbg!(s_max,s_min);
                // println!("{:#?}",&x_assoc);

                ( f0,x_assoc ) = f(s0,x_assoc);


                if f0>0.{
                    s_max = s0;
                }

                else if f0<0.{
                    s_min = s0;
                }

                df_at_s0 = dfds(s0,&x_assoc);
                s1 = s0 - f0/df_at_s0;
                

                // println!("Residuo f={} , Res. x= {}",f0.abs(),res);
                // usar esse criterio nao traz erros, e é mais (dependendo, economiza ate 1s)
                res = ((s1-s0)/s0).abs();
                
                // res = f0.abs();

                if (s1>=s_min) & (s1<=s_max) {continue;}

                else {
                    
                    s1 = (s_max+s_min)/2.;

                }


            }

            let vol = 1./(rhomax*s1);
            // println!("{vol}")
            vol
            // if it == it_max{

            //     Err(VolumeSolverError::MaxIterations { x: (vol), it: (it), tol: (tol) })

            // }

            // else if vol.is_nan() {
            //     Err(VolumeSolverError::NaNValue)
            // }

            // else {
            //     Ok(VolumeResult{x:vol,it})
            // }

        }
            
        pub fn which_phase_is(t:f64,p:f64,bm:f64,phase:&str)->Phase{

            if phase.to_uppercase()== "VAPOR".to_string(){

                let guess = bm/(bm + (R*t)/p);

                Phase::Vapor { initial_guess: (guess) } 


            }

            else if phase.to_uppercase()== "LIQUID".to_string(){
                
                let guess = 0.99;
                Phase::Liquid { initial_guess: (guess) } 
                
            }

            else {
                panic!("Phase must be 'Liquid' or 'Vapor'.")
            }



        }

        // fn aux_func_to_calc_volume(&self,s:f64,p:f64,rhomax:f64){

        //     let rho = s*rhomax;

        //     p_at_s =  self.pressure(t, rho, vx)


            
        // }

        // fn volume_solver(){}
        // possivelmente vou ter calc_pres_and_phi aq.
    }
    pub enum Phase{
        Vapor{initial_guess:f64,},
        Liquid{initial_guess:f64},
    }
    // #[derive(Debug)]
    // pub struct VolumeResult{ pub x:f64, pub it:i32,}
    // #[derive(Debug)]

    // pub enum VolumeSolverError{
    //     NaNValue,
    //     MaxIterations{x:f64,it:i32,tol:f64}

    // }

    // Carrega consigo o valor do chute inicial pro solver de volume da cpa 

}
