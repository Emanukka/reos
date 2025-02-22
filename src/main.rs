// Imports
mod eos;
mod core;
mod parameters;
mod stability;
mod tools;
/// TODO:
/// 
/// [x]:Extrair parâmetros de um arquivo JSON (ou de um txt mesmo?)
///    "Operating on untyped JSON values , untyped objcts"
///    Pelo visto é o meu caso, onde não é possível usar structs
///    pra extrair os valores no arquivo
/// 
/// [x]: Equações de estado (SRK,ASC,CPA)
/// 
/// [x]: Cálculo de volume
/// [x]: TPD
/// [x]: Cálculo de saturação
/// 
/// 
/// [x]: Programa simples que pede do usuario o nome dos componentes, t,p, composição
/// []: Otimizar (com agua+co2 ta dando aproximadamente 4 segundos em altas pressoes);
///     ou seja, tenho que verificar as partes do código que dão mais gargalo:
///     
///     - a leitura dos parametros (os parametros podem estar declarados ja no binário
///       uma vez que o usuário não vai mexer nos parametros do nada);
///       (testei aqui, e nao acrescenta nenhum centésimo de segundo - talvez isso mude
///       pra +10 comps, mas provavelmente não vai ser significativo)
/// 
/// 
///     - Cálculo da matriz Xassoc , Cálculo de volume, TPD, watcon:
///
///        evitar clone() no calculo de xassoc e da TPD ; evitar loop no xassoc ... (mudou nada)
/// 
///     1. Consegui acelerar um pouco modificando o solver de volume, utilizando
///        xassoc da iteração anterior como chute; no entanto, o calculo de sat.
///        continua lento pros 15 componentes
/// 
///     2. fiz pequenas modificaçãos, mas continuou lento
///        
///     
/// 
///     

use crate::eos::eos::CPAeos;
use ndarray::Array1;
use serde::de::value::Error;
use crate::stability::watcon;
use std::{result, time::Instant};
// const  VEC:Vec<&str> = vec["H2S","N2", "CO2","CH4","C2","C3","C4",
// "iC4","C5","iC5","C6","C7","C8","C9" ,"C10","water"];

// #[derive(Debug)]
// struct  ComponentNotPresent{}
// fn read_user_input()->Result<Vec<String>,ComponentNotPresent>{
    
//     let before = Instant::now();
        
//     let mut user_input = String::new();
//     // let mut comp_names: Vec<&str> = Vec::new(); 
//     // user_input.inser

//     let r_result = std::io::stdin().read_line(&mut user_input).unwrap();
//     let gb = String::new();
//     let vec_result: Vec<&str> = user_input.trim()
//     .split(|c| c == ',' || c == ' ')
//     .collect();
//     let mut comps_names = Vec::<String>::new();
//     for i in vec_result{
//         if i !=""{
//             comps_names.push(i.to_string());
            
//         }
//     }


//     // verifca se os componentes estão
//     let comparison:Vec<&str> = vec!["h2s","n2", "co2","ch4","c2","c3","c4",
//     "ic4","c5","ic5","c6","c7","c8","c9" ,"c10","water"];
//     for comp in &comps_names{

//         if !comparison.contains(&comp.as_str()) {
            
//             return Err(ComponentNotPresent{});
            
//         }else {
//             continue;
//         }
//     }

//     Ok(comps_names)

// }

fn read_pressure()->Result<f64,bool>{

    let mut user_input = String::new();
    // let mut comp_names: Vec<&str> = Vec::new(); 
    // user_input.inser

    println!("Insira o valor da pressão (bar)");
    let r_result = std::io::stdin().read_line(&mut user_input).unwrap();
    let gb = String::new();
    
    let p:f64 = match user_input.trim().parse::<f64>(){

        Ok(val)=>val,
        Err(e)=>{
            println!("O valor digitado não é um número.");
            return Err(false)
        } 

    };

    if (p.is_sign_positive()) & (p!=0.0){
        Ok(p*1e5)
    }else {
        println!("Erro: Pressão não pode ser negativa ou igual à 0");
        Err(false)
    }


}

fn read_temperature()->Result<f64,bool>{
    let mut user_input = String::new();
    // let mut comp_names: Vec<&str> = Vec::new(); 
    // user_input.inser
    println!("Insira o valor da temperatura (K):");
    let r_result = std::io::stdin().read_line(&mut user_input).unwrap();

    let t:f64 = match user_input.trim().parse::<f64>(){

        Ok(val)=>val,
        Err(e)=>{
            println!("O valor digitado não é um número.");
            return Err(false)
        } 

    };

    if (t.is_sign_positive()) & (t!=0.0){
        Ok(t)
    }else {
        println!("Error: Temperatura não pode ser negativa ou igual à 0");
        Err(false)
    }

}


fn read_co2()->f64{
    let mut user_input = String::new();
    // let mut comp_names: Vec<&str> = Vec::new(); 
    // user_input.inser
    println!("Insira o valor da composição de CO2:");
    let r_result = std::io::stdin().read_line(&mut user_input).unwrap();
    let val:f64 = user_input.trim().parse::<f64>().unwrap();

    if (val.is_sign_positive()){
        val
    }else {
        println!("Composição não pode ser negativa");
        panic!()
    }

}

fn read_ch4()->f64{
    let mut user_input = String::new();
    // let mut comp_names: Vec<&str> = Vec::new(); 
    // user_input.inser
    println!("Insira o valor da composição de CH4:");
    let r_result = std::io::stdin().read_line(&mut user_input).unwrap();
    let val:f64 = user_input.trim().parse::<f64>().unwrap();

    if (val.is_sign_positive()){
        val
    }else {
        println!("Composição não pode ser negativa");
        panic!()
    }

}

fn read_user_input()->String{

    let mut user_input = String::new();
    std::io::stdin().read_line(&mut user_input).unwrap();
    user_input.trim().to_string()

}


fn main() {

    let raw = r#"
    [Opções]:

    - [1] ==> [Determina novo estado termodinâmico (T,P)]

    - [2] ==> [Calcula saturação baseada na composição em base seca de CO2 e CH4 (Alcanos,N2,H2S serão traços)]

    - [exit]
    "#;
    let comps_names:Vec<&str> = vec!["h2s","n2", "co2","ch4","c2","c3","c4",
    "ic4","c5","ic5","c6","c7","c8","c9" ,"c10","water"];

    let eos = CPAeos::new(comps_names);
    let ncomp = eos.asc.parameters.ncomp;
    let cond = true;

    let mut opp: Option<f64> = None ;
    let mut opt: Option<f64> = None;

    println!("{raw}");
    // let mut watcon_result:f64; 
    while cond{

        println!("[Insira uma opção]:");
        let s = read_user_input();
        let input = s.as_str();


        match input{

            // definir novo estado termodinamico
            "1"=>{

                let res_t = match read_temperature() {

                    Ok(value)=> value,
                    Err(e)=> continue,
                    
                };
                let res_p = match read_pressure() {

                    Ok(value)=> value,
                    Err(e)=> continue,
                    
                }; 
                (opp,opt) = (Some(res_p),Some(res_t));
            }

            // calcular saturação, especificando CO2 e CH4 ( o resto será traço)

            "2"=>{

                if opp.is_none(){ println!("Erro: (T,P) não definidos. Defina um estado termodinâmico."); continue; }

                let (p,t) = (opp.unwrap(),opt.unwrap());

                let (co2,ch4) = (read_co2(),read_ch4());


                let sum: f64 = co2+ch4;
                let mut trace = 1e-15;
                if sum>1.0{
                    println!("Erro: A soma das composições de CO2 e CH4 são maiores do que 1. Tente novamente.");

                    continue;
                }else if sum!=1.0 {
                    trace = (1.-sum)/((ncomp-3)as f64);
                    
                }
                

                let mut y_dry_gas:Array1<f64> = Array1::from_elem(ncomp,trace );
                y_dry_gas[ncomp-1] = 0.0; 
                y_dry_gas[2] = co2; 
                y_dry_gas[3] = ch4; 



                // println!("{}",y_dry_gas.sum());
                // println!("{}",&y_dry_gas);
                let before = Instant::now();
                let watcon_result = watcon(&eos, t , p , &y_dry_gas);
                let sat_in_ppm = watcon_result*1e6;

                let exec_time = before.elapsed();
                // println!("Tempo de execução: {:.2?}", before.elapsed());


                // println!("Composição traço = {} \n",trace);

                // let result = r#"
                // {}
                // "#;

                println!(
                    r#"

                    [Condição Termodinâmica] = (T= {}K, P= {}bar)
                    [Saturação de água] = {}ppm 

                    [Composição]:
                        [CO2] = {}
                        [CH4] = {}
                        [Alcanos,N2,H2S] = {}


                    [Tempo de Execução] = {:.2?}

                    "#
                    ,t,(p*1e-5).round(),sat_in_ppm,co2,ch4,trace,exec_time);



            }

            "exit"=>{
                break;
            }
            _=>{

                
            }
            


        }


    }


    // let p = read_pressure() *1e5;
    // let t= read_temperature();

    // let composition_of_all_componentes = 1./(ncomp-1)as f64;
    // let mut y_dry_gas:Array1<f64> = Array1::from_elem(ncomp,composition_of_all_componentes );
    // y_dry_gas[ncomp-1] = 0.0; 


    // let before = Instant::now();

    // let watcon_result = watcon(&eos, t, p, &y_dry_gas);

    // println!("Composição em base seca de todos os componentes = {} \n",composition_of_all_componentes);

    // println!("Saturação de água = {} ( T= {}K, P= {}bar )",watcon_result,t,(p*1e-5).round());
    // let before = Instant::now();

    // println!("Tempo de execução: {:.2?}", before.elapsed())



}


#[cfg(test)]
#[allow(dead_code)]
mod tests {
    use core::prelude::v1;
    use std::collections::HashMap;
    use std::iter;
    use ndarray::{Array, Array1, Array2, Axis, Slice};
    use serde_json::value::Index;
    use crate::parameters::srk::{SRKParameters,RawSRK};
    use crate::parameters::asc::{ASCParameters, SchemeType};
    use crate::eos::eos::{ASCeos, CPAeos, SRKeos};
    use crate::core::core_equations::{ASCEquations, SRKEquations};
    use crate::stability::{phi_phi_guess, tpd, watcon};
    use crate::tools::newton;
    use ndarray::s;
    use crate::Instant;

    // #[test]
    fn test_srk(){
        // PARÂMETRIZAÇÃO:
        
        let vx = Array::from_vec(vec![0.5,0.5]);
        let vb_vec = vec![ 0.01556*1e-3  ,  0.0272*1e-3  ];
        let va0_vec = vec![ 2.2519*1e5*1e-6 , 3.5079 *1e5*1e-6 ];
        let vkappa_vec = vec![0.6108 ,   0.76 ];
        let vtc_vec = vec![  647.140  ,   304.12   ] ;

        let mut map: HashMap<(usize,usize), f64> = HashMap::new();

        let tup = ( (0,1) , -0.147);
        let res = map.insert(tup.0, tup.1);        
        let ncomp:usize   = 2;

        // Inicializa classe dos parametros
        let raw = RawSRK{
            va0:va0_vec,vb:vb_vec,vtc:vtc_vec,vkappa:vkappa_vec,kij:map
        };

        let p = SRKParameters::from_raw(ncomp, raw);
        // Inic. eos
        let rho = 200.0;
        let t = 300.0;
        let srk = SRKeos{
            parameters:p,
        };
        // CÁLCULOS:
        // bmix:
        let bmix = srk.bmix(&vx);
        // aSoave:
        let soave = srk.asoave(t);
        // aij matrix:
        let aij = SRKEquations::calc_sqrt_aij_matrix(t,&srk.parameters);
        // amix:
        let amix = SRKEquations::calc_amix(t,&vx,&srk.parameters);
        // dadn:
        let dadn = SRKEquations::calc_dadn(t,&vx,&srk.parameters);
        // phi:
        let phi = srk.phi(t,rho,&vx);
        // pres:
        let pres = srk.pressure(t,rho,&vx);

        println!("SOAVE = {soave}");
        println!("bMIX= {bmix}");
        println!("aijMAT= {aij}");
        println!("aMIX= {amix}");
        println!("dadn= {dadn}");
        println!("phi= {phi}");
        println!("p = {pres}")

    }


    // #[test]
    fn test_cpa(){
        let before = Instant::now();

        // let comps = vec!["water","co2"];
        let comps = vec!["co2","ch4","water"];
        let ncomp = comps.len();
        let cpa = CPAeos::new(comps);

        let rho = 200.0;
        let t= 300.0;
        let p = 700e5;
        // let vx = Array1::from_elem(ncomp,(1 as f64/ncomp as f64) );
        let vx = Array1::from_elem(ncomp,(0.5) );

        // println!("Associative comps = {:#?}",&cpa.asc.parameters.associative_components);
        // let (pcpa,phicpa)= (cpa.pressure(t, rho, &vx)*1e-5,cpa.phi(t, rho, &vx));

        // println!("Pcpa = {pcpa} bar, phicpa = {phicpa}")

        // Testando perfomance

        // println!("{}",&vx);
        // let volume = cpa.volume(t, None, p, &vx, "Liquid");
        // println!("Volume = {}",volume);
        // println!("Tempo de execução: {:.2?}", before.elapsed());

        // Antes de modificar
        // CO2,H20,CH4
        // Volume = 0.00003672465799580874
        // Tempo de execução: 5.81ms


        // 4.616872110461483 [0.83922743 1.02740494]

        let mut y_dry_gas:Array1<f64> = Array1::from_elem(ncomp, 1./(ncomp-1)as f64);
        let mut y_dry_gas:Array1<f64> = Array1::from_elem(ncomp, 1./(ncomp-1)as f64);
        y_dry_gas[ncomp-1] = 0.0; 

        println!("ydg ch4 e co2 = {}",&y_dry_gas);

        let watcres = watcon(&cpa, t, p, &y_dry_gas);

        // println!("{}",watcres-0.0038148252909579548);
        println!("{}",watcres*1e6);

        // 0.77s


    }

    // #[test]
    fn test_json_read(){
        
        // ok
        let comps_names = vec!["co2","water"];
        let ncomp = comps_names.len();

        let eos = CPAeos::new(comps_names);
        // let rho = 200.0;
        // let t= 300.0;
        // let vx = Array1::from_elem(ncomp,0.5);
        
        // let (pcpa,phicpa)= (eos.pressure(t, rho, &vx)*1e-5,eos.phi(t, rho, &vx));

        // println!("Pcpa = {pcpa} bar, phicpa = {phicpa}");

        // let v = vec![Array1::from_vec(vec![1,2,3])];

    }
    
    // #[test]
    fn test_cpa_volume_solver(){


        let comps_names = vec!["co2","water"];
        let ncomp = comps_names.len();

        let eos = CPAeos::new(comps_names);
        // let rho = 200.0;
        let t= 300.0;
        let vx = Array1::from_elem(ncomp,0.5);

        let p = 600e5;
        let phase = "Liquid";
        let vp = Array1::linspace(1e5, 100e5, 1);

        // Tentando otimizar:

        // rho = 200, t=300K , p =  vp = Array1::linspace(1e5, 100e5, 1);
        // 1. 0.02s
        // como deu 0.02s com 6 iterações, então 300 iterações gastam 1 segundo
        // dentro do solver, a cada iteração, são feitos 3 calculos pra xassoc
        // que levam 70 iterações; logo, no total, 1260 iterações
        // 
        // o que posso fazer para otimizar seria utilizar o valor da primeira iteração, para 
        // reduzir o tempo das iterações seguintes 


        // itaração dentro do solver=1
        // Iterações Xassoc = 78
        // Iterações Xassoc = 78
        // Iterações Xassoc = 78
        // itaração dentro do solver=2
        // Iterações Xassoc = 67
        // Iterações Xassoc = 67
        // Iterações Xassoc = 67
        // itaração dentro do solver=3
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // itaração dentro do solver=4
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // itaração dentro do solver=5
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // itaração dentro do solver=6
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // Iterações Xassoc = 64
        // VolumeResult {
        //     x: 2.7596162773194837e-5,
        //     it: 6,
        // }

        // ta dando em media 0.01s pra calcular 1 volume
        // para cada iteração ( supondo 7 por eexemplo),
        // e que supondo que cada iteração gasta 30 em média
        // pra convergir xassoc,
        // pra calculos mais extensos, pode aumentar (o quanto eu nao sei)
        // uma solução é, pra cada iteração no metodo de volume,
        // utilizar o xassoc da iteração passada,
        // pra que haja uma convergencia mais rapida
        // let res = eos.volume(t, None,p, &vx,&phase);
        //     match res {
        //         Ok(val) => println!("{:#?}",val),
        //         Err(e)=> println!("{:#?}",e)
        //     }
        // }
    
        // for p in vp{

        // let res = eos.volume(t, p, &vx,&phase);
        //     match res {
        //         Ok(val) => println!("{:#?}",val),
        //         Err(e)=> println!("{:#?}",e)
        //     }


        }

    #[test]
    fn test_sat_with_more_comps(){

        let comps_names = vec!["H2S","N2", "CO2","CH4","C2","C3","C4",
           "iC4","C5","iC5","C6","C7","C8","C9" ,"C10","water"];

        // let comps_names = vec!["CO2","water"];
        let ncomp = comps_names.len();

        let eos = CPAeos::new(comps_names);
        // let rho = 200.0;
        let t= 300.;
        // let vy = Array1::from_elem(ncomp,1/ncomp);

        // TPD
        let p = 700e5;
        // let mother_phase = "Vapor";
        // let wguess= 0.9;
        // let it_max = Some(100);
        // let tol = Some(1e-6);

        // let val = tpd(&eos, t, p, wguess, vy, mother_phase, it_max, tol);

        // println!("{:#?}",val);

        
        
        // phi-phi
        // dbg!(&eos.asc.parameters.associative_components);

        let mut y_dry_gas:Array1<f64> = Array1::from_elem(ncomp, 1e-16 );

        y_dry_gas[ncomp-1] = 0.0; 
        y_dry_gas[2] = 0.5; 
        y_dry_gas[3] = 0.5; 

        // // let phi_phi_res = phi_phi_guess(&eos, t, p, &y_dry_gas);

        // // println!("{:#?}",phi_phi_res);
        // dbg!(&y_dry_gas);
        let watcon_result = watcon(&eos, t, p, &y_dry_gas);


        // ta convergindo, mas o erro ta proximo (tentar aumentar tol)

        // println!("{}",watcon_result-0.0016066701899477218);
        // println!("{}",((watcon_result-0.001897154430418538).abs()/0.001897154430418538 )*100.);
        println!("{}",(watcon_result*1e6));

        // 0.0030140686156676787 (co2 e agua )

        // 0.0016066701899477218
        // ao adicioar ch4, o erro deixou de ser da orde de 1-9 pra 1e-6;
        // o tempo de iteração aumentou aproximadamente 1s

        // 0.001878775304104048


    }
        //

    fn rascunho(){

        // let iter = iter
        // let map: HashMap<, _> = HashMap::new();

        for i in 0..1{

        }
    }

}