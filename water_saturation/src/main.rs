use std::sync::Arc;

use reeos::Array1;
use water_saturation::{co2_water, WaterSaturation};



fn read_input(head:&str)->String{
        // Cria uma variável mutável para armazenar o input
    let mut entrada = String::new();


    println!("{}",head);
    // Lê da entrada padrão (stdin) e coloca na variável `entrada`
    std::io::stdin()
        .read_line(&mut entrada) // lê uma linha (inclui o '\n')
        .expect("Erro ao ler a entrada");

    // Remove o '\n' ou '\r\n' do final (opcional, mas comum)
    let entrada = entrada.trim();
    entrada.to_string()
}


fn main() {
    let texto = 
    r#"

    Calculo de solubilidade de água em gás natural (apenas CO2)

    Opções:
    [1]- Calcular yW, dados (T,P)
    [2]- Calcular Vetor_yW, dados (T,vetor_P)
    [9]- Sair
    "#;

    let eos=Arc::new(co2_water());

    let wsat=WaterSaturation{eos};


    println!("{}",texto);
    let mut ok=true;

    let bar_to_pa=1e5;
    let c_to_k=273.15;

    let mut p=bar_to_pa;
    let mut t=c_to_k;

    while ok{

        let input=read_input("Digite:");

        match input.as_str(){
            
            "1"=>{
                p=read_input("P (bar)").parse::<f64>().unwrap()*bar_to_pa;
                t=read_input("T (°C)").parse::<f64>().unwrap()+c_to_k;

                let y_dry_gas=Array1::from_vec(vec![1.0,0.0]);

                let result=wsat.watcon(t, p, &y_dry_gas);

                match result {
                    
                    Ok(y)=>{
                        println!("Saturação de água (ppm): {}",y*1e6)
                    }
                    Err(e)=>{
                        println!("Erro de cálculo: {}",e)
                    }
                }
            }

            "2"=>{

            }

            "9"=>{
                ok=false
            }
            _=>{
                println!("input errado!")
            }

        }

    }
    // println!("Você digitou: {}", entrada);
}
