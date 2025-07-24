use std::sync::Arc;

use reeos::Array1;

use water_saturation::data::Input;


// fn read_input(head:&str)->String{
//         // Cria uma variável mutável para armazenar o input
//     let mut entrada = String::new();


//     println!("{}",head);
//     // Lê da entrada padrão (stdin) e coloca na variável `entrada`
//     std::io::stdin()
//         .read_line(&mut entrada) // lê uma linha (inclui o '\n')
//         .expect("Erro ao ler a entrada");

//     // Remove o '\n' ou '\r\n' do final (opcional, mas comum)
//     let entrada = entrada.trim();
//     entrada.to_string()
// }


fn main() {
    let json = r#"

    {
        "composição":{

            "CO2":0.5,
            "CH4":0.5,
            "N2":0.0,
            "H2S":0.0

        },
        
        "temperatura": 25,
        "pressão":[100,200,300]

    }
    "#;

    let map:  Input= serde_json::from_str(json).unwrap();
    // let pessoas: Vec<Pessoa> = map.into_values().collect();
    let comp=map.composição.in;
    println!("{:?}", map);

}
