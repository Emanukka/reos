use std::error::Error;
use std::env;

fn build_file(name: String, json: &str) -> Result<(), Box< dyn Error>> {

    let mut current = env::current_dir()?;
    current.push(name);

    let _ = std::fs::write(current, json)?;

    Ok(())

}

pub fn to_json_vec<S:serde::Serialize>(name:&str, vec:Vec<S>, build:bool) -> Result<String, Box<dyn Error>> {
    
    let json = serde_json::to_string_pretty(&vec)?;
    let name = name.to_string() + ".json";
    if build{
        build_file(name, &json)?;
    }
    Ok(json)
}

#[cfg(test)]
mod tests {
    use crate::{models::cpa::parameters::readyto::{water4c,water4c_co2}};


    #[test]
    fn water_to_json(){

        let r1 = water4c();

        let json = super::to_json_vec("water_record", vec![r1], false).unwrap();

        println!("{}",json);
    }

    #[test]
    fn binary_to_json(){

        let r1: crate::parameters::BinaryRecord<crate::models::cpa::parameters::CPABinaryRecord> = water4c_co2();

        let json = super::to_json_vec("binary_record", vec![r1], false).unwrap();

        println!("{}",json);
    }
}
