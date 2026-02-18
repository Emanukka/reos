mod pure;
mod mixture;
mod json;
use crate::{models::cubic::{Cubic, alpha::Alpha, mixing_rule::MixingRule, models::PR78, options::CubicOptions, parameters::{CubicBinaryRecord, CubicParameters, CubicPureRecord}}, parameters::{BinaryRecord, Parameters, PureRecord}};

#[cfg(test)]
pub mod recipes {

    use super::*;
    pub fn water() -> Cubic {
        // tc = 647.1, pc = 220.55e5, w = 0.345
        let mr = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let pr = PureRecord::new(18.02, "water", mr);
        // let options = CubicOptions::new(PR78.into(), Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic_soave(PR78.into());
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from_parameters(p)

        
    }

    pub fn water_co2() -> Cubic {
        
        let mr1 = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let mr2 = CubicPureRecord::classic_soave(304.2, 73.83e5, 0.224, None);

        let pr1 = PureRecord::new(18.02, "water", mr1);
        let pr2 = PureRecord::new(0., "carbon dioxide", mr2);

        // let options = CubicOptions::new(PR78.into(), Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic_soave(PR78.into());
        
        let p = CubicParameters::new(vec![pr1, pr2], vec![], options).unwrap();
        
        Cubic::from_parameters(p)

    }
    pub fn water_co2_bip() -> Cubic {
        
        let mr1 = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let mr2 = CubicPureRecord::classic_soave(304.2, 73.83e5, 0.224, None);

        let pr1 = PureRecord::new(18.02, "water", mr1);
        let pr2 = PureRecord::new(0., "carbon dioxide", mr2);

        let mbr1 = CubicBinaryRecord{kij:0.1, lij:0.05};
        // println!("{}",mbr1);
        let br1 = BinaryRecord::new(mbr1, "water", "carbon dioxide");
        // let options = CubicOptions::new(PR78.into(), Alpha::soave(), MixingRule::default());
        
        let options = CubicOptions::classic_soave(PR78.into());

        let p = CubicParameters::new(vec![pr1, pr2], vec![br1], options).unwrap();
        // println!("{}",p);
        // 
        Cubic::from_parameters(p)

    }
    }
    pub fn nhexane_dehlouz() -> Cubic {
        
        unimplemented!()
        // let pr = CubicPureRecord::Twu91 { 
        //     tc: 507.60, 
        //     pc: 30.25 * 1e5, 
        //     l: 0.28726, 
        //     n: 2.01991, 
        //     m: 0.83405, 
        //     volt: Some(0.81628 / 1e6)
        // };
        // let pr = PureRecord::new(86.17848, "n-hexane".into(), pr);
        // Cubic::from_parameters(CubicParameters::new(vec![pr], vec![], PR78.into()))
    }
    pub fn water_dehlouz() -> Cubic {
        
        unimplemented!()
        // let pr = CubicPureRecord::Twu91 { 
        //     tc: 507.60, 
        //     pc: 220.60e5, 
        //     l: 0.38720, 
        //     n: 1.96692, 
        //     m: 0.87198, 
        //     volt: Some(5.27106e-6)
        // };
        // let pr = PureRecord::new(18.01528, "water".into(), pr);
        // Cubic::from_parameters(CubicParameters::new(vec![pr], vec![], PR78.into()))
    }


#[test]
fn alpha_function_unmatch_records(){

    let mr = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let pr = PureRecord::new(18.02, "water", mr);

        let options = CubicOptions::classic(PR78.into(), Alpha::twu91());
        
        let boxx = CubicParameters::new(vec![pr], vec![], options).unwrap_err();
        if let Ok(e) = boxx.downcast::<super::alpha::AlphaError>(){

            assert_eq!(e, Box::new(super::alpha::AlphaError::AssociatedRecord("water".into())))

        }

        // Cubic::from_parameters(p)

}