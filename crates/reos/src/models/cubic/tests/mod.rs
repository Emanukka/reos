mod pure;
mod mixture;
mod from_json;

use crate::{models::cubic::{Cubic, alpha::{Alpha, AlphaOption}, mixing_rule::MixingRule, models::{CubicModelOption, PR78}, options::CubicOptions, parameters::{CubicBinaryRecord, CubicParameters, CubicPureRecord}}, parameters::{BinaryRecord, Parameters, PureRecord}};

#[cfg(test)]
pub mod recipes {

    use crate::models::cubic::{alpha::AlphaRecord, models::CubicModelOption, parameters::Kij};

    use super::*;
    pub fn water() -> Cubic {
        // tc = 647.1, pc = 220.55e5, w = 0.345
        let mr = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let pr = PureRecord::new(18.02, "water", mr);
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from(p)

        
    }
    pub fn water_vt() -> Cubic {
        // tc = 647.1, pc = 220.55e5, w = 0.345
        let mr = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, Some(0.5 / 1e6));
        let pr = PureRecord::new(18.02, "water", mr);
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from(p)

        
    }

    pub fn water_vt_jaubert() -> Cubic {
        // tc = 647.1, pc = 220.55e5, w = 0.345
        let mr = CubicPureRecord::classic(647.1, 220.55e5, AlphaRecord::Twu91 { l: 0.3865, m: 0.8720, n: 1.9693 }, Some(5.3041 / 1e6));

        let pr = PureRecord::new(18.02, "water", mr);
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from(p)

        
    }
    
    pub fn octane_vt() -> Cubic {

        let mr = CubicPureRecord::classic(568.70, 24.90 * 1e5, AlphaRecord::Twu91 { l: 0.3385, m: 0.8185, n: 2.0747 }, Some(6.4134 / 1e6));

        let pr = PureRecord::new(18.02, "water", mr);
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from(p)

    }

    pub fn water_co2() -> Cubic {
        
        let mr1 = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let mr2 = CubicPureRecord::classic_soave(304.2, 73.83e5, 0.224, None);

        let pr1 = PureRecord::new(18.02, "water", mr1);
        let pr2 = PureRecord::new(0., "carbon dioxide", mr2);

        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);
        
        let p = CubicParameters::new(vec![pr1, pr2], vec![], options).unwrap();
        
        Cubic::from(p)

    }

    pub fn water_co2_vt() -> Cubic {
        
        let mr1 = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345,Some(0.5 / 1e6));
        let mr2 = CubicPureRecord::classic_soave(304.2, 73.83e5, 0.224,Some(0.8 / 1e5) );

        let pr1 = PureRecord::new(18.02, "water", mr1);
        let pr2 = PureRecord::new(0., "carbon dioxide", mr2);

        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);
        
        let p = CubicParameters::new(vec![pr1, pr2], vec![], options).unwrap();
        
        Cubic::from(p)

    }
    pub fn water_co2_bip() -> Cubic {
        
        let mr1 = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let mr2 = CubicPureRecord::classic_soave(304.2, 73.83e5, 0.224, None);

        let pr1 = PureRecord::new(18.02, "water", mr1);
        let pr2 = PureRecord::new(0., "carbon dioxide", mr2);
        
        let mbr1 = CubicBinaryRecord{kij:Kij{a:0.01, b: 0.0005}};
        // println!("{}",mbr1);
        let br1 = BinaryRecord::new(mbr1, "water", "carbon dioxide");
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);

        let p = CubicParameters::new(vec![pr1, pr2], vec![br1], options).unwrap();
        // println!("{}",p);
        // 
        Cubic::from(p)

    }
    pub fn water_co2_bip_regressed() -> Cubic {
        
        let mr1 = CubicPureRecord::regressed_soave(0.12277, 0.0145e-3, 647.14, 0.6736, None);
        let mr2 = CubicPureRecord::regressed_soave(0.35079, 0.0272e-3, 304.12, 0.7602, None);

        let pr1 = PureRecord::new(18.02, "water", mr1);
        let pr2 = PureRecord::new(0., "carbon dioxide", mr2);
        
        let mbr1 = CubicBinaryRecord{kij:Kij{a:0.01, b: 0.0005}};
        // println!("{}",mbr1);
        let br1 = BinaryRecord::new(mbr1, "water", "carbon dioxide");
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);

        let p = CubicParameters::new(vec![pr1, pr2], vec![br1], options).unwrap();
        // println!("{}",p);
        // 
        Cubic::from(p)

    }

    pub fn nhexane_dehlouz() -> Cubic {
        
        // let pr = CubicPureRecord::Twu91 { 
        //     tc: 507.60, 
        //     pc: 30.25 * 1e5, 
        //     l: 0.28726, 
        //     n: 2.01991, 
        //     m: 0.83405, 
        //     volt: Some(0.81628 / 1e6)
        // };
        // let pr = PureRecord::new(86.17848, "n-hexane".into(), pr);
        // Cubic::from(CubicParameters::new(vec![pr], vec![], CubicModelOption::PR78))
        let mr = CubicPureRecord::classic(507.60, 30.25 * 1e5, AlphaRecord::Twu91 { l: 0.28726, m: 0.83405, n: 2.01991 }, Some(0.81628 / 1e6));

        let pr = PureRecord::new(86.17848, "hexane", mr);
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from(p)
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
        // Cubic::from(CubicParameters::new(vec![pr], vec![], CubicModelOption::PR78))
    }

}
#[test]
fn alpha_function_unmatch_records(){

    let mr = CubicPureRecord::classic_soave(647.1, 220.55e5, 0.345, None);
        let pr = PureRecord::new(18.02, "water", mr);

        let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);
        
        let boxx = CubicParameters::new(vec![pr], vec![], options).unwrap_err();
        if let Ok(e) = boxx.downcast::<super::alpha::AlphaError>(){

            assert_eq!(e, Box::new(super::alpha::AlphaError::AssociatedRecord("water".into())))

        }

        // Cubic::from(p)

}