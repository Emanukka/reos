mod pure;
mod mixture;
mod from_json;
mod parameters;

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
    pub fn hexane_pr78_twu91() -> Cubic {
        // species;Mw;Tc;Pc;acentricfactor;L;M;N;v_shift
// 
// hexane~|~n-hexane;86.17536;507.6;3.025e6;0.301261;0.287261195704866;0.834050543996051;2.01990552851337;8.1627905356359e-7

        let mr = CubicPureRecord::classic(507.60, 30.25 * 1e5, AlphaRecord::Twu91 { l: 0.287261195704866, m:0.834050543996051, n:2.01990552851337 }, Some(8.1627905356359e-7));

        let pr = PureRecord::new(86.17848, "hexane", mr);
        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);
        
        let p = CubicParameters::new(vec![pr], vec![], options).unwrap();
        Cubic::from(p)
    }

    pub fn hexane_octane_tcpr78_twu91() -> Cubic {
        
        let mr1 = CubicPureRecord::classic(507.60, 30.25 * 1e5, AlphaRecord::Twu91 { l: 0.287261195704866, m:0.834050543996051, n:2.01990552851337 }, Some(8.1627905356359e-7));
        // octane~|~n-octane;114.22852;568.7;2.49e6;0.399552;0.357003591909609;0.817165588467549;1.99996173995883;6.42842482682732e-6

        let mr2 = CubicPureRecord::classic(568.7, 2.49e6, AlphaRecord::Twu91 { l: 0.357003591909609, m:0.817165588467549, n:1.99996173995883 }, Some(6.42842482682732e-6));

        let pr1 = PureRecord::new(1., "hexane", mr1);
        let pr2 = PureRecord::new(0., "octane", mr2);

        // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let options = CubicOptions::classic(CubicModelOption::PR78, AlphaOption::Twu91);
        
        let p = CubicParameters::new(vec![pr1, pr2], vec![], options).unwrap();
        
        Cubic::from(p)

    }
    pub fn oil_iuri()-> Cubic{

        // let mr1 = CubicPureRecord::regressed_soave(0.12277, 0.0145e-3, 647.14, 0.6736, None);
        // let mr2 = CubicPureRecord::regressed_soave(0.35079, 0.0272e-3, 304.12, 0.7602, None);

        // let pr1 = PureRecord::new(18.02, "water", mr1);
        // let pr2 = PureRecord::new(0., "carbon dioxide", mr2);
        
        // let mbr1 = CubicBinaryRecord{kij:Kij{a:0.01, b: 0.0005}};
        // // println!("{}",mbr1);
        // let br1 = BinaryRecord::new(mbr1, "water", "carbon dioxide");
        // // let options = CubicOptions::new(CubicModelOption::PR78, Alpha::soave(), MixingRule::default());
        let names = vec!["asph", "resin", "c13_c19", "c6_c12", "c2_c5", "c1", "co2"];
        let ppath = "../../../scaling/dl/dl/oil_scaling.json";
        let options = CubicOptions::classic_soave(CubicModelOption::PR78);

        let p = CubicParameters::from_json(&names, &ppath, None, options).unwrap();
        // println!("{}",p);

        Cubic::from(p)
        // let p = CubicParameters::new(vec![pr1, pr2], vec![br1], options).unwrap();
        // println!("{}",p);
        // 
        // unimplemented!()
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