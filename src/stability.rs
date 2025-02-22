
use ndarray::Array1;
use crate::{eos::eos::{CPAeos, Phase}, tools::newton};
const XW: f64 = 0.999999;
const YWGUESS:f64 = 2000e-6;


// define as funções necessárias para realizar o cálculo de saturação de água em um gás (cuja informação será sua composição em base seca)

pub fn tpd(eos:&CPAeos,t:f64,p:f64,wguess:f64,vy:Array1<f64>, it_max:Option<i32>,tol:Option<f64>)->f64{
    // wguess: chute da fração de água na fase líquida; 
    // ~0.99

    // vy: vetor composição da fase mãe (no caso do cálculo da saturação de água,
    //     a fase mãe será a vapor, e a filha (incipiente), liquida);
    // let (vapor,liquid) = ("VAPOR".to_string(),"LIQUID".to_string());

    // let daughter_phase: &str = match mother_phase.to_uppercase() {
    //     vapor => "LIQUID",
    //     liquid=> "VAPOR"
    // };


    let it_max = if let Some(val) = it_max {val} else {100};
    let tol = if let Some(val) = tol {val} else {1e-6};
    let ncomp = eos.asc.parameters.ncomp;

    // let one_minus_n= 1.0 - ncomp as f64;
    let water_index = ncomp-1;
    let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-wguess)/(ncomp-1) as f64  ) ;

    vx[water_index]= wguess;

    // constantes
    let rho_mother_phase = 1./eos.volume(t, None,p, &vy, "VAPOR");
    let phiy = eos.phi(t, rho_mother_phase, &vy); 

    let hy = (&vy*phiy).ln();
    // let mut vw_old:&Array1<f64>;
    let mut vw_old: Array1<f64>;
    let mut vw:Array1<f64> = vx.clone();
    let mut phi_x = Array1::<f64>::zeros(ncomp);
    let mut res = 1.0;
    let mut it = 0;

    // println!("{:p}",&vw);
    // dbg!(&hy,&vy,&vx);
    let bm = eos.srk.bmix(&vx);
    let rho_max = 1./bm;
    let mut rho_x = if let Phase::Liquid { initial_guess } = CPAeos::which_phase_is(t, p, bm, "LIQUID") {
        initial_guess*rho_max
    }else{

        // println!("erro");
        1.
    };
    // println!("rhox = {rho_x}");

    while (res>tol) & (it<it_max) {

    // println!("TPD it = {}",it);

        
        // vw_old = &vw;

        // vw_old = &vw;
        vw_old = vw;
        // println!("{:p}",vw);


        // usar rho antigo?
        // println!("rhox = {rho_x}");

        // rho_x = 1./eos.volume(t, Some(rho_x),p, &vx, daughter_phase).unwrap().x;
        rho_x = 1./eos.volume(t, Some(rho_x),p, &vx, "LIQUID");
        phi_x = eos.phi(t,rho_x,&vx);

        vw = ( &hy-phi_x.ln() ).exp();
        
        // se eu fizesse com let vw =..., o valor o copiado no inicio do proximo loop
        // fosse do vw do escopo maior; então o loop nao iria prosseguir

        // let vw= ... cria um novo endereço na memoria, mas, ao terminar o loop, esse valor é sempre apagado
        // e entao o valor utilizado é o do escopo maior
        // println!("{:p}",&vw);


        // analisar dps o comportamento em relação a referencias (de vw_old e vw)

        // dbg!(&vw_old,&vw,"second");

        vx = &vw/vw.sum();

        res = (&vw - vw_old).pow2().sum().sqrt();
        // dbg!(&rho_x);
        // dbg!(it,res,rho_x);
        it+=1
    }
    
    let hx = (vx*phi_x).ln();
    // dbg!(it);

    let dg = (vw*(hx-hy)).sum();

    // pedir pra retornar um TPDResult, contendo dg,x,it
    dg

}

pub fn phi_phi_guess(eos:&CPAeos,t:f64,p:f64,y_dry_gas:&Array1<f64>)->f64{

    // nunca foi necessário fazer um antoine pra estimar ywguess, mas talvez
    // seja util em algum momento

    
    let ncomp = eos.asc.parameters.ncomp;
    let water_index = ncomp-1;


    let res = |log10_yw0:f64|{

        let yw0 = 10.0_f64.powf(log10_yw0);

        let mut vy_wet = (y_dry_gas/y_dry_gas.sum())*(1.-yw0);
        vy_wet[water_index] = yw0;

        let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;
        vx[water_index]= XW;


        // fase vapor (mistura)

        let rho_vap = 1./eos.volume(t, None,p, &vy_wet, "Vapor"); 
        // let rho_vap = 1./eos.volume(t, None,p, &vy_wet, "Vapor").x; 
        let phi_water_vap = eos.phi(t, rho_vap, &vy_wet)[water_index];

        // fase liq (~agua pura)
        // let rho_liq = 1./eos.volume(t, None,p, &vx, "Liquid").unwrap().x; 
        let rho_liq = 1./eos.volume(t, None,p, &vx, "Liquid"); 
        let phi_water_liq = eos.phi(t, rho_liq, &vx)[water_index];
        
        ((yw0*phi_water_vap)/(phi_water_liq)).ln()
        // yw0.ln() + phi_water_vap.ln() - phi_water_liq.ln()

    };
    
    let log10yw0 = YWGUESS.log10();
    let result = newton(res, log10yw0, None, None).unwrap();

    // println!("{result}");
    10.0_f64.powf(result.x)


}

pub fn watcon(eos:&CPAeos,t:f64,p:f64,y_dry_gas:&Array1<f64>)->f64{

    let log10_ywguess_from_phiphi = phi_phi_guess(eos, t, p, y_dry_gas).log10();
    let ncomp = eos.asc.parameters.ncomp;
    let water_index = ncomp-1;

    let res = |log10_yw0:f64|{
        
        let yw0 = 10.0_f64.powf(log10_yw0);
        let mut vy_wet = (y_dry_gas/ y_dry_gas.sum())*(1.-yw0);
        vy_wet[water_index] = yw0;

        let mut vx = Array1::<f64>::ones(ncomp) * ( (1.-XW)/(ncomp-1) as f64  ) ;
        vx[water_index]= XW;

        tpd(eos,t,p,XW,vy_wet,None,None)

    };

    let newton_result = newton(res, log10_ywguess_from_phiphi, None, None).unwrap();

    // println!("{newton_result}");
    let result = 10.0_f64.powf(newton_result.x);
    // println!("{result}")
    result
    


}