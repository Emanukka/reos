use super::R;

pub fn dimensionless_delta_jl(t:f64,epsilon:f64,kappa:f64)->f64{
    (epsilon / R / t).exp_m1() * kappa
}

pub fn elliot_dimensionless_delta_jl(
    t:f64,
    epsilon_j:f64,
    kappa_j:f64,
    epsilon_l:f64,
    kappa_l:f64,
    )->f64{
    let dj = dimensionless_delta_jl(t, epsilon_j, kappa_j);
    let dl = dimensionless_delta_jl(t, epsilon_l, kappa_l);
    
    (dj * dl).sqrt()
}

pub fn cr1_factor_ik(f_ii:f64,f_kk:f64)->f64{
    0.5 * (f_ii + f_kk)
}
pub fn elliot_factor_ik(f_ii:f64,f_kk:f64)->f64{
    (f_ii * f_kk).sqrt()
}

pub fn cr1_association_strength_jl(
    t:f64,
    f_ii:f64,
    f_kk:f64,
    epsilon:f64,
    kappa:f64)->f64 {

        let fik = cr1_factor_ik(f_ii, f_kk);
        let djl = dimensionless_delta_jl(t, epsilon, kappa);

        fik * djl
        
}

pub fn ecr_association_strength_jl(
    t:f64,
    f_ii:f64,
    f_kk:f64,
    epsilon_j:f64,
    epsilon_l:f64,
    kappa_j:f64,
    kappa_l:f64)->f64 {

        let fik = elliot_factor_ik(f_ii, f_kk);
        let djl = elliot_dimensionless_delta_jl(t, epsilon_j, kappa_j, epsilon_l, kappa_l);

        fik * djl

}

