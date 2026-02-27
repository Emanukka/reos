use crate::models::IDEAL_GAS_CONST as R;
use crate::models::associative::parameters::AssociativeParameters;
use crate::state::eos::{EosError, EosResult};

use ndarray::{Array1, Array2};


pub mod sites;
pub mod sites_interaction;
pub mod parameters;
pub mod strength;
pub mod combining_rule;

// #[cfg(test)]
// mod tests;

#[derive(Clone)]
pub struct Associative{
  pub parameters:AssociativeParameters,
}



impl Associative {

    pub fn x_ab_analytic(m1:f64,m2:f64,mult1:f64,mult2:f64,k:f64)->Array1<f64>{

        
        let x2 = (2.0 * m1 * m2) / (
                    m1 * m2
                    + k * m1 * mult1
                    + (4.0 * k * m2 * m1 * m1 * mult1
                    + (k * m2 * mult2 - k * m1 * mult1 + m1 * m2).powi(2)
                    ).sqrt()
                    -k * m2 * mult2
                );
        let x1 = m1 / (m1 + k * mult2 * x2);

        Array1::from_vec(vec![x1, x2])
    }

    pub fn x_cc_analytic(m:f64, mult:f64, k:f64)->Array1<f64>{
        
        let x1 = 2.0 * m / (m + (4. * k * mult * m + m.powi(2)).sqrt());
        Array1::from_vec(vec![x1])
    }   

    pub fn unbonded_sites_fraction(&self,x:&Array1<f64>, k:&Array2<f64>)-> Array1<f64>{
        
        let p = &self.parameters;
        let m = self.sites_mole_frac(x);
        let sites = &p.sites;

        
        match (p.na * p.nb, p.nc) {
            
            (1, 0) => Self::x_ab_analytic(m[0], m[1], sites[0].mul, sites[1].mul, k[(0,1)] ),

            (0, 1) => Self::x_cc_analytic(m[0],sites[0].mul, k[(0,0)]), 
            
            (_, _) => self.x_tan(&m, &k).unwrap_or(Array1::ones(sites.len()) * f64::NAN),

        }
        
    }

    pub fn sites_mole_frac(&self, x:&Array1<f64> )->Array1<f64>{
        
        let sites = &self.parameters.sites;
        Array1::from_shape_fn(sites.len(), |j| x[sites[j].owner])  

    }

    pub fn h(&self, x:&Array1<f64>, unbonded:&Array1<f64> )->f64{
        
        let mut h = 0.0;
        let sites = &self.parameters.sites;
        
        for sj in sites {
            
            let i = sj.owner;
            let j = sj.idx;
            let m = sj.mul;


            h += x[i] * (1. - unbonded[j]) * m
        }
        h
    }

    pub fn association_constants(&self, t:f64, rho:f64, x:&Array1<f64>, volf:&Array1<f64>)->Array2<f64>{

        let sites = &self.parameters.sites;
        let interactions = &self.parameters.interactions;
        let s = sites.len();
        let mut kmat = Array2::zeros((s,s));

        for interaction in interactions{
            
            let j = interaction.site_j;
            let l = interaction.site_l;

            let i = sites[j].owner;
            let k = sites[l].owner;

            let f_ii = volf[i];
            let f_kk = volf[k];

            kmat[(j,l)] = interaction.association_strength_jl(t, f_ii, f_kk, sites) * x[i] * x[k] * rho;
            kmat[(l,j)] = kmat[(j,l)]
            
        }
        kmat

    }

    pub fn dk_dt(&self, t:f64, rho:f64, x:&Array1<f64>, volf:&Array1<f64>)->Array2<f64>{

        let sites = &self.parameters.sites;
        let interactions = &self.parameters.interactions;
        let s = sites.len();
        let mut dk_dt = Array2::zeros((s,s));

        for interaction in interactions{
            
            let j = interaction.site_j;
            let l = interaction.site_l;

            let i = sites[j].owner;
            let k = sites[l].owner;

            let f_ii = volf[i];
            let f_kk = volf[k];

            dk_dt[(j,l)] = interaction.association_strength_jl_dt(t, f_ii, f_kk, sites) * x[i] * x[k] * rho;
            dk_dt[(l,j)] = dk_dt[(j,l)]
            
        }
        dk_dt

    }

    pub fn x_tan(&self, m:&Array1<f64>, k:&Array2<f64>) -> EosResult<Array1<f64>>{

        let s = m.len();
        let sites = &self.parameters.sites;
        let mult = Array1::from_shape_fn(s, |j| sites[j].mul);
        // let mult= &self.parameters.multiplicity;

        let mut unbonded =  Array1::from_elem(s,0.2);

        let tol = 1e-12;
        let max_iter = 1000;  

        for i in 0..max_iter {
            if self.tan_step(
                &mult,
                &mut unbonded,
                &m,
                k,
                tol,
            )? {
                break;
            }
            if i == max_iter - 1 {
                return Err(EosError::NotConverged("Unbonded fraction".into()));
            }
        }
        
        Ok(unbonded)


    }

    pub fn tan_step(
        &self,
        mult:&Array1<f64>,
        unbonded:&mut Array1<f64>,
        m:&Array1<f64>,
        k:&Array2<f64>,
        tol:f64)->EosResult<bool>{

        
        let xm =  &(& * unbonded * mult);
        let s = xm.len();

        let mut g = Array1::zeros(s);

        for j in 0..s {
            
            let dot = k.row(j).dot(xm);

            let correc = m[j] / (m[j] + dot);

            unbonded[j] = 0.75 * correc + 0.25 * unbonded[j];
            // xm[j] = unbonded[j] * mult[j];
            g[j] = m[j] * (1. / unbonded[j] - 1.) - dot;
        }

        let err = g.iter().map(|&g| g.powi(2)).sum::<f64>().sqrt();

        Ok( err < tol)

    }
}

// derivatives
impl Associative {

    pub fn df_dv(&self, h:f64, d:f64, dlng_drho:f64)->f64 {

        d * (1. + d * dlng_drho) * h / 2.
        
    }

    pub fn df_dn(&self, h:f64, ndlng_dni:&Array1<f64>, unbonded:&Array1<f64>, )->Array1<f64>{
        
        let sites = &self.parameters.sites;
        let mut df_dn = - 0.5 * h * ndlng_dni;

        for sj in sites{

            let i = sj.owner;
            let j = sj.idx;
            let m = sj.mul;

            df_dn[i] += unbonded[j].ln() * m 
        }   
        df_dn
        


    }
    
    pub fn helmholtz(&self, x:&Array1<f64>, unbonded:&Array1<f64> )->f64{

        let sites = &self.parameters.sites;
        
        sites.iter().fold(0.0, |acc, sj| {

            let i = sj.owner;
            let j = sj.idx;
            let m = sj.mul;

            acc + x[i] * (unbonded[j].ln() - 0.5 * unbonded[j] + 0.5) * m
        })

    }

    pub fn df_dt(&self, _:f64, x:&Array1<f64>, k:&Array2<f64>, dk_dt:&Array2<f64>)->f64{

        let sites = &self.parameters.sites;
        let unbonded = &self.unbonded_sites_fraction(x, k);
        let mut df_dt = 0.0;

        let s = sites.len();

        for j in 0..s{
            for l in 0..s {

                let xmj = unbonded[j] * sites[j].mul;
                let xml = unbonded[l] * sites[l].mul; 

                df_dt += dk_dt[(j, l)] * xmj * xml 

            }
        }
        df_dt *= - 0.5 ;

        df_dt


    }


    // // pub fn x_michelsen(
    // //     &self,
    // //     m:&Array1<f64>,
    // //     k:&Array2<f64>) -> EosResult<Array1<f64>>{

    // //     let mult= &self.parameters.multiplicity;
    // //     // let m= self.m(x);

    // //     let mut unbonded =  Array1::from_elem(mult.len(),0.2);

    // //     let tol = 1e-8;
    // //     let max_iter = 100;  
        
        
    // //     for _ in 0..5 {
        
    // //         let xm = & unbonded * mult;


    // //         let correc = m / (
    // //             m + k.dot(&xm)
    // //         );

    // //         unbonded.iter_mut().zip(&correc).for_each(|(x, &correc)| {

    // //             *x = 0.75 * correc + 0.25 * *x

    // //         });

    // //     }

    // //     for i in 0..max_iter {
    // //         if self.newton_step(
    // //             &mut unbonded,
    // //             &m,
    // //             k,
    // //             tol,
    // //         )? {
    // //             // dbg!(i);
    // //             break;
    // //         }
    // //         if i == max_iter - 1 {
    // //             return Err(EosError::NotConverged("Unbonded fraction".into()));
    // //         }
    // //     }
        
    // //     Ok(unbonded)


    // // }
    // // pub fn newton_step(
    // //     &self,
    // //     unbonded:&mut Array1<f64>,
    // //     m:&Array1<f64>,
    // //     k:&Array2<f64>,
    // //     tol:f64)->EosResult<bool>{
        
    // //     let s = unbonded.len();
    // //     let mult = &self.parameters.multiplicity;
    // //     // let q_old = Self::q_michelsen(mult,m, unbonded, k);
        

    // //     // let xm =  &(& * unbonded * mult);
    // //     let xm =  Zip::from(&mut *unbonded).and(mult).map_collect(|x,y| *x * y);
    // //     // let h = Self::hessian(&xm,m, unbonded,k);
    // //     // let g = Self::grad(&xm,m, unbonded, k);
    // //     let mut h= -k.clone();
    // //     let mut g= Array1::zeros(s);
            
    // //     for j in 0..s {
            
    // //         // let xm = unbonded[j] * mult[j];

    // //         let dot = k.row(j).dot(&xm);

    // //         g[j] = m[j] * (1. / unbonded[j] - 1.) - dot;


    // //         h[(j, j)] -=  ( m[j] + dot ) / unbonded[j];
    // //         // h[(j, j)] -=  m[j]  / unbonded[j].powi(2);
    // //     }
    // //     let err = g.norm();
    // //     // h
    // //     // } 
    // //     // let hinv = h.inv()?;
    // //     // let mut dx = h.solve(&(-&g))?;

    // //     // println!("{}",&h);
    // //     // println!("{}",&g);
    // //     let dx= - h.solve(&g)?;
    // //     unbonded.iter_mut().zip(dx).for_each(|(x,dx)|{
            
    // //         *x += dx
    // //     });

    // //     // dbg!(&err);
    // //     // dbg!(&unbonded);
    // //     Ok( err < tol)

    // // }

    // pub fn x_michelsen(
    //     &self,
    //     mult:&DVector<f64>,
    //     m:&DVector<f64>,
    //     k:&DMatrix<f64>) -> EosResult<DVector<f64>>{

    //     let s = m.len();

    //     // let mult= self.parameters.multiplicity.as_slice().unwrap();
    //     // let mult = DVector::from_row_slice(mult);
    //     // let m= self.m(x);
    //     // let m = m.as_slice().unwrap();
    //     // let m = &DVector::from_row_slice(m);

    //     let mut unbonded =  DVector::from_element(s,0.2);

    //     let tol = 1e-8;
    //     let max_iter = 100;  
        
    //     // let k_slice = k.as_slice().unwrap();
    //     // let k = &DMatrix::from_row_slice(s,s,k_slice);
        
    //     for _ in 0..5 {
        
    //         let xm = & unbonded.component_mul(&mult);


    //         let correc = m.component_div(
    //             &(m + k * xm)
    //         );

    //         unbonded.iter_mut().zip(&correc).for_each(|(x, &correc)| {

    //             *x = 0.75 * correc + 0.25 * *x

    //         });

    //     }

    //     for i in 0..max_iter {
    //         if self.newton_step(
    //             &mult,
    //             &mut unbonded,
    //             m,
    //             k,
    //             tol,
    //         )? {
    //             // dbg!(i);
    //             break;
    //         }
    //         if i == max_iter - 1 {
    //             return Err(EosError::NotConverged("Unbonded fraction".into()));
    //         }
    //     }
        
    //     Ok(unbonded)


    // }
    // pub fn newton_step(
    //     &self,
    //     mult:&DVector<f64>,
    //     unbonded:&mut DVector<f64>,
    //     m:&DVector<f64>,
    //     k:&DMatrix<f64>,
    //     tol:f64)->EosResult<bool>{
        
    //     let s = unbonded.len();
    //     // let q_old = Self::q_michelsen(mult,m, unbonded, k);
        

    //     // let xm =  &(& * unbonded * mult);
    //     // let xm =  Zip::from(&mut *unbonded).and(mult).map_collect(|x,y| *x * y);
    //     // let xm =  Zip::from(&mut *unbonded).and(mult).map_collect(|x,y| *x * y);
    //     let xm = unbonded.component_mul(mult);

    //     // let h = Self::hessian(&xm,m, unbonded,k);
    //     // let g = Self::grad(&xm,m, unbonded, k);
    //     let mut h= -k.clone();
    //     let mut g= DVector::zeros(s);
            
    //     for j in 0..s {
            
    //         // let xm = unbonded[j] * mult[j];
    //         // dbg!(k.row(j));
    //         // dbg!(&xm);
    //         let rowj = k.row(j).transpose();

    //         // dbg!(&rowj);
    //         let dot = rowj.dot(&xm);

    //         g[j] = m[j] * (1. / unbonded[j] - 1.) - dot;


    //         h[(j, j)] -=  ( m[j] + dot ) / unbonded[j];
    //         // h[(j, j)] -=  m[j]  / unbonded[j].powi(2);
    //     }
        
    //     let err = g.norm();

    //     // h
    //     // } 
    //     // let hinv = h.inv()?;
    //     // let mut dx = h.solve(&(-&g))?;

    //     // println!("{}",&h);
    //     // println!("{}",&g);

    //     let lu = LU::new(h);
    //     // let lu = h.lu();

    //     let dx= - lu.solve(&g).unwrap();

    //     *unbonded = &*unbonded + &dx;        
    //     // unbonded.iter_mut().zip(&dx).for_each(|(x,dx)|{
            
    //     //     *x += dx

    //     // });

    //     // dbg!(&err);
    //     // dbg!(&unbonded);
    //     Ok( err < tol)

    // }

}

