use std::fmt;

pub struct NewtonResult{
    pub x:f64,
    pub it:i32,
    pub it_backtracking:i32
}

impl fmt::Display for NewtonResult {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "x = {}, it = {}, it backtracking = {}",self.x,self.it,self.it_backtracking)
    }
}
#[derive(Debug)]
pub enum ErrorAtFindRoot {
    NaNValue(i32),
    DerivativeIsZero(i32),
    MaxIterations
}
// 
impl std::error::Error for ErrorAtFindRoot {}

impl fmt::Display for ErrorAtFindRoot {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::MaxIterations=> write!(f, "Max number of iterations was reached"),
            Self::NaNValue(i)=> write!(f, "NaN value in result of root's find at {} iteration",i),
            Self::DerivativeIsZero(i)=> write!(f, "Derivative result was 0 at {} iterations.
             Then, it was not possible to perfome any calculation",i)
                
            }
            
        }
    }



pub fn  newton<F>(f:F,x:f64,tol:Option<f64>,it_max:Option<i32>)->Result<NewtonResult,ErrorAtFindRoot>
    where F: Fn(f64)->f64 {


        if x.is_nan(){
            return Err(ErrorAtFindRoot::NaNValue(0));
        }
        let it_max: i32 = if let Some(val) = it_max {val} else {100};
        let tol = if let Some(val) = tol {val} else {1e-6};

        let fn_dfdx = dfdx(&f);
        let mut res: f64 = 1.0;

        let mut x0:f64;
        let mut x1:f64 = x;
        let mut t:f64;
        let mut h:f64;
        let mut ft:f64;
        let mut alpha:f64;
        let mut it: i32 = 0;
        let mut d :f64 = 1.;
        let mut backtrack_it: i32 = 0;

        while (res.abs()>tol) & (it<it_max){


            // x1 = x0 - f(x0)/dfdx(x0)
            alpha=1.;

            x0 = x1;
            // dbg!(res,x0,it);

            res = f(x0);
            d = fn_dfdx(x0);
            h = - res/d;
            t = x0 + h;

            ft = f(t);

            while (ft.abs())>(res.abs()) {

                backtrack_it+=1;
                alpha=alpha/2.;
                t = x0 + alpha*h;
                ft = f(t);
            }
            // dbg!(res.abs()>tol);

            x1 = t; 
            it+=1;
            // dbg!(res,x0,it);


        }

        // println!("it={}, backtrack it={}",it,backtrack_it);
        if x1.is_nan(){

            if d==0.0{
                Err(ErrorAtFindRoot::DerivativeIsZero(it))
            }else{
                Err(ErrorAtFindRoot::NaNValue(it))
            }
        }else if it==it_max {

            Err(ErrorAtFindRoot::MaxIterations)
        }else {
            Ok(NewtonResult{x:x1,it:it,it_backtracking:backtrack_it})
        }

    }

fn dfdx<F>(f:F)->impl Fn(f64)-> f64
    where F: Fn(f64)-> f64 {

        move |x:f64|( f(x+1e-5) - f(x-1e-5) )/(2.*1e-5)

    }