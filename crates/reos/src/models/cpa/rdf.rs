use std::marker::PhantomData;

use ndarray::{Array1, Array2};

#[derive(Clone)]
pub struct Rdf<T>{
  /// Co-volume of components of mixture
  pub b:Array1<f64>,
  /// 0.5 * (bi + bk) 
  pub bij:Array2<f64>,
  model: PhantomData<T>
}

impl<T:RdfModel> Rdf<T> {

    // pub fn new(b:Array1<f64>)->Self{

    //   Rdf { b, bij, model: PhantomData }
    // }
    pub fn detadrho(&self,x:&Array1<f64>)->f64{
        self.b.dot(x)/4.0
    }
    pub fn detadni(&self,rho:f64)->Array1<f64>{
        rho * &self.b/4.0
    }
    pub fn dlngdrho(&self,rho:f64,x:&Array1<f64>)->f64{
      let dlngdeta = T::dlngdeta(rho, x,&self.b);
      let detadrho = self.detadrho(x);
      dlngdeta * detadrho
    }
    pub fn ndlngdni(&self,rho:f64,x:&Array1<f64>)->Array1<f64>{
      let dlngdeta = T::dlngdeta(rho, x, &self.b);
      let detadni = self.detadni(rho);
      dlngdeta * detadni
    }
    pub fn g(&self,rho:f64,x:&Array1<f64>)->f64{
      T::g(rho, x, &self.b)
    }

}
pub trait RdfModel {
    
    fn model()->Self where Self: Sized;

    fn new<T:RdfModel>(b:Array1<f64>)->Rdf<T>{

      let n = b.len();
      let bij = Array2::from_shape_fn((n,n), |(i,j)| {
        0.5 * (b[i] + b[j])
      });    
      
      Rdf { b, bij, model: PhantomData }
    
    }

    fn eta(
      rho:f64,
      x:&Array1<f64>,
      b:&Array1<f64>)->f64{
        rho*b.dot(x)/4.0
      }
    fn dlngdeta(
      rho:f64,
      x:&Array1<f64>,
      b:&Array1<f64>)->f64;
    
    fn g(
      rho:f64,
      x:&Array1<f64>, 
      b:&Array1<f64>)->f64;

    fn which(&self)->String;

}

#[derive(Clone)]
pub struct CS;
#[derive(Clone)]
pub struct Kontogeorgis;

impl RdfModel for Kontogeorgis {


  fn model()->Self where Self: Sized {
      Self
  }
  fn g(
    rho:f64,
    x:&ndarray::Array1<f64>,
    b:&ndarray::Array1<f64>)->f64 {
    1.0 / (1.0 - 1.9 * Self::eta(rho, x, b))
  }

  fn dlngdeta(
    rho:f64,
    x:&Array1<f64>,
    b:&Array1<f64>)->f64 {
    let gmix = Self::g(rho, x, b);
    1.9 * gmix 
  }

  fn which(&self)->String{
    "Kontogeorgis RDF (sCPA)".to_string()
  }
}

impl RdfModel for CS {

  fn model()->Self where Self: Sized {
      Self
  }

  fn g(
    rho:f64,
    x:&Array1<f64>,
    b:&Array1<f64>)->f64 {
    let eta = Self::eta(rho, x, b);
    (1. - 0.5*eta) / (1.0 - eta).powi(3)
  }

  fn dlngdeta(
    rho:f64,
    x:&Array1<f64>,
    b:&Array1<f64>)->f64 {
    let eta = Self::eta(rho, x, b);
    (5. - 2.*eta)/(2. - eta)/(1. - eta)
  }

  fn which(&self)->String{
    "Carnahan-Starling RDF (CPA)".to_string()
  }
}