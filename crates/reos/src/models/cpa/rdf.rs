use ndarray::{Array1};

pub trait RDF {
    fn model()->Self where Self: Sized;
    fn eta(&self,rho:f64,x:&Array1<f64>,vb:&Array1<f64>)->f64{
        rho*vb.dot(x)/4.0
    }
    fn detadrho(&self,x:&Array1<f64>,vb:&Array1<f64>)->f64{
        vb.dot(x)/4.0
    }
    fn detadni(&self,rho:f64,vb:&Array1<f64>)->Array1<f64>{
        rho*vb/4.0
    }
    fn dlngdrho(&self,rho:f64,x:&Array1<f64>,vb:&Array1<f64>)->f64{
      let dlngdeta = self.dlngdeta(rho, x, vb);
      let detadrho = self.detadrho(x, vb);
      dlngdeta * detadrho
    }
    fn ndlngdni(&self,rho:f64,x:&Array1<f64>,vb:&Array1<f64>)->Array1<f64>{
      let dlngdeta = self.dlngdeta(rho, x, vb);
      let detadni = self.detadni(rho, vb);
      dlngdeta * detadni
    }

    fn dlngdeta(&self,rho:f64,x:&Array1<f64>,vb:&Array1<f64>)->f64;
    fn rdf(&self,rho:f64,x:&Array1<f64>, vb:&Array1<f64>)->f64;


    fn which(&self)->String;

}

#[derive(Clone)]
pub struct CarnahanStarlingRDF;
#[derive(Clone)]
pub struct ElliotRDF;

impl RDF for ElliotRDF {

  fn model()->Self where Self: Sized {
    ElliotRDF
  }
  fn rdf(&self,rho:f64,x:&ndarray::Array1<f64>,vb:&ndarray::Array1<f64>)->f64 {
    1.0 / (1.0 - 1.9 * self.eta(rho, x, vb))
  }

  fn dlngdeta(&self,rho:f64,x:&ndarray::Array1<f64>,vb:&ndarray::Array1<f64>)->f64 {
    let gmix = self.rdf(rho, x, vb);
    1.9 * gmix 
  }

  fn which(&self)->String{
    "Elliot RDF (sCPA)".to_string()
  }
}

impl RDF for CarnahanStarlingRDF {

  fn model()->Self where Self: Sized {
    CarnahanStarlingRDF
  }
  fn rdf(&self,rho:f64,x:&ndarray::Array1<f64>,vb:&ndarray::Array1<f64>)->f64 {
    let eta = self.eta(rho, x, vb);
    (1. - 0.5*eta) / (1.0 - eta).powi(3)
  }

  fn dlngdeta(&self,rho:f64,x:&ndarray::Array1<f64>,vb:&ndarray::Array1<f64>)->f64 {
    let eta = self.eta(rho, x, vb);
    (5. - 2.*eta)/(2. - eta)/(1. - eta)
  }

  fn which(&self)->String{
    "Carnahan-Starling RDF (CPA)".to_string()
  }
}