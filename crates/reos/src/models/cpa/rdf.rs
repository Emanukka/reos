
use enum_dispatch::enum_dispatch;
use ndarray::{Array1, Array2};

#[enum_dispatch]
pub trait RDFmodel {

    fn eta(&self, rho:f64, x:&Array1<f64>, b:&Array1<f64>) -> f64{
        rho*b.dot(x) / 4.0
    }

    fn dlngdeta(&self, rho:f64, x:&Array1<f64>, b:&Array1<f64>)->f64;
    
    fn g(&self, rho:f64, x:&Array1<f64>, b:&Array1<f64>)->f64;

	
	fn to_string(&self,) -> String;

}

#[derive(Clone)]
pub struct CS;
#[derive(Clone)]
pub struct KG;

impl RDFmodel for KG {
	
	fn g(&self, rho:f64, x:&ndarray::Array1<f64>, b:&ndarray::Array1<f64>)->f64 {
		1.0 / (1.0 - 1.9 * self.eta(rho, x, b))
  	}

  	fn dlngdeta(&self, rho:f64, x:&Array1<f64>, b:&Array1<f64>) -> f64 {
    	let gmix = self.g(rho, x, b);
    	1.9 * gmix 

  	}

  	fn to_string(&self,) -> String {

		"KG".into()

  	}

}

impl RDFmodel for CS {


  	fn g(&self, rho:f64, x:&Array1<f64>, b:&Array1<f64>) -> f64 {
  		let eta = self.eta(rho, x, b);
  	  	(1. - 0.5*eta) / (1.0 - eta).powi(3)
	
  	}

  	fn dlngdeta(&self, rho:f64, x:&Array1<f64>, b:&Array1<f64>) -> f64 {
  	  	let eta = self.eta(rho, x, b);
  	  	(5. - 2.*eta)/(2. - eta)/(1. - eta)
  	}

	fn to_string(&self,) -> String {
		"CS".into()
	}

}
#[enum_dispatch(RDFmodel)]
#[derive(Clone)]
pub enum RDF {
	CS,
	KG
}



#[derive(Clone)]
pub struct RDFcpa{
  pub b:Array1<f64>,
  pub bij:Array2<f64>,
  pub model: RDF

}



impl RDFcpa {

    pub fn detadrho(&self,x:&Array1<f64>)->f64{

        self.b.dot(x)/4.0

    }
    pub fn detadni(&self,rho:f64)->Array1<f64>{

        rho * &self.b/4.0

    }
    pub fn dlngdrho(&self,rho:f64,x:&Array1<f64>)->f64{

      let dlngdeta = self.model.dlngdeta(rho, x,&self.b);
      let detadrho = self.detadrho(x);

      dlngdeta * detadrho
    }
    pub fn ndlngdni(&self,rho:f64,x:&Array1<f64>)->Array1<f64>{

      let dlngdeta = self.model.dlngdeta(rho, x, &self.b);
      let detadni = self.detadni(rho);

      dlngdeta * detadni

    }
    pub fn g(&self,rho:f64,x:&Array1<f64>)->f64{

      self.model.g(rho, x, &self.b)

    }


}

#[derive(serde::Deserialize, serde::Serialize, Clone, Debug)]
#[serde(rename_all="lowercase")]
pub enum RDFmodelOption{
	CS,
	KG
}

impl From<RDFmodelOption> for RDF {
	fn from(value: RDFmodelOption) -> Self {
		match value {
			RDFmodelOption::KG => RDF::KG(KG),
			RDFmodelOption::CS => RDF::CS(CS)
		}
	}
}

impl Default for RDFmodelOption {

	fn default() -> Self {
		Self::KG
	}
} 