

pub trait PureRecord{

}

pub trait Parameters<R: PureRecord>{

  fn from_records(records:Vec<R>)->Self;
  
}