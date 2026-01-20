use crate::{impl_pure_record,impl_binary_record,impl_parameters,impl_eos};
use crate::parameters::records::*;
use reos::models::cubic::Cubic;
use reos::models::cubic::models::CubicModels;
use reos::models::cubic::{ 
    parameters::{CubicParameters,CubicBinaryRecord, CubicPureRecord},
};
use pyo3::pymethods;

use crate::{contribution::PyContribution, eos::PyEquationOfState};
use reos::state::eos::EquationOfState;
use reos::parameters::{Parameters,PureRecord,BinaryRecord};



impl_pure_record!(PyCubicPureRecord,CubicPureRecord,"CubicPureRecord", "../../docs/cubic/pr.md");

impl_binary_record!(PyCubicBinaryRecord,CubicBinaryRecord,"CubicBinaryRecord", "../../docs/cubic/br.md");

// type PR78Parameters = CubicParameters<PR78_MARKER>;
// type PR76Parameters = CubicParameters<PR76>;

impl_parameters!(PyCubicParameters,CubicParameters,PyCubicPureRecord,PyCubicBinaryRecord,"CubicParameters", CubicModels,"../../docs/cubic/parameters.md",);

impl_eos!(Cubic, Cubic, PyCubicParameters, cubic, "../../docs/cubic/eos.md");
