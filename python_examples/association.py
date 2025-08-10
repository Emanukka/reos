#%%
from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
from auxiliary_functions.parameters import *
from auxiliary_functions.association_functions import *


parameters = CPAParameters.from_records(
    cubic = [c_w,c_co2],
    assoc = [a_w,a_co2])

parameters.set_cubic_binary(
    j = 0,
    i = 1,
    kij_a =  0.000877,
    kij_b = -0.15508 )

parameters.set_assoc_binary(
    j = 0,
    i = 1,
    rule = "mcr1",
    beta = 0.1836)

eos = EquationOfState.cpa(parameters)



#%%

assoc=eos.get_association()

F=assoc.get_fmap()
T=assoc.get_tmat()
