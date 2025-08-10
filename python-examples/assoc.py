#%%
from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from auxiliary_functions.parameters import *
# from auxiliary_functions.vle_functions import *
# from auxiliary_functions.data import *
from auxiliary_functions.wsat_data import *
from auxiliary_functions.water_sat import *
from auxiliary_functions.association_functions import *

#%%

c_co2 = CubicRecord(
    a0 = 0.35079,
    b  = 0.0272e-3,
    c1 = 0.7602,
    tc = 304.12)

c_w = CubicRecord(
    a0 = 0.12277,
    b  = 0.0145e-3,
    c1 = 0.6736,
    tc = 647.14)

a_co2 = AssociationRecord.solvate(
    b  = 0.0272e-3,
    na = 0,
    nb = 1,
    nc = 0)

a_w = AssociationRecord.associative(
    b    = 0.0145e-3,
    eps  = 166.55e2,
    beta = 0.0692,
    na   = 2,
    nb   = 2,
    nc   = 0)

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

T=298.15
P=1e5
x=np.array([0.5,0.5])

state= State.tpx(eos,T,P,x,density_initialization="vapor")

ln_phi=state.ln_phi()


#%%

assoc=eos.get_association()

F=assoc.get_fmap()
T=assoc.get_tmat()
