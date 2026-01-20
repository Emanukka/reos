#%%

import numpy as np
from math import isclose


from reos.cpa import CPAPureRecord,CPABinaryRecord, CPAParameters
from reos.eos import EquationOfState
from reos.consts import Consts
# from reos.cubic import CubicBinaryRecord
R = Consts.ideal_gas_const()
#%%instance


inert = CPAPureRecord.new("methane", 16.04, a0=1, b=2, c1=3, tc=4)
solvate = CPAPureRecord.new("carbon_dioxide", 44.01, a0=1, b=2, c1=3, tc=4, nb=1) #1 electron acceptor
associative = CPAPureRecord.new("water", 18.01528, a0=1, b=2, c1=3, tc=4, epsilon=1.0,kappa=2.0, na=2, nb=2) #2 electron donors and 2 acceptors



# %%

data = {"aij":-0.15508, "bij":0.000877, "kappa":0.1836,"rule":"cr1"}
bin = CPABinaryRecord.new("water", "carbon_dioxide", **data)


print(bin)

#%% scpa eos
name = "water"
molar_weight = 18.01528  # g/mol

cubic = {"a0":0.12277, "b":0.0145e-3, "c1":0.6736, "tc":647.14}
assoc = {"epsilon":166.55e2, "kappa":0.0692, "na":2, "nb":2}

pr1 = CPAPureRecord.new(name, molar_weight, **cubic, **assoc)

p = CPAParameters.from_records([pr1])

#%%

name = "water"
molar_weight = 18.01528 

cubic = {"a0":0.12277, "b":0.0145e-3, "c1":0.6736, "tc":647.14}
assoc = {"epsilon":166.55e2, "kappa":0.0692, "na":2, "nb":2}

pr1 = CPAPureRecord.new(name, molar_weight, **cubic, **assoc)

p = CPAParameters.from_records([pr1])

eos = EquationOfState.scpa(p)
t = 298.15
d = 1000.0
x = np.array([1.0])

P = eos.pressure(t, d, x)
Pres = P - eos.ideal_gas_pressure(t, d)
Pres_reduced = Pres / R / t
chem_pot_res_reduced = eos.chem_pot(t, d, x)[0] / R / t 
helmholtz_res_reduced = eos.helmholtz(t, d, x) / R / t
entropy_res_reduced = eos.entropy(t, d, x) / R

assert(isclose(Pres_reduced, -57.5159551979349 + -945.9409464127781, rel_tol=1e-9))
assert(isclose(chem_pot_res_reduced, -0.115660251059 + -2.54386196979185, rel_tol=1e-9))
assert(isclose(helmholtz_res_reduced, -0.058144295861 + -1.597921023379, rel_tol=1e-9))
assert(isclose(entropy_res_reduced, -0.041951593945 + -4.713659269705, rel_tol=1e-9))

#%%