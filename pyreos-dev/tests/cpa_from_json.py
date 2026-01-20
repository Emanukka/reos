#%%

import numpy as np
from math import isclose


from reos.cpa import CPAPureRecord,CPABinaryRecord, CPAParameters
from reos.eos import EquationOfState
from reos.consts import Consts
# from reos.cubic import CubicBinaryRecord
R = Consts.ideal_gas_const()

# %%

p = CPAParameters.from_json(["water"],"../parameters/cpa/pure.json")


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