#%%
import numpy as np
from math import isclose
from reos.cpa import CPAParameters
from reos.eos import EquationOfState
from reos.state import State
# from reos.consts import Consts

#%%

p = CPAParameters.from_json(["water"],"../parameters/cpa/pure.json")
eos = EquationOfState.scpa(p)
t = 298.15
d = 1000.0
x = np.array([1.0])
X = eos.unbonded_sites_fraction(t, d, x)
X_rust = [0.530287110928 , 0.530287110928]

assert (all([isclose(X[i], X_rust[i], rel_tol=1e-9) for i in range(len(X))]))

#%%
p = CPAParameters.from_json(names = ["water","carbon dioxide"], 
                            ppath = "../parameters/cpa/pure.json",
                            bpath = "../parameters/cpa/binary.json")

eos = EquationOfState.scpa(p)
t = 298.15
d = 55_000
x = np.array([0.5, 0.5])
X = eos.unbonded_sites_fraction(t, d, x)
X_rust = [0.038741216735623495, 0.20484626829381528, 0.6677898967806235]

assert (all([isclose(X[i], X_rust[i], rel_tol=1e-9) for i in range(len(X))]))

#%%
import numpy as np
from reos.cpa import CPAParameters
from reos.eos import EquationOfState
from reos.state import State

parameters = CPAParameters.from_json(["water"], "../../parameters/cpa/kontogeorgis2006.json")
eos = EquationOfState.scpa(parameters)

t = 298.15
p = 1e5
x = np.array([1.0])

s = State.tpx(eos, t, p, x) 

print(s)

X = eos.unbonded_sites_fraction(t, s.density, x)

print(f"Unbonded sites fraction = {X}")
