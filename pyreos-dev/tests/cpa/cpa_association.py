#%%
import numpy as np
from math import isclose
from reos.cpa import CPAParameters
from reos.eos import EquationOfState
from reos.state import State
# from reos.consts import Consts


# %%

#%%

opt = {
    "rdf_model": "kg",
    "cubic_model": "srk"
}

parameters = CPAParameters.from_json(["water"], "../../parameters/cpa/pure.json", **opt)

eos = EquationOfState.cpa(parameters)


#%%
t = 298.15
d = 58_000.
x = np.array([1.0])


assoc_calcs = eos.get_assoc_calcs(t, d, x)

X, Δ, K = assoc_calcs["X"], assoc_calcs["Delta"], assoc_calcs["K"]

X_rust = [0.07588168172556269, 0.07588168172556281]
Δ_rust = 6.02214076e23 *  np.array([[0., 2.2974432513007877e-27],[2.2974432513007877e-27, 0.]])

assert (np.allclose(X, X_rust, rtol = 1e-15))
assert(np.allclose(Δ, Δ_rust, rtol = 1e-15))
assert(np.allclose(K, d * Δ_rust , rtol = 1e-15))


#%%
parameters = CPAParameters.from_json(names = ["water","carbon dioxide"], 
                            ppath = "../../parameters/cpa/pure.json",
                            bpath = "../../parameters/cpa/binary.json",
                            **opt)

eos = EquationOfState.cpa(parameters)
t = 298.15
d = 1_000
x = np.array([0.7, 0.3])

assoc_calcs = eos.get_assoc_calcs(t, d, x)

X, Δ, K, m = assoc_calcs["X"], assoc_calcs["Delta"], assoc_calcs["K"], assoc_calcs["m"]


# X = eos.unbonded_sites_fraction(t, d, x)
X_rust = [0.5786343099065254, 0.5957666698192122, 0.9200489870741287]

Δ_rust = 6.02214076e23 * np.array([[0.0, 1.390825988133717e-27, 1.7812705976946612e-28],
                  [1.390825988133717e-27, 0.0, 0.0],
                  [1.7812705976946612e-28, 0.0, 0.0]])

mjml = m * m[:, None]

assert (np.allclose(X, X_rust, rtol = 1e-15))
assert(np.allclose(Δ, Δ_rust, rtol = 1e-15))
assert(np.allclose(K, d * mjml * Δ_rust , rtol = 1e-15))


#%%

import numpy as np
from reos.cpa import CPAParameters, CPAPureRecord
from reos.eos import EquationOfState
from reos.state import State

parameters = CPAParameters.from_records([
    CPAPureRecord.new(
        name = "water",
        molar_weight = 18.01528,
        a0 = 0.12277,
        b = 0.014515e-3,
        c1 = 0.67359,
        tc = 647.29,
        epsilon = 166.55e2,
        kappa = 0.0692,
        na = 2,
        nb = 2)],
    rdf_model = "kg",
    cubic_model = "srk"
    )

eos = EquationOfState.cpa(parameters)

t = 298.15
p = 1e5
x = np.array([1.0])

s = State.tpx(eos, t, p, x) 

print(s)

assoc_calcs = eos.get_assoc_calcs(t, d, x)

X, Δ = assoc_calcs["X"], assoc_calcs["Delta"]

# X = eos.unbonded_sites_fraction(t, s.density, x)

print(f"Unbonded sites fraction = {X}")
# print(f"Association strength = {Δ}")
