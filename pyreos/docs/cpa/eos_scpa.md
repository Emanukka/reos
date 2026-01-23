# SCPA Initializer

Initializer for **SCPA EoS**.
**SCPA** combines the cubic EoS and the association term, which has a simplification at the Radial Distribution Function `g` (see [Kontogeorgis et al, 1999](https://doi.org/10.1016/S0378-3812(99)00060-6)).

The cubic model is selected when initializing the **CPAParameters** (see [CPAParameters](./parameters.md)). If none is provided, the default is the **SRK** model.

$$ g(\rho) = \frac{1}{1-0.475b \rho}$$

Example:

```py
name = "water"
molar_weight = 18.01528 

cubic = {"a0":0.12277, "b":0.0145e-3, "c1":0.6736, "tc":647.14}
assoc = {"epsilon":166.55e2, "kappa":0.0692, "na":2, "nb":2}

pr1 = CPAPureRecord.new(name, molar_weight, **cubic, **assoc)

p = CPAParameters.from_records([pr1]) # default to SRK model

eos = EquationOfState.scpa(p)
t = 298.15
d = 1000.0
x = np.array([1.0])

P = eos.pressure(t, d, x)
Pres = P - eos.ideal_gas_pressure(t, d)
Pres_reduced = Pres / R / t
chem_pot_res_reduced = eos.chem_pot(t, d, x)[0] 
helmholtz_res_reduced = eos.helmholtz(t, d, x) / R / t
entropy_res_reduced = eos.entropy(t, d, x) / R

assert(isclose(Pres_reduced, -57.5159551979349 + -945.9409464127781, rel_tol=1e-9))
assert(isclose(chem_pot_res_reduced, -0.115660251059 + -2.54386196979185, rel_tol=1e-9))
assert(isclose(helmholtz_res_reduced, -0.058144295861 + -1.597921023379, rel_tol=1e-9))
assert(isclose(entropy_res_reduced, -0.041951593945 + -4.713659269705, rel_tol=1e-9))

```
