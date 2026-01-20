# Cubic Initializer


Initializer for **Cubic EoS**. The implementation of Cubic derive from the generic implementation of a Cubic EoS.
(See [Martin, 1979](https://doi.org/10.1021/i160070a001))

The specific cubic model is selected when initializing the parameters (see [CubicParameters](./parameters.md)).

Example:
```py
name = "water"
molar_weight = 18.01528  # g/mol

pr1 = CubicPureRecord.new(name, molar_weight, a0=0.12277, b=0.0145e-3, c1=0.6736, tc=647.14)
p = CubicParameters.from_records([pr1]) # default to SRK model

srk = EquationOfState.cubic(p)
t = 298.15
d = 1000.0
x = np.array([1.0])

P = srk.pressure(t, d, x)
Pres = P - srk.ideal_gas_pressure(t, d)
Pres_reduced = Pres / R / t
chem_pot_res_reduced = srk.lnphi(t, d, x)[0] + np.log(srk.compressibility(t, d, x))
helmholtz_res_reduced = srk.helmholtz(t, d, x) / R / t
entropy_res_reduced = srk.entropy(t, d, x) / R

assert(isclose(Pres_reduced, -57.5159551979349, rel_tol=1e-10))
assert(isclose(chem_pot_res_reduced, -0.115660251059, rel_tol=1e-10))
assert(isclose(helmholtz_res_reduced, -0.058144295861, rel_tol=1e-10))
assert(isclose(entropy_res_reduced, -0.041951593945, rel_tol=1e-10))


```