
## A Thermodynamic State at `(T,P,x)` conditions

Units in SI.

Create a State object, trying to solve for density.


### Parameters
----------

- `eos` : **EquationOfState**
- `temperature` : float
- `pressure` : float
- `x` : numpy.ndarray[float]
- `phase` : Optional[str]
    - "stable" (default): choose the stable phase
    - "liquid": try liquid phase
    - "vapor": try vapor phase

### Returns
-------
State at (T,P,x) 

### Example
```py
name = "water"
molar_weight = 18.01528 

cubic = {"a0":0.12277, "b":0.0145e-3, "c1":0.6736, "tc":647.14}
assoc = {"epsilon":166.55e2, "kappa":0.0692, "na":2, "nb":2}

pr1 = CPAPureRecord.new(name, molar_weight, **cubic, **assoc)

p = CPAParameters.from_records([pr1])

eos = EquationOfState.scpa(p)

t = 298.15
press = 1e5
x = np.array([1.0])

s = State.tpx(eos, t, press, x) # default: most stable phase

d = s.density

assert(isclose(d, 55848.93240895648, rel_tol= 1e-9))

print(s)
# State(t = 298.150 K, p = 100000.000000 Pa, ρ = 55848.932409 mol/m³)
```