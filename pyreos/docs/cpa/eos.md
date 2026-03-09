# CPA Initializer

Initializer for **CPA EoS**.
**CPA** combines the cubic EoS and the association term(see [Kontogeorgis et al, 1999](https://doi.org/10.1016/S0378-3812(99)00060-6)).

The specific cubic model and others model's specifications are selected when initializing the parameters.

<!-- Example:

```py
name = "water"
molar_weight = 18.01528 

cubic = {"a0":0.12277, "b":0.0145e-3, "c1":0.6736, "tc":647.14}
assoc = {"epsilon":166.55e2, "kappa":0.0692, "na":2, "nb":2}

pr1 = CPAPureRecord.new(name, molar_weight, **cubic, **assoc)

p = CPAParameters.from_records([pr1], )

eos = EquationOfState.cpa(p) -->


<!-- ``` -->
