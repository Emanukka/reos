

```python
# 1. Initialize the model
parameters = CPAParameters.from_json(
    ["water"], 
    ppath="../../../parameters/cpa/kontogeorgis2006.json", 
    rdf_model="kg", 
    cubic_model="srk")

eos = EquationOfState.cpa(parameters)

# 2. Initialize the state
T = 298.15; P = 1e5
state = State.pure_tp(eos, T, P)

# 3. Calculate properties
sites = parameters.get_sites()
sites_fraction = eos.get_assoc_calcs(T, state.density, state.composition)["X"]
```