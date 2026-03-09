# CPAPureRecord Initializer

Parameters received as `kwargs`.

---------------

## General fields

- name: str
    component's name
- molar_weight: float
    component's molar weight

## Cubic fields

- tc: float
    critial temperature
- pc: float, optional
    critical pressure
- a: float, optional
    attractive parameter
- b: float, optional
    volumetric parameter
- w: float, optional
    acentric factor soave alpha function,
- c1: float, optional
    soave alpha function parameter,
- l,m,n: float, optional
    twu91 alpha function parameters

## Associative fields

- epsilon: float, default (0.0)
    associative energy parameter
- kappa: float, default (0.0)
    associative binding volume
- na: int, default(0)
    multiplicity of A sites
- nb: int, default(0)
    multiplicity of B sites
- nc: int, default(0)
    multiplicity of C sites

Example:

```py
from reos.cpa import CPAPureRecord

inert = CPAPureRecord.new(name="methane", molar_weight=16.04, a0=1, b=2, c1=3, tc=4)
solvate = CPAPureRecord.new(name="carbon_dioxide", molar_weight=44.01, a0=1, b=2, c1=3, tc=4, nb=1) #1 electron acceptor
associative = CPAPureRecord.new(name="water", molar_weight=18.01528, a0=1, b=2, c1=3, tc=4, na=2, nb=2) #2 electron donors and 2 acceptors
```
