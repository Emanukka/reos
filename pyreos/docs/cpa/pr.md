# CPAPureRecord Initializer
### Cubic fields:

Can receive 2 types of pure cubic parameters as `kwargs`:
- `(a0, b, c1, tc)`
- `(tc, pc, w)`

    where: 
    - `a0` is the energy parameter (Pa·m^6/mol^2),
    - `b` is the co-volume parameter (m^3/mol),
    - `tc` is the critical temperature (K),
    - `pc` is the critical pressure (Pa),
    - `w` is the acentric factor (-).
    
    If both sets are provided, the first one will be used.
    If none is provided, an error will be raised.

### Associative fields:

- `(epsilon=0.0, kappa=0.0, na=0, nb=0, nc=0)`
    
    where: 
    - `epsilon` is the association energy (J/mol),
    - `kappa` is the association volume (-),
    - `na, nb, nc` are the number of association sites of type A (-), B (+) and C (±) respectively.

Example:
```py
from reos.cpa import CPAPureRecord

inert = CPAPureRecord.new("methane", 16.04, a0=1, b=2, c1=3, tc=4)
solvate = CPAPureRecord.new("carbon_dioxide", 44.01, a0=1, b=2, c1=3, tc=4, nb=1) #1 electron acceptor
associative = CPAPureRecord.new("water", 18.01528, a0=1, b=2, c1=3, tc=4, na=2, nb=2) #2 electron donors and 2 acceptors
```