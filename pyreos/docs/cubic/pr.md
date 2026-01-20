# CubicPureRecord Initializer



This method can receive 2 types of pure cubic parameters as `kwargs`:
- `a0, b, c1, tc`
- `tc, pc, w,`

    where: 
    
    - `a0` is the energy parameter (Pa·m^6/mol^2),
    - `b` is the co-volume parameter (m^3/mol),
    - `tc` is the critical temperature (K),
    - `pc` is the critical pressure (Pa),
    - `w` is the acentric factor (-).

If both sets are provided, the first one will be used.
If none is provided, an error will be raised.


Example:

```py
from reos.cubic import CubicPureRecord
name = "water"
molar_weight = 18.01528  # g/mol
set1 = CubicPureRecord.new(name, molar_weight, a0=0.0, b=0.0, c1=0.0, tc=0.0)
set2 = CubicPureRecord.new(name, molar_weight, tc=0.0, pc=0.0,w=0.0)
```
