# Cubic Pure Record Initializer

Parameters received as `kwargs`.

---------------

- `name`: str
    component's name
- `molar_weight`: float
    component's molar weight
- `tc`: float
    critial temperature
- `pc`: float, optional
    critical pressure
- `a`: float, optional
    attractive parameter
- `b`: float, optional
    volumetric parameter
- `w`: float, optional
    acentric factor soave alpha function
- `c1`: float, optional
    soave alpha function parameter
- `l,m,n`: float, optional
    twu91 alpha function parameters

```py
from reos.cubic import CubicPureRecord

set1 = CubicPureRecord.new(name="water", molar_weight=18.01528, a0=0.0, b=0.0, c1=0.0, tc=0.0)
set2 = CubicPureRecord.new(name="water", molar_weight=18.01528, tc=0.0, pc=0.0,w=0.0)
```
