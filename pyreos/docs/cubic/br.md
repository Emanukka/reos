# CubicBinaryRecord initializer

Cubic Binary Record Initializer

Parameters

---------------
kwargs:
id1: str,
    name of component 1
id2: str,
    name of component 2
kij: dict[str, float], optional
    can have 2 fields: `a`, `b`, such that `kij = a + bT`

Example:

```py

set1 = CubicBinaryRecord.new(id1="water", id2="carbon dioxide", aij=0.1, bij=0.01)
set2 = CubicBinaryRecord.new(id1="water", id2="carbon dioxide", kij=0.2)

```
