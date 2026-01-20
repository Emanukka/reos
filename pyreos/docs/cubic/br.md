# CubicBinaryRecord initializer

This method can receive 2 types of binary cubic parameters as `kwargs`:

- Temperature Dependent: `aij,bij`, where: `kij = aij + bij * T`

- Temperature Independent: `kij`

Example:

```py
id1 = "water"
id2 = "co2"
set1 = CubicBinaryRecord.new(id1, id2, aij=0.1, bij=0.01)
set2 = CubicBinaryRecord.new(id1, id2, kij=0.2)

```