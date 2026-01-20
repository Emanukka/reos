# CPABinaryRecord initializer


### Cubic fields

Can receive 2 types of binary cubic parameters as kwargs:

- Temperature Dependent: `aij,bij`, where: `kij = aij + bij * T`

- Temperature Independent: `kij`

    The default value is `None`.

### Associative fields

The following parameters will be attached to all interactions between sites `j` and `l`,
where `j` belongs to component `id1` and `l` belongs to component `id2`.
If none is provided, then the default value is the `cr1` combining rule: `epsilon` and `kappa`
from arithmetic and geometric mean extrapoleted from the pure components parameters, respectively.

- `epsilon`: cross association energy (J/mol) between `id1` and `id2`. Default is `None`.

- `kappa`: cross association volume (-) between `id1` and `id2`. Default is `None`.

- `rule`: combining rule for the sites inte.
        Options are: `"cr1"`, `"ecr"`. Default is `"cr1"`.

Example:

```py
bin = CPABinaryRecord.new("water", "carbon_dioxide", aij=-0.15508, bij=0.000877, kappa=0.1836)
print(bin)

BinaryRecord(id1="water", id2="carbon_dioxide", 
                model_record=CPABinaryRecord(
                    cubic=CubicBinaryRecord(aij=-0.15508, bij=0.000877),
                    assoc=AssociationBinaryRecord(epsilon=None, kappa=Some(0.1836), combining_rule=cr1)
                    )
                )
```