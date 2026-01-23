# CubicParameters Initializer

Initializer for **Cubic EoS parameters**.
The `opt` argument enables to choose the cubic model to be used.
This method must receive 3 arguments:

- `pure_records`: **[CubicPureRecord]**,
- `binary_records`: **[CubicBinaryRecord]**, optional
- `opt`: optional[str], cubic model option,
    - Default : **SRK model**
    - options are:
        - `"srk"`: Soave-Redlich-Kwong
        - `"pr76"`: Peng-Robinson 1976
        - `"pr78"`: Peng-Robinson 1978

