# CPABinaryRecord initializer

Parameters received as `kwargs`.

## General fields

id1: str,
    name of component 1
id2: str,
    name of component 2

### Cubic fields

- kij: dict[str, float], optional
    can have 2 fields: `a`, `b`, such that `kij = a + bT`

### Associative fields

The following parameters will be attached to all interactions between sites `j` and `l`,
where `j` belongs to component `id1` and `l` belongs to component `id2`.

- `assoc_rule`
    - `cr1`: str,
        use the following rule: `ϵⱼₗ = 0.5 (ϵⱼⱼ + ϵₗₗ)`, `κⱼₗ = √(κⱼⱼ κₗₗ)`
    - `ecr`: str
        uses the following rule: `Δⱼₗ = √(ΔⱼⱼΔₗₗ)`
    - `mcr1`: dict[str, float]
        fields: `kappa`
- `set`: dict[str, float],
    fields: `epsilon`, `kappa`