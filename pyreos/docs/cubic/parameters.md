# CubicParameters Initializer

Initializer for **Cubic EoS parameters**.

- `pure_records`: **[CubicPureRecord]**,
- `binary_records`: **[CubicBinaryRecord]**, optional
- `opt`: kwargs,
    - fields:
        - `cubic_model`
            - `srk`
            - `pr76`
            - `pr78`
        - `alpha_model`
            - `soave`: uses α(Tᵣ)= [1 + c₁(1 - √Tᵣ)]²
            - `twu91`: uses α(Tᵣ)= Tᵣᴺ⁽ᴹ⁻¹⁾exp[ L(1 - Tᵣᴹᴺ)]
            - (default `soave`)