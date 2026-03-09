# CPAParameters Initializer

Initializer for **CPA EoS** parameters.

This method must receive 3 arguments:

- `pure_records`: **[CPAPureRecord]**,

- `binary_records`: **[CPABinaryRecord]**, optional

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