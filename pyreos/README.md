# Welcome to ℜ*eos*

[![docs](https://img.shields.io/badge/docs-github--pages-blue?logo=github)](https://emanukka.github.io/reos/)

ℜ*eos* is a thermodynamic library written in **Rust** with a **Python interface**. It provides tools for calculating thermodynamic properties and phase equilibria.

```python
import numpy as np
from reos.cpa import CPAParameters, CPAPureRecord
from reos.eos import EquationOfState
from reos.state import State

parameters = CPAParameters.from_records([ 
    CPAPureRecord.new(
        name = "water",
        molar_weight = 18.01528,
        a0 = 0.12277,
        b = 0.014515e-3,
        c1 = 0.67359,
        tc = 647.29,
        epsilon = 166.55e2,
        kappa = 0.0692,
        na = 2,
        nb = 2
    )
])

# or  CPAParameters.from_json(["water"], "./parameters/cpa/kontogeorgis2006.json")

eos = EquationOfState.scpa(parameters)

t = 298.15
p = 1e5
x = np.array([1.0])

s = State.tpx(eos, t, p, x) 

print(s)

X = eos.unbonded_sites_fraction(t, s.density, x)

print(f"Unbonded sites fraction = {X}")

```

```bash
State(t = 298.150 K, p = 100000.000000 Pa, ρ = 55784.919890 mol/m³)
Unbonded sites fraction = [0.07825138 0.07825138]
```

## Models

| Model | Description |
|:-----:|:-----------:|
| CPA   | Cubic Plus Association (srk, pr76, pr78)|
| cubic | Soave-Redlich-Kwong, Peng-Robinson 1976, Peng-Robinson 1978|

Each model implement its analytical expressions of derived properties of Helmholtz potential.

## Installation

### Python

#### From pip

```bash
pip install reos
```

#### From source

You must have rust compiler and [maturin/PyO3](https://pyo3.rs/v0.27.2/) installed.

To install the python package in your machine (such `pip install -e .`) and use the optimized version of the package:

```bash
cd pyreos
maturin develop --release 
```

To build the the python wheel, just use

```bash
maturin build --release 
```

## Project layout

### Repository

```bash
├── crates
├── parameters
├── pyreos
├── pyreos-dev
└── pyreos-examples
```

- `crates`: Contains `reos`, which is the Rust package that implement all the core functionalites
- `pyreos`: Rust package that create the Python interface from `reos`
- `pyreos-dev`: Directory used to test functionalites of `reos` in Python
- `pyreos-examples`: Python examples of how to use `reos`

### Python package

- `reos.eos`: Enables the creation and manipulation of equations of state with different models.
- `reos.state`: Provides tools for working with thermodynamic states, including property calculations and phase equilibria.
- `reos.{model_name}`: Each current model has its own submodule with its name, which contains the **pure model record**, the **binary model record** and the **parameters objects** .
