# REOS

ℜ*eos* is a thermodynamic library written in **Rust** with a **Python interface**. It provides tools for calculating thermodynamic properties and phase equilibria.

```python
import numpy as np
from reos.cpa import CPAParameters
from reos.eos import EquationOfState
from reos.state import State

parameters = CPAParameters.from_json(["water"], "../../parameters/cpa/kontogeorgis2006.json")
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

|name|description|
|:-:|:-:|
|`cpa`|Cubic Plus Association (srk, pr76 or pr78)|
|`cubic`|Soave-Redlich-Kwong, Peng-Robinson 1976 or Peng-Robinson 1978|

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
├── docs
├── mkdocs.yml
├── parameters
├── pyreos
├── pyreos-dev
├── pyreos-examples
├── README.md
└── site
```

- `crates`: Contains `reos`, which is the Rust package that implement all the core functionalites
- `pyreos`: Rust package that create the Python interface from `reos`
- `pyreos-dev`: Directory used to test functionalites of `reos` in Python
- `pyreos-examples`: Python examples of how to use `reos`

### Python package

- `reos.eos`: Enables the creation and manipulation of equations of state with different models.
- `reos.state`: Provides tools for working with thermodynamic states, including property calculations and phase equilibria.
- `reos.{model_name}`: Each current model has its own submodule with its name, which contains the **pure model record**, the **binary model record** and the **parameters objects** .

## About
