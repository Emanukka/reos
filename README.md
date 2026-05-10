# Welcome to ℜ*eos*

[![docs](https://img.shields.io/badge/docs-github--pages-blue?logo=github)](https://emanukka.github.io/reos/)

ℜ*eos* is a thermodynamic library written in **Rust** with a **Python interface**. It provides tools for calculating thermodynamic properties and phase equilibria.

```python
import numpy as np
from reos.cpa import CPAParameters, CPAPureRecord
from reos.eos import EquationOfState
from reos.state import State

parameters = CPAParameters.from_json(["water", "carbon dioxide"], 
                             rdf_model="kg",
                             cubic_model="srk",
                             ppath="../../../parameters/cpa/tsivintzelis2011.json", 
                             bpath= "../../../parameters/cpa/tsivintzelis2011_binary.json")

eos = EquationOfState.cpa(parameters)

# 2. Initialize the state
T = 323.15; P = 1e5; z = np.array([0.75, 0.25])
state = State.tpx(eos, T, P, z)

# 3. Calculate properties
unbonded_sites_fraction = eos.get_assoc_calcs(state.temperature, state.density, state.composition)["X"]

```

## Models

| Model | Description |
|:-----:|:-----------:|
| CPA   | CPA, SCPA |
| cubic | SRK, PR76, PR78, Twu-91, volume-translation|

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

- `parameters`: Directory with pure and binary parameters for the current models implemented in Reos

- `crates`: Contains `reos`, which is the Rust package that implement all the core functionalites

- `pyreos`: Rust package that create the Python interface from `reos`

- `pyreos-dev`: Directory used to test functionalites of `reos` in Python

- `pyreos-examples`: Python examples of how to use `reos`

### Python package

- `reos.eos`: Enables the creation and manipulation of equations of state with different models.

- `reos.state`: Provides tools for working with thermodynamic states, including property calculations and phase equilibria.

- `reos.{model_name}`: Each current model has its own submodule with its name, which contains the **pure model record**, the **binary model record** and the **parameters objects** .
