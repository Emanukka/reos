# Welcome to ℜ*eos*

ℜ*eos* is a thermodynamic library written in **Rust** with a **Python interface**. It provides tools for calculating thermodynamic properties and phase equilibria.


# Installation

## Python

```bash
pip install reos
```

## From source


# Project layout

## Python package

- `reos.eos`: Enables the creation and manipulation of equations of state with different models.
- `reos.state`: Provides tools for working with thermodynamic states, including property calculations and phase equilibria.
- `reos.{model_name}`: Each current model has its own submodule with its name, which contains the **pure model record**, the **binary model record** and the **parameters objects** 
