# REOS 



# The Project's Sctructre

```
├── crates
│   ├── reos
│   └── reos_py
├── python-examples
└── README.md
```

## reos

The rust core-base of project.

## reos_py 
Python Interface for the features inside reos-crate, made with pyo3 and maturin.

To install reos package in your python enviromnent, if you have cargo and maturin installed, just do the following command inside ```/reos_py```:

```
maturin develop --release
```
This command will install reos in your machine, and you can simply import it like any other python module:

```py
from reos.reos import EquationOfState,State,CPAParameters,CubicRecord,AssociationRecord

```
Otherwise, if you want to use reos inside a google colab notebook, you can use the ```reos_colab.ipynb```, where will be installed cargo and maturin, to compile the reos in a ```.wheel``` file, and finnaly import it .


## python_examples

This folder contains the usage of principal functionalities of the reos crate in Python through some examples:

- Phase Diagram of VLE for mixtures with binary mixtures with associative + inert mixtures and cross-associating mixtures
- Water solubility in CO2 vapor phase

