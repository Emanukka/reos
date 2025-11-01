# REOS 

```python
import numpy as np
from reos.reos import EquationOfState,State,CPAParameters,CubicRecord,AssociationRecord



#1) Cubic pure parameters
c_water = CubicRecord(a0 = 0.12277, b = 0.0145e-3, c1 = 0.6736, tc = 647.14)
c_co2 = CubicRecord(a0 = 0.35079, b = 0.0272e-3, c1 = 0.7602, tc = 304.12)


#2) Associative pure parameters (Water 4C and CO2 1ea, data from  Tsivintzelis,2011)  
a_water  = AssociationRecord.associative(b = 0.0145e-3, eps = 166.55e2, beta = 0.0692, na = 2, nb = 2, nc = 0)
a_co2  = AssociationRecord.solvate(b = 0.0272e-3, na = 0, nb = 1, nc = 0)

#3) Build a CPAParameters object 
parameters = CPAParameters.from_records(
    cubic = [c_water,c_co2],
    assoc = [a_water,a_co2])

# Set the combination's rule 
parameters.set_cubic_binary(j = 0, i = 1, kij_a =  0.000877, kij_b = -0.15508 ) #kij = aT + b
parameters.set_assoc_binary(j = 0, i = 1, rule = "mcr1", beta = 0.1836) # modified CR1

#4) Create a EoS object w\ CPA contribution
eos = EquationOfState.cpa(parameters)

#5) Instantiate a State object from the EoS
T=298.15
P=500e5
x=np.array([0.5,0.5])
state= State.tpx(eos,T,P,x,density_initialization="vapor")

# Compute thermodynamic properties
phi=np.exp(state.ln_phi())
print(f"phi = {phi}")

# Get association variables
assoc = eos.get_association()
sites_map =assoc.get_sites_map()

X=state.non_bonded_sites()

print(f"Sites = {sites_map}") 

print(f"Non-Bonded Sites fraction = {X}")
```
   
```
phi = [2.14385748e-04 5.65853284e-01]
Sites = ((0, 0), (1, 0), (1, 1)) <span style="color:green">#(type, owner)</span>
Non-Bonded Sites fraction = [0.07333045 0.22590468 0.69485154]
```

```


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
Otherwise, if you want to use reos inside a google colab notebook, you can use the ```reos_colab.ipynb``` ( located in ```/python_examples```), where will be installed cargo and maturin, to compile the reos in a ```.wheel``` file, and finnaly import it .


## python_examples

This folder contains the usage of principal functionalities of the reos crate in Python through some examples:

- Phase Diagram of VLE for mixtures with binary mixtures with associative + inert mixtures and cross-associating mixtures
- Water solubility in CO2 vapor phase

