# REOS 

```py
#Array with [A,B,C] for each component
antoine=np.array([metoh_antoine,octanol_antoine])

p=CPAParameters.from_records(
    cubic=[c_methanol_3b,c_octanol_3b],
    assoc=[a_methanol_3b,a_octanol_3b])

p.set_cubic_binary(0,1,0.0,-0.025)
p.set_assoc_binary(0,1,"ecr")
eos_3b=EquationOfState.cpa(p)
P=101.32e3

VAR,LIQUID,VAPOR=linspace_bubble_t(eos_3b,P,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   title="Methanol(1) and Octanol(2) 3B+ECR ",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
   text=f"{P/1e3}kPa",
   factor=1.0,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=300,y_sup=480,
   exp_data=exp_data)

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]

#Get the non-bonded sites fraction from each state
XL=np.array([LIQUID[i].non_bonded_sites() for i in range(100)])
XV=np.array([VAPOR[i].non_bonded_sites() for i in range(100)])

sites=["-MeOH","+MeOH","-OcOH","+OcOH"]

plt.figure(figsize=(5, 2.5))

for (i,s) in enumerate(sites):

    plt.plot(xL,XL[:,i],label=s)

plt.title("Methanol(1) and Octanol(2) 3B+ECR Liquid Phase")
plt.ylabel("Fraction of Non-Bonded Sites")

plt.xlabel(r"$x_1$") 
plt.legend()

```
![alt text](<python_examples/phase_diagram_plot/Methanol(1) and Octanol(2) 3B+ECR .png>)
![alt text](python_examples/xassoc_plot/Methanol(1)_and_Octanol(2)_3B+ECR_VaporPhase.png)

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

