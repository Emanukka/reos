#%%
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from os import listdir
from reos.eos import EquationOfState
from reos.state import State
from reos.cpa import CPAParameters

# import water_saturation
# from water_saturation.methods import *
from wsat.methods import linspace_wsat
from si_units import KELVIN, BAR, PASCAL, CELSIUS
# import si_units

#%% 

# 1. Initialize the model
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

#%% 1. Getting the water saturation exp. data
sheet_name = str(T) 
excel = pd.ExcelFile("data/carbon_dioxide@1_0.xlsx")
exp_data = pd.read_excel(excel, sheet_name = sheet_name,)

y_dry_gas = np.array([0., 1.])

#%% 2. Getting the CPA parameters for pure water
pw = CPAParameters.from_json(["water"],
                             rdf_model="kg",
                             cubic_model="srk",
                             ppath="../../../parameters/cpa/tsivintzelis2011.json")
eosw = EquationOfState.cpa(pw)

#%% 3. Calculate saturation
N = 100
vpressure, yw, vdeveloped, vincipient = linspace_wsat(
    eos, 
    eosw, 
    T, 
    y_dry_gas, 
    p0=1e5, 
    pf=600e5, 
    N=N)

#%% 4. Obtain unbonded sites fraction for the water poor phase
Xmat = np.zeros((N, 3))

for i,state in enumerate(vdeveloped):

    X = eos.get_assoc_calcs(state.temperature, state.density, state.composition)["X"]
    Xmat[i, :] = X
#%% 5. Plot
plt.rcParams.update({
    "text.usetex": True,
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 30,
    "font.size": 46,
    "legend.fontsize": 25,

    "xtick.labelsize": 30,
    "ytick.labelsize": 30,

})

fig1, axes1 = plt.subplots(1, 2, figsize=(14,5), constrained_layout=True)

axes1[0].plot(vpressure/1e5, yw*1e6, "black")

axes1[0].scatter(exp_data["p"], exp_data["y"]*1e6, marker="o", label = f"{T} K", s = 300,facecolors="None", edgecolors="black")

axes1[0].set_xlabel(r"$ P / \mathrm{bar}$")
axes1[0].set_ylabel(r"$ \mathrm{Water \ content / ppm}$")
axes1[0].set_ylim(0, 15_000)
axes1[0].set_xlim(0, 600)
axes1[0].legend()

Xlabels = [r"$\rm{H_2O \ -}$",r"$\rm{H_2O \ +}$", r"$\rm{CO_2 \ +}$", ]
colors = ["#3b5b92", "#5fe6d6", "#c44536"]

for k in range(3):

    axes1[1].plot(vpressure/1e5, Xmat[:, k], label = Xlabels[k], color = colors[k])

axes1[1].set_xlabel(r"$ P / \mathrm{bar}$")
axes1[1].set_ylabel(r"$ X$")
axes1[1].set_ylim(0, 1)
axes1[1].set_xlim(0, 600)
axes1[1].legend()
