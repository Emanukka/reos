#%%
from reos.scaling import Dehlouz
from reos.state import State
from reos.eos import EquationOfState as E
from reos.cubic import CubicParameters
from reos.consts import Consts
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt

R = Consts.ideal_gas_const()

param_scaling = [0.819374, 0.159219, 0.602420, 1.278829, 0.682719, 0.681497]


x = np.array([1.])
# %%

# t = 198.15
# p = 1e5
parameters = CubicParameters.from_json(["water"], "./parameters.json", opt = "pr78" )
eos = E.cubic(parameters)
state_critical = State.pure_tp(eos, 647.096, 220.6e5)
sc = state_critical.entropy_isov() / R
# d = state.density
# Dehlouz(state, param_scaling)
#%%

# T = np.linspace(300., 650., num = 10)
T = [416, 455.6, 494.44]

P = np.linspace(0.1, 1000., num = 100)

states = np.zeros((len(P), len(T)))
visc = np.zeros((len(P), len(T)))

#%%
for j, tj in enumerate(T):

    for i, pi in enumerate(P):

        state = State.pure_tp(eos, tj, pi * 1e5)
        
        # sres = state.entropy_isov()
        scaling = Dehlouz(state, param_scaling)
        visc[i,j] = scaling.viscosity(sc)
        

#%%
import pandas as pd
#%%


#%%
import matplotlib.cm as cm
import matplotlib.colors as colors

# normalização da temperatura → cor
# norm = colors.Normalize(vmin=T.min(), vmax=T.max())
# cmap = cm.jet   # parecido com o da figura

plt.figure()
fig, ax = plt.subplots()

for j,t in enumerate(T):
    
    plt.plot(P, visc[:,j] * 1000, label=f"{t} K")


list = [416, 455, 494]
for texp in list:
    s = "water" + str(texp) + ".csv"
    df = pd.read_csv(s, delimiter = "\t")
    pexp = df["Pressure (bar)"]
    viscexp = df["Viscosity (uPa*s)"] * 1e-3
    ax.scatter(pexp, viscexp)# barra de cores (o "termômetro")
# sm = cm.ScalarMappable(norm=norm, cmap=cmap)
# sm.set_array([])

# fig.colorbar(sm, ax=ax, label="T (K)")

ax.set_xlabel("P (bar)")
ax.set_ylabel("Viscosity (mPa.s)")

plt.title("Water")
plt.legend()
plt.show()
# %%
