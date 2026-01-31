#%%

from reos.scaling import Dehlouz
from reos.state import State
from reos.eos import EquationOfState as E
from reos.cubic import CubicParameters
import numpy as np
import scipy.optimize as opt
from reos.consts import Consts
import pandas as pd

R = Consts.ideal_gas_const()
param_scaling = [-0.174896, 0.023392, 0.569631, 0.187410, 0.512607, 0.585792]
x = np.array([1.])


#%%


list = [250, 271, 313, 507]

#%%


# %%

t = 198.15
p = 1e5
parameters = CubicParameters.from_json(["n-hexane"], "./parameters.json", opt = "pr78" )
eos = E.cubic(parameters)
# state = State.pure_tp(eos, t, p)
# d = state.density
# Dehlouz(state, param_scaling)

#%%
state_critical = State.pure_tp(eos,  507.60, 30.25e5)
sc = -1.290042565
# sc = state_critical.entropy_isov() / R
# d = state.density
# Dehlouz(state, param_scaling)
#%%

T = np.linspace(250., 650., num = 20)

P = np.linspace(0.1, 700., num = 100)

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




#%%
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors

#%%
plt.scatter(pexp, viscexp)

#%%
# normalização da temperatura → cor
norm = colors.Normalize(vmin=T.min(), vmax=T.max())
cmap = cm.jet   # parecido com o da figura

plt.figure()
fig, ax = plt.subplots()

for j,t in enumerate(T):
    
    plt.plot(P, visc[:,j] * 1000, color=cmap(norm(t)))

# barra de cores (o "termômetro")
sm = cm.ScalarMappable(norm=norm, cmap=cmap)
sm.set_array([])

fig.colorbar(sm, ax=ax, label="T (K)")

ax.set_xlabel("P (bar)")
ax.set_ylabel("Viscosity (mPa.s)")

for texp in list:
    s = "hexane" + str(texp) + ".csv"
    df = pd.read_csv(s, delimiter = "\t")
    pexp = df["Pressure (bar)"]
    viscexp = df["Viscosity (uPa*s)"] * 1e-3
    ax.scatter(pexp, viscexp)

plt.show()

# def dp(P, d, h = 1e-5):

#     return (P(d + h) - P(d - h)) / (2 * h)

# def d2p(P, d, h = 1e-5):

#     return (P(d + h) - 2. * P(d) + P(d - h)) / (h**2.)


# def F(X):

#     t, p = np.exp(X)
#     print(t, p)
#     state = State.pure_tp(eos, t, p)
#     d = state.density
#     alias = lambda d: eos.pressure(t, d, x)

#     dpdr = dp(P = alias, d = d)
#     dpdr2 = d2p(P = alias, d = d)

#     dpdv = - dpdr * d**2
#     dpdv2 = dpdr2 * d**3
#     return np.array([dpdv, dpdv2])

# # %%
# tcexp =  507.60 
# pcexp = 30.25 * 1e5 

# X0 = np.array([tcexp, pcexp])
# sol = opt.root(F, np.log(X0))

# #%%

# t, p = np.exp(sol.x)

# # %%

# %%
