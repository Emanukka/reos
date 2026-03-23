#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from reos.state import State
from reos.consts import *
# from reos.cpa import CPAParameters
from reos.cubic import CubicParameters, CubicPureRecord

from reos.eos import EquationOfState

from si_units import RGAS, MOL, JOULE, KELVIN
# import argparse
# parser = argparse.ArgumentParser()
# parser.add_argument("-s","--save", action="store_true")
# parser.add_argument("-p","--print", action="store_true")
# args = parser.parse_args()
# SAVE = args.save
# PRINT = args.print

xsize = 3.15 * 1.2
ysize = 3.15 * 1.2

plt.rcParams.update({
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 14,
    "font.size": 16,
    "legend.fontsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "figure.figsize": (xsize, ysize),
    # "figure.dpi": 300,  
    "text.usetex": True,

})

data = pd.read_csv("octane.csv", sep="\t")


temperature = data["Temperature(K)"].to_numpy()[::2]
pressure = data["Pressure(bar)"].to_numpy()[::2]
liquid_density = data["Density(l,mol/m3)"].to_numpy()[::2]
vapor_density = data["Density(v,mol/m3)"].to_numpy()[::2]

# SL = data["Entropy(l,J/mol*K)"].to_numpy()
# SV = data["Entropy(v,J/mol*K)"].to_numpy()
# HL = data["Enthalpy(l,kJ/mol)"].to_numpy()
# HV = data["Enthalpy(v,kJ/mol)"].to_numpy()

# dH = HV - HL
# dS = SV - SL

# dH = dH[::2]
# dS = dS[::2]
R = RGAS * MOL * KELVIN / JOULE

#%%Initializing


comp1 = {
    "name":"octane",
    "tc":568.70,
    "pc": 24.90 * 1e5,
    "l": 0.3385,
    "m": 0.8185, 
    "n": 2.0747,
    "c":6.4134 / 1e6
    
}
comp2 = {
    "name":"octane",
    "tc":568.70,
    "pc": 24.90 * 1e5,
    "l": 0.3385,
    "m": 0.8185, 
    "n": 2.0747,
}

   
parameter1 = CubicParameters.from_records(
    [
        CubicPureRecord.new(**{
            "name":"octane",
            "tc":568.70,
            "pc": 24.90 * 1e5,
            "l": 0.3385,
            "m": 0.8185, 
            "n": 2.0747,
            "c":6.4134 / 1e6
        }) 
    ], 
    cubic_model = "pr78", alpha_model = "twu91")

parameter2 = CubicParameters.from_records(
    [
        CubicPureRecord.new(**{
            "name":"octane",
            "tc":568.70,
            "pc": 24.90 * 1e5,
            "l": 0.3385,
            "m": 0.8185, 
            "n": 2.0747,
        }) 
    ], 
    cubic_model = "pr78", alpha_model = "twu91")

models = [
    EquationOfState.cubic(parameter1),
    EquationOfState.cubic(parameter2)
]
# eos = EquationOfState.cubic(parameters)

#%%

v_min = (0.00014773239768095766 - 6.4134 / 1e6 )
# v0 = 1.8970510424824202e-5



vv = np.logspace(np.log10(1.01 * v_min), np.log10(R*250/1e7) ,100)
# vv = np.linspace(v0, R*300.0/1e5,100)

vp = np.zeros_like(vv)
for model in [models[0]]:
    for i, v in enumerate(vv):
        vp[i] = model.pressure(250, 1 / v , np.array([1.]))
        print(i, vp[i])

plt.hlines(1, 1/vv[0], 1/vv[-1], linestyles="dashed", color="black")
# plt.vlines(1.1 *1.8970510424824202e-5, -4000, 6000, linestyles="dashed", color="black")
plt.vlines(1 / (1.01*v_min), -1000, 1000, linestyles="dashed", color="b", label="1.01 * v_min")
plt.vlines(1 / (1.1*v_min), -1000, 1000, linestyles="dashed", color="r", label="1.1 * v_min")
plt.vlines(0.99 * (0.9/v_min), -1000, 1000, linestyles="dashed", color="g", label="0.99 * ( 0.9/v_min)")

# plt.hlines(1, 1e-5, 1e-3, linestyles="dashed", color="black")

plt.ylim(-10, 10)
plt.xlim(6000, 7 * 1e3)

plt.semilogx(1/vv, vp / 1e5, "-")
plt.legend()
# plt.plot(vv, vp / 1e5, "-")

#%%

v_min = (0.00014773239768095766 - 6.4134 / 1e6 )
# v0 = 1.8970510424824202e-5
rho_max = 1 / v_min


vv = np.logspace(np.log10(1.01 * v_min), np.log10(R*250/1e5) ,100)
# vv = np.linspace(v0, R*300.0/1e5,100)

rhov = 1 / vv

s = rhov/rho_max
vp = np.zeros_like(vv)
for model in [models[0]]:
    for i, v in enumerate(vv):
        vp[i] = model.pressure(250, 1 / v , np.array([1.]))
        print(i, vp[i])

# plt.hlines(1, 1/vv[0], 1/vv[-1], linestyles="dashed", color="black")
# # plt.vlines(1.1 *1.8970510424824202e-5, -4000, 6000, linestyles="dashed", color="black")
# plt.vlines(1 / (1.01*v_min), -1000, 1000, linestyles="dashed", color="b", label="1.01 * v_min")
# plt.vlines(1 / (1.1*v_min), -1000, 1000, linestyles="dashed", color="r", label="1.1 * v_min")

plt.vlines(0.99 , -1e7, 1e7, linestyles="dashed", color="g")

# plt.hlines(, 1e-5, 1e-3, linestyles="dashed", color="black")

# plt.ylim(-0.5e7, 0.5e7)
plt.hlines(0, 0, 1, linestyles="dashed", color="black")
plt.xlim(0.5, 1)

plt.semilogx(s, (1-s) * (vp - 1e5), "-")
plt.legend()

#%%

# plt.scatter(range(len(vv)), vv)
#%% Psat calc

def calc_psat(eos, t, p0):
    
    e = 1
    it = 0

    while abs(e) > 1e-8 and it < 100:

        s1 = State.tpx(eos, t, p0, np.array([1.0]), 'vapor')
        s2 = State.tpx(eos, t, p0, np.array([1.0]), 'liquid')
        phiv = np.exp(s1.lnphi()[0])
        phil = np.exp(s2.lnphi()[0])
        r = phil / phiv
        e = 1.0 - r
        p0 = p0 * r
        it += 1

    # print(it)
    return p0, s1, s2



#%% Calc properties

N = 100
T=np.linspace(230, 553,N)

PRES= np.zeros((N, len(models)))
VAP = np.zeros_like(PRES,dtype=object)
LIQ = np.zeros_like(PRES,dtype=object)


rhoV = np.zeros_like(VAP)
rhoL = np.zeros_like(LIQ)

# entropyV = np.zeros_like(VAP)
# entropyL = np.zeros_like(LIQ)

for j,model in enumerate(models):
    for (i,t) in enumerate(T):
        
    
        try:
            PRES[i,j], VAP[i,j], LIQ[i,j] = calc_psat(model, t, PRES[i-1, j] if i>0 else 1e-5)

            liq = LIQ[i, j]
            vap = VAP[i, j]

            # entropyL[i] = liq.entropy() 
            # entropyV[i] = vap.entropy() 
            # entropyL[i] = liq.tp_entropy() 
            # entropyV[i] = vap.tp_entropy() 

            rhoV[i, j] = vap.density
            rhoL[i, j] = liq.density

        except Exception as e:
            print(f"Error at T={t} K: {e}")
            PRES[i, j] = np.nan
            VAP[i, j] = None
            LIQ[i, j] = None
            rhoV[i, j] = np.nan
            rhoL[i, j] = np.nan

        # XL[i, j] =  eos.get_assoc_calcs(t, rhoL[i, j], np.array([1.0]))["X"][0]
        # XV[i, j] =  eos.get_assoc_calcs(t, rhoV[i, j], np.array([1.0]))["X"][0]

# DS = entropyV - entropyL
# DH = T * DS
to_kgm3 = 18.0153 / 1000

#%%

temperature_label = r"$\mathrm{T / K}$"
pressure_label = r"$\mathrm{P / bar}$"
density_label = r"$\mathrm{Density / (kg \; m^{-3})}$"
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(8,4), sharex=True)

plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace= None)
ax = axs[0]

ax.set_ylabel(pressure_label)
ax.set_xlabel(temperature_label)
ax.plot(T, PRES[:, 0] / 1e5, "-", color="b")
ax.plot(T, PRES[:, 1] / 1e5, "--", color="r")

ax.scatter(
    temperature[:],pressure[:],
    marker='o', facecolors='none', edgecolors='black'
)

ax = axs[1]

ax.set_ylabel(density_label)
ax.set_xlabel(temperature_label)
ax.plot(T, rhoL[:, 0] * to_kgm3, "-", color="r", label="tPR78")
ax.plot(T, rhoL[:, 1] * to_kgm3, "--", color="b", label="PR78")


ax.plot(T, rhoV[:, 1] * to_kgm3, "-", color="r", label="PR78")
ax.plot(T, rhoV[:, 0] * to_kgm3, "--", color="b", label="tPR78")

ax.scatter(
    temperature[:],liquid_density[:] * to_kgm3,
    marker='o', facecolors='none', edgecolors='black'
)

plt.legend()

#%% Plot
# fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,8), sharex=True)
# fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
# plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace= None)


# # entropy_label = r"$\mathrm{\Delta S_{vap} / J \; (mol^{-1}} K^{-1})$"
# # enthalpy_label = r"$\mathrm{\Delta H_{vap} / kJ \; mol^{-1}}$"
# # X_label = r"$X_{A,B}$"


# #%%
# ax = axs[0, 0]

# ax.set_ylabel(pressure_label)
# ax.set_xlabel(temperature_label)

# ax.plot(T[2:], PRES[2:] / 1e5 , "-", color="black")

# ax.scatter(
#     temperature,pressure,
#     marker='o', facecolors='none', edgecolors='black'
# )

# ax = axs[0, 1]

# ax.set_ylabel(density_label)
# ax.set_xlabel(temperature_label)

# ax.plot(T, rhoL * to_kgm3, "-", color="black")
# ax.plot(T[2:], rhoV[2:] * to_kgm3, "-.", color="black")

# ax.scatter(
#     temperature[::2],liquid_density[::2] * to_kgm3,
#     marker='o', facecolors='none', edgecolors='black'
# )

# ax.scatter(
#     temperature[::2], vapor_density[::2] * to_kgm3,
#     marker='o', facecolors='none', edgecolors='black'
# )


# ax = axs[1, 0]

# ax.set_ylabel(enthalpy_label)
# ax.set_xlabel(temperature_label)

# ax.plot(T[2:], DH[2:] / 1000 , "-", color="black")

# ax.scatter(
#     temperature,dH,
#     marker='o', facecolors='none', edgecolors='black'
# )

# ax = axs[1, 1]

# ax.set_ylabel(entropy_label)
# ax.set_xlabel(temperature_label)

# ax.plot(T[2:], DS[2:]  , "-", color="black")

# ax.scatter(
#     temperature,dS,
#     marker='o', facecolors='none', edgecolors='black'
# )



# if SAVE:
    # fig.savefig("plots/pdes.png", bbox_inches='tight', dpi = 300)

# %%
# plt.figure(figsize=(xsize,ysize))
# plt.ylabel(r"$X_{A,B}$")
# plt.xlabel(r"$\mathrm{T/K}$")
# plt.scatter(Tliq, Xliq2005,marker='o',facecolors='none',edgecolors='black')

# plt.plot(T[2:],XV[2:],linestyle="-",label="vap",color="black")
# plt.plot(T, XL,linestyle="-.",label="liq" ,color="black")

# plt.scatter(Tvap, Xvap2005,label="Dufal 2015",marker='o',facecolors='none',edgecolors='black')

# plt.legend()

# if SAVE:

    # plt.savefig("plots/xsat.png", bbox_inches='tight', dpi = 300)

# if PRINT:
if True:
    fig1 = plt.figure(1)
    # fig2 = plt.figure(2)
    wfig2 = 1.2 * fig1.get_figwidth() * fig1.get_dpi()
    fig1.canvas.manager.window.wm_geometry("+%d+%d" % (100, 100))
    fig2.canvas.manager.window.wm_geometry("+%d+%d" % (wfig2, 400))
    fig1.tight_layout()
    fig2.tight_layout()
    plt.show()