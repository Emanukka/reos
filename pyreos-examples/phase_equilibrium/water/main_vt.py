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

data = pd.read_csv("watersat_NIST.csv")


temperature = data["Temperature(K)"].to_numpy()[::2]
pressure = data["Pressure(bar)"].to_numpy()[::2]
liquid_density = data["Density(l;mol/m3)"].to_numpy()[::2]
vapor_density = data["Density(v;mol/m3)"].to_numpy()[::2]

# SL = data["Entropy(l;J/mol*K)"].to_numpy()
# SV = data["Entropy(v;J/mol*K)"].to_numpy()
# HL = data["Enthalpy(l;kJ/mol)"].to_numpy()
# HV = data["Enthalpy(v;kJ/mol)"].to_numpy()

# dH = HV - HL
# dS = SV - SL

# dH = dH[::2]
# dS = dS[::2]

R = RGAS * MOL * KELVIN / JOULE

# print(R)
#%%

Tliq=np.array([

252,
270.00,
290.00,
300.00,
310.00,
330.00,
350.00,
370.00,
390.0,
410.00,
430.00,
450.00,
470.00,
490.00,
510.00,
530.00,
550.00,
570.00,
590.00,
600.00,
610.00,
620.00,
])

Tvap=np.array([
370.0,
390.0,
430,
470.0,
550,
590,
600,
610,
620
])

Xvap2005=np.array([
0.995,
0.988,
0.977,
0.961,
0.902,
0.863,
0.836,
0.821,
0.784

])

Xliq2005=np.array([

0.025,
0.044,
0.063,
0.069,
0.077,

0.100,

0.120,
0.141,

0.162,
0.189,
0.210,
0.237,
0.262,
0.288,
0.315,
0.345,
0.374,
0.407,
0.444,
0.464,
0.483,
0.506,
])


tsat=np.array([
273.16,
283.16,
293.16,
303.16,
313.16,
323.16,
333.16,
343.16,
353.16,
363.16,
373.16,
383.16,
393.16,
403.16,
413.16,
423.16,
433.16,
443.16,
453.16,
463.16,
473.16,
483.16,
493.16,
503.16,
513.16,
523.16,
533.16,
543.16,
553.16,
563.16,
573.16,
583.16,
593.16,
603.16,
613.16,
623.16,
633.16,
643.16,])

psat=np.array([
0.0061165,
0.01229,
0.023408,
0.042494,
0.073889,
0.12358,
0.19956,
0.31214,
0.47434,
0.70208,
1.0145,
1.4343,
1.9874,
2.7036,
3.6164,
4.7629,
6.1839,
7.9238,
10.03,
12.555,
15.553,
19.081,
23.2,
27.976,
33.475,
39.768,
46.93,
55.038,
64.176,
74.429,
85.891,
98.664,
112.86,
128.6,
146.03,
165.31,
186.68,
210.46,
])

#%%Initializing


water1 = {
    "name":"water",
    "tc":647.10,
    "pc": 220.64 * 1e5,
    "l":0.3865,
    "m":0.8720, 
    "n":1.9693, 
    "c":5.3041 / 1e6
    
}
water2 = {
    "name":"water",
    "tc":647.10,
    "pc": 220.64 * 1e5,
    "l":0.3865,
    "m":0.8720, 
    "n":1.9693, 
    
}

parameters1 = CubicParameters.from_records([CubicPureRecord.new(**water1)], cubic_model = "pr78", alpha_model = "twu91")
parameters2 = CubicParameters.from_records([CubicPureRecord.new(**water2)], cubic_model = "pr78", alpha_model = "twu91")

# model1 = EquationOfState.cubic(parameters1)

# pr = CubicPureRecord.new(**water2)
# parameters2 = CubicParameters.from_records([pr], cubic_model = "pr78", alpha_model = "twu91")

# model2 = EquationOfState.cubic(parameters2)

# parameters = CPAParameters.from_json(["water"], 
#                                      ppath="../../../parameters/cpa/kontogeorgis2006.json", 
#                                      rdf_model="kg", 
#                                      cubic_model="srk")

# eos = EquationOfState.cpa(parameters)

models = [
    EquationOfState.cubic(parameters1),
    EquationOfState.cubic(parameters2),

    
]

#%%

v0 = 1.8970510424824202e-5 - 5.3041 / 1e6 
# v0 = 1.8970510424824202e-5



vv = np.logspace(np.log10(v0), np.log10(R*300.0/1e5) ,100)
# vv = np.linspace(v0, R*300.0/1e5,100)

vp = np.zeros_like(vv)
# for model in [eos]:
for i, v in enumerate(vv):
    vp[i] = models[0].pressure(300.0, 1 / v , np.array([1.]))
    print(i, vp[i])

#%%
# plt.hlines(vp[0], v0, vv[-1], linestyles="dashed", color="black")
# plt.vlines(1/(1.1 *1.8970510424824202e-5), -4000, 6000, linestyles="dashed", color="black")
plt.vlines(40.971563, -4000, 6000, linestyles="dashed", color="black")
plt.vlines(52713.394506, -4000, 6000, linestyles="dashed", color="black")

plt.hlines(1, 10, 1e5, linestyles="dashed", color="black")

plt.ylim(-10, 10)
plt.xlim(1e1, 1e5)
plt.xlim(1e1, 1e2)

plt.semilogx(1/vv, vp / 1e5, "-")
# plt.plot(vv, vp / 1e5, "-")

#%%
max_density = 1 / v0


s = (1 / vv) / max_density

f_obj = (1 - s) * (vp - 1e5)

plt.plot(s, f_obj / 1e8)
plt.hlines(0, 0, 1, linestyles="dashed", color="black")
plt.ylim(-1.75, 0.5)
# plt.scatter(range(len(vv)), vv)
#%% Psat calc

def calc_psat(eos, t, p0):
    
    e = 1
    it = 0
    print("psat")
    while abs(e) > 1e-8 and it < 100:

        print(it, p0, e)
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

T=np.linspace(300.0,640,N)

PRES= np.zeros((N, len(models)))
VAP = np.zeros_like(PRES,dtype=object)
LIQ = np.zeros_like(PRES,dtype=object)

rhoV = np.zeros_like(VAP)
rhoL = np.zeros_like(LIQ)

# entropyV = np.zeros_like(VAP)
# entropyL = np.zeros_like(LIQ)

for j,model in enumerate(models):
    for (i,t) in enumerate(T):
    
        print(i, t, PRES[i-1,j] if i>0 else 1e5)
        PRES[i,j], VAP[i,j], LIQ[i,j] = calc_psat(model, t, PRES[i-1,j] if i>0 else 1e5)
        
        # liq = LIQ[i,j]
        # vap = VAP[i,j]

        # entropyL[i] = liq.entropy() 
        # entropyV[i] = vap.entropy() 
        # entropyL[i] = liq.tp_entropy() 
        # entropyV[i] = vap.tp_entropy() 

        rhoV[i,j] = VAP[i,j].density
        rhoL[i,j] = LIQ[i,j].density
   

# DS = entropyV - entropyL
# DH = T * DS
to_kgm3 = 18.0153 / 1000

#%% Plot
fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(6,5), sharex=True)
fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.5, hspace= None)

temperature_label = r"$\mathrm{T / K}$"
pressure_label = r"$\mathrm{P / bar}$"
density_label = r"$\mathrm{Density / (kg \; m^{-3})}$"
# entropy_label = r"$\mathrm{\Delta S_{vap} / J \; (mol^{-1}} K^{-1})$"
# enthalpy_label = r"$\mathrm{\Delta H_{vap} / kJ \; mol^{-1}}$"
# X_label = r"$X_{A,B}$"

ax = axs[0]

ax.set_ylabel(pressure_label)
ax.set_xlabel(temperature_label)

ax.plot(T, PRES[:, 0] / 1e5 , "-", color="r",label="tPR78-twu91")
ax.plot(T, PRES[:, 1] / 1e5 , "--", color="b", label="PR78-twu91")


ax.scatter(
    temperature,pressure,
    marker='o', facecolors='none', edgecolors='black'
)

ax = axs[1]

ax.set_ylabel(density_label)
ax.set_xlabel(temperature_label)

ax.plot(T, rhoL[:, 0] * to_kgm3, "-", color="r")
ax.plot(T, rhoL[:, 1] * to_kgm3, "--", color="b")

ax.plot(T, rhoV[:, 0] * to_kgm3, "-", color="r")
ax.plot(T, rhoV[:, 1] * to_kgm3, "--", color="b")

ax.scatter(
    temperature,liquid_density * to_kgm3,
    marker='o', facecolors='none', edgecolors='black'
)

ylim=ax.get_ylim()

ax.vlines(0.8 * 647.1, ylim[0],ylim[1], linestyles="dashed", color="black")
# plt.legend()
# ax.set_ylim(0, 1200)
plt.savefig("figs/water_vt_psat_rhosat.png")
plt.show()
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