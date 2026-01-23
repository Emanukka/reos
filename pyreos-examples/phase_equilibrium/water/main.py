#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from reos.state import State
from reos.consts import *
from reos.cpa import CPAParameters
from reos.eos import EquationOfState


xsize=3.15 * 1.2
ysize=3.15 * 1.2

plt.rcParams.update({
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 14,
    "font.size": 16,
    "legend.fontsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "figure.figsize": (xsize, ysize),  
    "text.usetex": True,

})

data = pd.read_csv("watersat_NIST.csv")


temperature = data["Temperature(K)"].to_numpy()[::2]
pressure = data["Pressure(bar)"].to_numpy()[::2]
liquid_density = data["Density(l;mol/m3)"].to_numpy()[::2]
vapor_density = data["Density(v;mol/m3)"].to_numpy()[::2]

SL = data["Entropy(l;J/mol*K)"].to_numpy()
SV = data["Entropy(v;J/mol*K)"].to_numpy()
HL = data["Enthalpy(l;kJ/mol)"].to_numpy()
HV = data["Enthalpy(v;kJ/mol)"].to_numpy()

dH = HV - HL
dS = SV - SL

dH = dH[::2]
dS = dS[::2]
R = Consts.ideal_gas_const()

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
p = CPAParameters.from_json(["water"], ppath="../../../parameters/cpa/kontogeorgis2006.json",)
eos = EquationOfState.scpa(p)

# antoine = np.array([water_antoine])
#%% Psat calc
def calc_psat(t,p0):
    
    e = 1
    max = 100
    it = 0

    while abs(e) > 1e-8 and it < max:

        s1 = State.tpx(eos,t,p0,np.array([1.0]),'vapor')
        s2 = State.tpx(eos,t,p0,np.array([1.0]),'liquid')
        phiv = np.exp(s1.lnphi()[0])
        phil = np.exp(s2.lnphi()[0])

        r = phil / phiv

        p0 = p0 * r
        it += 1

    return p0, s1, s2


#%% Calc properties

N = 100
T=np.linspace(250.0,650,N)

PRES= np.zeros_like(T)
VAP = np.zeros_like(T,dtype=object)
LIQ = np.zeros_like(T,dtype=object)

XV = np.zeros_like(VAP,dtype=object)
XL = np.zeros_like(LIQ,dtype=object)

rhoV = np.zeros_like(VAP)
rhoL = np.zeros_like(LIQ)

entropyV = np.zeros_like(VAP)
entropyL = np.zeros_like(LIQ)

for (i,t) in enumerate(T):
    
    PRES[i], VAP[i], LIQ[i] = calc_psat(t, PRES[i-1] if i>0 else 1e5)
    
    liq = LIQ[i]
    vap = VAP[i]

    # entropyL[i] = liq.entropy() 
    # entropyV[i] = vap.entropy() 
    entropyL[i] = liq.entropy() 
    entropyV[i] = vap.entropy() 

    rhoV[i] = vap.density
    rhoL[i] = liq.density

    XL[i] = eos.unbonded_sites_fraction(t, rhoL[i], np.array([1.0]))[0]
    XV[i] = eos.unbonded_sites_fraction(t, rhoV[i], np.array([1.0]))[0]



#%%

# plt.title("W4C")
plt.ylabel(r"$X_{A,B}$")
plt.xlabel(r"$\mathrm{T/K}$")
plt.scatter(Tliq, Xliq2005,marker='o',facecolors='none',edgecolors='black')

plt.plot(T[2:],XV[2:],linestyle="-",label="vap",color="black")
plt.plot(T, XL,linestyle="-.",label="liq" ,color="black")

plt.scatter(Tvap, Xvap2005,label="Dufal 2015",marker='o',facecolors='none',edgecolors='black')

plt.legend()

# plt.savefig("plots/water/xsat.png", bbox_inches='tight')

#%%

plt.ylabel(r"$\mathrm{P/bar}$")
plt.scatter(tsat,psat,marker='o',facecolors='none',edgecolors='black')
plt.xlabel(r"$\mathrm{T/K}$")
plt.plot(T[2:],PRES[2:]/1e5,color="black")

# plt.savefig("plots/water/psat.png", bbox_inches='tight')

# %%
to_kgm3 = 18.0153 / 1000

plt.ylabel(r"$\mathrm{Density / (kg/m^3)}$")
plt.xlabel("T/K")
plt.plot(T,rhoL * to_kgm3,linestyle="-",label="liq" ,color="black")
plt.plot(T[2:],rhoV[2:] * to_kgm3,linestyle="-.",label="vap",color="black")
plt.scatter(temperature[::2], liquid_density[::2] * to_kgm3,marker='o',facecolors='none',edgecolors='black',alpha=1)
plt.scatter(temperature[::2], vapor_density[::2] * to_kgm3,marker='o',facecolors='none',edgecolors='black',alpha=1)
plt.legend()

# plt.savefig("plots/water/rhosat.png", bbox_inches='tight')
# plt.savefig("plots/water_psat_X.pdf",bbox_inches='tight')

#%%

plt.ylabel(r" $\mathrm{\Delta H _{vap}/ kJ \; mol^{-1}}$")
plt.xlabel(r"$\mathrm{T/K}$")

plt.scatter(temperature, dH, marker='o',s=40, facecolors='none', edgecolors='black',alpha=0.9)

DS = entropyV - entropyL
DH = T * DS

plt.plot(T[2:], DH[2:] / 1000 , color="black")
# plt.savefig("plots/water/dhsat.png", bbox_inches='tight')

# plt.show()
#%%

plt.ylabel(r" $\mathrm{\Delta S _{vap}/ J \; (mol^{-1}} K^{-1})$")
plt.xlabel(r"$\mathrm{T/K}$")

plt.scatter(temperature, dS, marker='o',s=40 ,facecolors='none', edgecolors='black',alpha=0.9)

DS = entropyV - entropyL

plt.plot(T[2:], DS[2:]  , color="black")

# plt.savefig("plots/water/dssat.png", bbox_inches='tight')

# plt.show()

#%%


fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(10,8), sharex=True)
fig.tight_layout() # Or equivalently,  "plt.tight_layout()"
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.2, hspace= None)

temperature_label = r"$\mathrm{T / K}$"
pressure_label = r"$\mathrm{P / bar}$"
density_label = r"$\mathrm{Density / (kg \; m^{-3})}$"
entropy_label = r"$\mathrm{\Delta S_{vap} / J \; (mol^{-1}} K^{-1})$"
enthalpy_label = r"$\mathrm{\Delta H_{vap} / kJ \; mol^{-1}}$"
X_label = r"$X_{A,B}$"

#%%
ax = axs[0, 0]

ax.set_ylabel(pressure_label)
ax.set_xlabel(temperature_label)

ax.plot(T[2:], PRES[2:] / 1e5 , "-", color="black")

ax.scatter(
    temperature,pressure,
    marker='o', facecolors='none', edgecolors='black'
)

# ax.legend()
ax = axs[0, 1]

ax.set_ylabel(density_label)
ax.set_xlabel(temperature_label)

ax.plot(T, rhoL * to_kgm3, "-", color="black")
ax.plot(T[2:], rhoV[2:] * to_kgm3, "-.", color="black")

ax.scatter(
    temperature[::2],liquid_density[::2] * to_kgm3,
    marker='o', facecolors='none', edgecolors='black'
)

ax.scatter(
    temperature[::2], vapor_density[::2] * to_kgm3,
    marker='o', facecolors='none', edgecolors='black'
)


ax = axs[1, 0]

ax.set_ylabel(enthalpy_label)
ax.set_xlabel(temperature_label)

ax.plot(T[2:], DH[2:] / 1000 , "-", color="black")

ax.scatter(
    temperature,dH,
    marker='o', facecolors='none', edgecolors='black'
)



ax = axs[1, 1]

ax.set_ylabel(entropy_label)
ax.set_xlabel(temperature_label)

ax.plot(T[2:], DS[2:]  , "-", color="black")

ax.scatter(
    temperature,dS,
    marker='o', facecolors='none', edgecolors='black'
)

#%%
fig.savefig("plots/pdes.png", bbox_inches='tight')

# %%
plt.figure(figsize=(xsize,ysize))
plt.ylabel(r"$X_{A,B}$")
plt.xlabel(r"$\mathrm{T/K}$")
plt.scatter(Tliq, Xliq2005,marker='o',facecolors='none',edgecolors='black')

plt.plot(T[2:],XV[2:],linestyle="-",label="vap",color="black")
plt.plot(T, XL,linestyle="-.",label="liq" ,color="black")

plt.scatter(Tvap, Xvap2005,label="Dufal 2015",marker='o',facecolors='none',edgecolors='black')

plt.legend()

plt.savefig("plots/xsat.png", bbox_inches='tight', dpi = 300)

# %%
