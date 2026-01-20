#%%
from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
from phase_equilibrium.antoine import *
from auxiliary_functions.association_functions import *

import matplotlib.pyplot as plt
xsize=3.15
ysize=3.15
from auxiliary_functions.association_functions import get_non_bondend_sites_from_states
yes_or_no=False
plt.rcParams.update({
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": (xsize, ysize),  

})

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

#%%
p=CPAParameters.from_records(
    cubic=[c_w],
    assoc=[a_w])


antoine=np.array([water_antoine])
#%%
def calc_psat(t,p0):
    
    e=1
    max=100
    it=0
    while abs(e)>1e-8 and it<max:

        s1=State.tpx(eos,t,p0,np.array([1.0]),'vapor')
        s2=State.tpx(eos,t,p0,np.array([1.0]),'liquid')
        phiv=np.exp(s1.ln_phi())
        phil=np.exp(s2.ln_phi())

        r=phil/phiv

        p0=p0*r
        it+=1


    return p0,s1,s2

def res1(t,p0):
    s1=State.tpx(eos,t,p0,np.array([1.0]),'vapor')
    s2=State.tpx(eos,t,p0,np.array([1.0]),'liquid')
    phiv=np.exp(s1.ln_phi())
    phil=np.exp(s2.ln_phi())

    r=phil/phiv

    return phiv,phil
#%%
eos=EquationOfState.cpa(p)

T=np.linspace(250.0,650,50)

PRES=np.zeros_like(T)
VAP=np.zeros_like(T,dtype=object)
LIQ=np.zeros_like(T,dtype=object)


for (i,t) in enumerate(T):

    PRES[i],VAP[i],LIQ[i]=calc_psat(t,PRES[i-1] if i>0 else 1e5)
    # print(i,p)
    # print(res[i]/1e5)


#%%

#%%

plt.ylabel("P/bar")
plt.scatter(tsat,psat,marker='o',facecolors='none',edgecolors='black',label="NIST")
plt.xlabel("T/K")
plt.legend()
plt.plot(T[1:],PRES[1:]/1e5,color="black")
plt.savefig("plots/wpsat/water_psat.pdf",bbox_inches='tight')
# VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)
# plt.vlines(647.096,0,max(PRES/1e5),color="red",linestyle="--")

#%%

XV=np.zeros_like(VAP,dtype=object)
XL=np.zeros_like(LIQ,dtype=object)
rhoV=np.zeros_like(VAP)
rhoL=np.zeros_like(LIQ)
for i in range(len(LIQ)):

    XV[i]=VAP[i].non_bonded_sites()[0]
    XL[i]=LIQ[i].non_bonded_sites()[0]
    rhoV[i]=VAP[i].density()
    rhoL[i]=LIQ[i].density()
#%%

plt.title("W4C")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel("T/K")
plt.scatter(Tliq,Xliq2005,marker='o',facecolors='none',edgecolors='black',)

plt.plot(T,XL,linestyle="-",label="liq" ,color="black")
plt.scatter(Tvap,Xvap2005,label="Mol.Sim.",marker='o',facecolors='none',edgecolors='black')

plt.plot(T[1:],XV[1:],linestyle="--",label="vap",color="black")
plt.legend()
plt.savefig("plots/wpsat/water_psat_X_4C.pdf",bbox_inches='tight')
#%%

liq_monomer_4c=XL**4
vap_monomer_4c=XV**4
plt.title("W4C")
plt.ylabel("Monomer Fraction")
plt.xlabel("T/K")
plt.plot(T,liq_monomer_4c,linestyle="-",label="liq" ,color="black")
plt.plot(T,vap_monomer_4c,linestyle="--",label="vap",color="black")


plt.legend()
# plt.savefig("plots/wpsat/monomer_4c.pdf",bbox_inches='tight')

# %%

plt.xlabel("Density (m3/mol)")
plt.ylabel("T/K")
plt.plot(rhoL,T,linestyle="-",label="liq" ,color="black")
plt.plot(rhoV[1:],T[1:],linestyle="--",label="vap",color="black")
plt.legend()
# plt.savefig("plots/water_psat_X.pdf",bbox_inches='tight')

#%%
#%%

p=CPAParameters.from_records(
    cubic=[c_w3b],
    assoc=[a_w3b])

eos=EquationOfState.cpa(p)

T=np.linspace(298.15,660,50)

PRES=np.zeros_like(T)
VAP=np.zeros_like(T,dtype=object)
LIQ=np.zeros_like(T,dtype=object)


for (i,t) in enumerate(T):

    PRES[i],VAP[i],LIQ[i]=calc_psat(t,PRES[i-1] if i>0 else 1e5)
    # print(i,p)
    # print(res[i]/1e5)

#%%

plt.title("W3B ")
plt.ylabel("P/bar")
plt.xlabel("T/K")
plt.plot(T,PRES/1e5,color="black")
# plt.savefig("plots/water_psat.pdf",bbox_inches='tight')
# VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)
plt.vlines(647.096,0,max(PRES/1e5),color="red",linestyle="--")

#%%

XV=np.zeros_like(VAP,dtype=object)
XL=np.zeros_like(LIQ,dtype=object)
rhoV=np.zeros_like(VAP)
rhoL=np.zeros_like(LIQ)
for i in range(len(LIQ)):

    XV[i]=VAP[i].non_bonded_sites()
    XL[i]=LIQ[i].non_bonded_sites()
    rhoV[i]=VAP[i].density()
    rhoL[i]=LIQ[i].density()

XL=np.stack(XL).T
XV=np.stack(XV).T


#%%

plt.title("W3B")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel("T/K")
plt.plot(T,XL[0],linestyle="-",label="liq (-) site" )
plt.plot(T,XV[0],linestyle="--",label="vap (-) site")
plt.plot(T,XL[1],linestyle="-",label="liq (+) site" ,)
plt.plot(T,XV[1],linestyle="--",label="vap (+) site",)

plt.legend()
# plt.savefig("plots/water_psat_X.pdf",bbox_inches='tight')

#%%
liq_monomer_3b=XL[0]*(XL[1]**2)
vap_monomer_3b=XV[0]*(XV[1]**2)

plt.title("W3B")
plt.ylabel("Monomer Fraction")
plt.xlabel("T/K")
plt.plot(T,liq_monomer_3b,linestyle="-",label="liq" ,color="black")
plt.plot(T,vap_monomer_3b,linestyle="--",label="vap",color="black")
plt.legend()
# plt.savefig("plots/wpsat/monomer_4c.pdf",bbox_inches='tight')

#%%4c and 3b

plt.ylabel("Monomer Fraction")
plt.xlabel("T/K")
plt.plot(T,liq_monomer_3b,linestyle="-",label="liq3B" )
plt.plot(T,liq_monomer_4c,linestyle="-",label="liq4C" )

plt.plot(T,vap_monomer_3b,linestyle="--",label="vap3B")
plt.plot(T,vap_monomer_4c,linestyle="--",label="vap4C")

plt.legend()

# %%

plt.xlabel("Density (m3/mol)")
plt.ylabel("T/K")
plt.plot(rhoL,T,linestyle="-",label="liq" ,color="black")
plt.plot(rhoV,T,linestyle="--",label="vap",color="black")
plt.legend()

#%%
eos=EquationOfState.cpa(p)

# T=np.linspace(298.15,660,50)

PRES=np.linspace(1e5,1e6,50)
VAP=np.zeros_like(T,dtype=object)
LIQ=np.zeros_like(T,dtype=object)
# 

for (i,p) in enumerate(PRES):

    VAP[i],LIQ[i]=res1(400.15,p)
    # print(i,p)
    # print(res[i]/1e5)
# %%
# plt.plot(PRES,VAP)
# plt.plot(PRES,LIQ)
# plt.plot(PRES,abs(LIQ-VAP),"-")
# plt.plot(PRES,abs((LIQ/VAP)-1),"--")

