#%%

from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

from auxiliary_functions.parameters import *
from auxiliary_functions.vle_functions import *
from auxiliary_functions.data import *


#%%
pACOH_OCT=CPAParameters.from_records(
    cubic=[c_acoh,c_octane],
    assoc=[a_acoh,a_octane])

pACOH_OCT.set_cubic_binary(0,1,0.0, 0.064)

ACOH_OCT=EquationOfState.cpa(pACOH_OCT)# print(pWATER_ACETIC.as_string())

# peq=PhaseEquilibrium(ACOH_OCT)

antoine=np.array([acoh_antoine,octane_antoine])
# p,y,vx=vle_diagram(T,peq,factor=1e3)

xorv,porv=acoh_octane["orv"]
xbol,pbol=acoh_octane["bol"]
exp_data=[xorv,porv,xbol,pbol]

T=343.2

PRES,LIQUID,VAPOR=linspace_bubble_p(ACOH_OCT,T,antoine,N=500)

bubble_diagram(
   PRES,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=15,
   y_sup=30,
   text=f"{T}K",
   title="Acetic Acid 1A(1) and Octane(2)",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   exp_data=exp_data)

#%%
xascL=np.zeros_like(LIQUID,dtype=object)
xascV=np.zeros_like(VAPOR,dtype=object)
for i,t in enumerate(PRES):

    xascL[i]=LIQUID[i].non_bonded_sites()
    
    xascV[i]=VAPOR[i].non_bonded_sites()

plt.plot(xascL,PRES/1e3)
plt.plot(xascV,PRES/1e3)
# plt.yscale('log')
#%%

T=323.15
# VLE_DIAGRAM(("t",T),peq,antoine,y_label="P/kPa",x_label="x1,y1(AcOH)",y_lim=[15,31],x_figsize=8,factor=1e3,title="AcOH(1A)&Octane",exp_data=[xorv,porv,xbol,pbol],N_points=100)
pPROPANOIC_HEPTANE=CPAParameters.from_records(
    cubic=[c_propanoic,c_heptane],
    assoc=[a_propanoic,a_heptane])

pPROPANOIC_HEPTANE.set_cubic_binary(0,1,0.0, 0.017)

PROPANOIC_HEPTANE=EquationOfState.cpa(pPROPANOIC_HEPTANE)
# print(pWATER_ACETIC.as_string())

peq=PhaseEquilibrium(PROPANOIC_HEPTANE)
antoine=np.array([propanoic_antoine,heptane_antoine])
# p,y,vx=vle_diagram(T,peq,factor=1e3)
xorv,porv=propanoic_hep["orv"]
xbol,pbol=propanoic_hep["bol"]
exp_data=[xorv,porv,xbol,pbol]

VAR,LIQUID,VAPOR=linspace_bubble_p(PROPANOIC_HEPTANE,T,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=0,
   y_sup=25,
    text=f"{T}K",

   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   title="Propanoic Acid(1) and Heptane",
   exp_data=exp_data)



#%%
p=CPAParameters.from_records(
    cubic=[c_methanol_2b,c_octanol_2b],
    assoc=[a_methanol_2b,a_octanol_2b])

p.set_cubic_binary(0,1,0.0,0.0)
p.set_assoc_binary(0,1,"cr1")
eos_2b=EquationOfState.cpa(p)


antoine=np.array([metoh_antoine,octanol_antoine])
P=101.32e3

xorv,torv=metoh_otctanol["orv"]
xbol,tbol=metoh_otctanol["bol"]

exp_data=[xorv,torv,xbol,tbol]


VAR,LIQUID,VAPOR=linspace_bubble_t(eos_2b,P,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   title="Methanol(1) and Octanol(2) 2B+CR1 ",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
   text=f"{P/1e3}kPa",

   factor=1.0,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=300,y_sup=480,
   exp_data=exp_data)


#%%
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
#    orv_linestyle="--",
#    bol_linestyle="--",
   text=f"{P/1e3}kPa",

   factor=1.0,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=300,y_sup=480,
   exp_data=exp_data)

#%%
#Fig. 7, Application of the CPA equation of state to organic acids

pPROPANOIC_HEPTANE=CPAParameters.from_records(
    cubic=[c_propanoic,c_heptane],
    assoc=[a_propanoic,a_heptane])

pPROPANOIC_HEPTANE.set_cubic_binary(0,1,0.0,  0.029)

PROPANOIC_HEPTANE=EquationOfState.cpa(pPROPANOIC_HEPTANE)
# print(pWATER_ACETIC.as_string())

peq=PhaseEquilibrium(PROPANOIC_HEPTANE)
# p,y,vx=vle_diagram(T,peq,factor=1e3)

antoine=np.array([propanoic_antoine,heptane_antoine])
propanoic_hep_teb
xorv,torv=propanoic_hep_teb["orv"]
xbol,tbol=propanoic_hep_teb["bol"]
exp_data=[xorv,torv,xbol,tbol]
P=101.33e3



VAR,LIQUID,VAPOR=linspace_bubble_t(PROPANOIC_HEPTANE,P,antoine,N=100)




#%%
xascL=np.zeros_like(LIQUID,dtype=object)
xascV=np.zeros_like(VAPOR,dtype=object)
for i,t in enumerate(VAR):

    xascL[i]=LIQUID[i].non_bonded_sites()
    xascV[i]=VAPOR[i].non_bonded_sites()

plt.plot(VAR,xascL)
plt.plot(VAR,xascV)


#%%
# bubble_diagram(
#    VAR,
#    LIQUID,
#    VAPOR,
#    factor=1.0,
#    y_figsize=2.5,
#    x_figsize=5,
#    y_inf=360,
#    y_sup=420,
#    text=f"{P/1e3}kPa",
#    y_label="T/K",
#    x_label=r"$x_1,y_1$",
#    title="Propanoic Acid 1A(1) and Heptane(2)",
#    exp_data=exp_data)

# tsat propanoic estranho
#%%

pMETHANOL_ACETIC=CPAParameters.from_records(
    cubic=[c_methanol_2b,c_acoh],
    assoc=[a_methanol_2b,a_acoh])

pMETHANOL_ACETIC.set_cubic_binary(0,1,0.0,-0.04)
pMETHANOL_ACETIC.set_assoc_binary(0,1,"ecr")
METHANOL_ACETIC=EquationOfState.cpa(pMETHANOL_ACETIC)


antoine=np.array([metoh_antoine,acoh_antoine])
T=308.15

xorv,porv=metoh_acoh["orv"]
xbol,pbol=metoh_acoh["bol"]
exp_data=[xorv,porv,xbol,pbol]

PRES,LIQUID,VAPOR=linspace_bubble_p(METHANOL_ACETIC,T,antoine,N=100)

bubble_diagram(
   PRES,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=0.0,
   y_sup=30,
   text=f"{T}K",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   title="Methanol 2B(1) and AcOH 1A(2) (ECR)",
   exp_data=exp_data)


# BOL,_,linspaceZ,XASCL,XASCV=VLE_DIAGRAM(
#     ("t",T),
#     METHANOL_ACETIC,antoine,
#     y_label=r"P/kPa",
#     x_label=r"$x_1,y_1$",
#     y_lim=[0,30],
#     x_figsize=5,
#     y_figsize=5,
#     factor=1e3,
#     exp_data=exp_data,
#     title=r"Methanol 2B(1) and AcOH 1A(2) (ECR)",
#     save_fig=False,
#     N_points=100)

# #now testando modelos de grafico de X
# XASCV=np.stack(XASCV)
# XASCL=np.stack(XASCL)

# plt.figure()
# # plt.plot(BOL,XASC)
# # plt.plot(linspaceZ,XASCV[:,0],label="- MetOH VAP")
# # plt.plot(linspaceZ,XASCV[:,1],label="+ MetOH VAP")
# # plt.plot(linspaceZ,XASCV[:,2],label="+-AcOH  VAP")

# # plt.plot(linspaceZ,XASCL[:,0],label="-MetOH LIQ")
# # plt.plot(linspaceZ,XASCL[:,1],label="+MetOH LIQ")
# # plt.plot(linspaceZ,XASCL[:,2],label="+-AcOH LIQ")

# plt.plot(BOL,XASCV[:,0],label="- MetOH VAP")
# plt.plot(BOL,XASCV[:,1],label="+ MetOH VAP")
# plt.plot(BOL,XASCV[:,2],label="+-AcOH  VAP")
# plt.plot(BOL,XASCL[:,0],label="-MetOH LIQ")
# plt.plot(BOL,XASCL[:,1],label="+MetOH LIQ")
# plt.plot(BOL,XASCL[:,2],label="+-AcOH LIQ")

# plt.legend()
# # plt.xlim(-0.01,1.01)
# plt.ylim(-0.01,1.01)
