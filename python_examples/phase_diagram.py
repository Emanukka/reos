#%% main

from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

from auxiliary_functions.parameters import *
from auxiliary_functions.vle_functions import *
# from python_examples.auxiliary_functions.exp_data import *
from auxiliary_functions.exp_data import *
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




#%% dois
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

PRES,LIQUID,VAPOR=linspace_bubble_p(ACOH_OCT,T,antoine,N=100)

bubble_diagram(
   PRES,
   LIQUID,
   VAPOR,
   plot_dir="plots/acetic_octane",
   factor=1e3,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=15,
   y_sup=30,
   text=f"{T}K",
   title="(1) Acetic Acid 1A and (2) Octane",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   save_fig=True,
   

   exp_data=exp_data)

#%% plt x
N=100
print(N)
xascL=np.zeros_like(LIQUID,dtype=object)
xascV=np.zeros_like(VAPOR,dtype=object)
for i,t in enumerate(PRES):

    xascL[i]=LIQUID[i].non_bonded_sites()
    
    xascV[i]=VAPOR[i].non_bonded_sites()


plt.figure(figsize=(xsize,ysize))

sites=["-+ AcOH"]
plt.xlabel(r"$x_1,y_1$") 

plt.plot([LIQUID[i].composition()[0] for i in range(N)],xascL,label="Liquid",color="black",linestyle='-')
plt.plot([VAPOR[i].composition()[0] for i in range(N)],xascV ,label="Vapor",color="black",linestyle='--')
plt.ylabel(r"$X_{C,AcOH}$")
# plt.text(0.8,0.8,"T=343.2K")
plt.title("(1) Acetic Acid 1A and (2) Octane")

plt.legend()
#plt.tight_layout()
plt.savefig("plots/acetic_octane/X.pdf",bbox_inches='tight') 

plt.show()



# plt.yscale('log')
#

# T=323.15
# # VLE_DIAGRAM(("t",T),peq,antoine,y_label="P/kPa",x_label="x1,y1(AcOH)",y_lim=[15,31],x_figsize=8,factor=1e3,title="AcOH(1A)&Octane",exp_data=[xorv,porv,xbol,pbol],N_points=100)
# pPROPANOIC_HEPTANE=CPAParameters.from_records(
#     cubic=[c_propanoic,c_heptane],
#     assoc=[a_propanoic,a_heptane])

# pPROPANOIC_HEPTANE.set_cubic_binary(0,1,0.0, 0.017)

# PROPANOIC_HEPTANE=EquationOfState.cpa(pPROPANOIC_HEPTANE)
# # print(pWATER_ACETIC.as_string())

# antoine=np.array([propanoic_antoine,heptane_antoine])
# # p,y,vx=vle_diagram(T,peq,factor=1e3)
# xorv,porv=propanoic_hep["orv"]
# xbol,pbol=propanoic_hep["bol"]
# exp_data=[xorv,porv,xbol,pbol]

# VAR,LIQUID,VAPOR=linspace_bubble_p(PROPANOIC_HEPTANE,T,antoine,N=100)

# bubble_diagram(
#    VAR,
#    LIQUID,
#    VAPOR,
#    factor=1e3,
#    y_figsize=2.5,
#    x_figsize=5,
#    y_inf=0,
#    y_sup=25,
#     text=f"{T}K",

#    y_label="P/kPa",
#    x_label=r"$x_1,y_1$",
#    title="Propanoic Acid(1) and Heptane",
#     save_fig=yes_or_no,

#    exp_data=exp_data)



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
   plot_dir="plots/metanol_octanol",
   title="(1) Methanol 2B and (2) Octanol 2B CR1 ",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
   text=f"{P/1e3}kPa",

   factor=1.0,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=300,y_sup=480,
   save_fig=True,
   exp_data=exp_data)

#%%
xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]


#%%
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-MeOH","+MeOH","-OcOH","+OcOH"]
linestyle=['-','--','-.',':']

#%%

plt.figure(figsize=(xsize,ysize))

for (i,s) in enumerate(sites):

    plt.plot(xL,XL[:,i],label=s,linestyle=linestyle[i],color="black")

plt.ylabel("Fraction of Non-Bonded Sites")
plt.title("(1) Methanol 2B and (2) Octanol 2B CR1 Liquid Phase",fontsize=10)

# plt.savefig("metanol_octanol/X2B_liq.pdf") 
plt.xlabel(r"$x_1$") 

plt.legend()

#%%

plt.figure(figsize=(xsize,ysize))

for (i,s) in enumerate(sites):
    plt.plot(xV,XV[:,i],label=s)
plt.xlabel(r"$y_1$") 

plt.ylabel("Fraction of Non-Bonded Sites")
plt.title("(1) Methanol 2B and (2) Octanol 2B CR1  Vapor Phase",fontsize=10 )
plt.savefig("metanol_octanol/X2B_vap.pdf",bbox_inches='tight') 

plt.legend()

#%%
antoine=np.array([metoh_antoine,octanol_antoine])

p=CPAParameters.from_records(
    cubic=[c_methanol_3b,c_octanol_3b],
    assoc=[a_methanol_3b,a_octanol_3b])

p.set_cubic_binary(0,1,0.0,-0.025)
p.set_assoc_binary(0,1,"ecr")
eos_3b=EquationOfState.cpa(p)
P=101.32e3
xorv,torv=metoh_otctanol["orv"]
xbol,tbol=metoh_otctanol["bol"]

exp_data=[xorv,torv,xbol,tbol]

VAR,LIQUID,VAPOR=linspace_bubble_t(eos_3b,P,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   plot_dir="plots/metanol_octanol",
   title="(1) Methanol and (2) Octanol 3B ECR ",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
   text=f"{P/1e3}kPa",
   save_fig=True,

   factor=1.0,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=300,y_sup=480,
   
   exp_data=exp_data)

#%%

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]


#%%
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-MeOH","+MeOH","-OcOH","+OcOH"]
linestyle=['-','--','-.',':']

#%%

# liq_strength_metoh=np.zeros_like(LIQUID,dtype=object)
# liq_strength_oct=np.zeros_like(LIQUID,dtype=object)

# for (i,s) in enumerate(LIQUID):

#     liq_strength_metoh[i]=LIQUID[i].association_strength()[0]
#     liq_strength_oct[i]=LIQUID[i].association_strength()[2]


#%%

# plt.plot(xL,[k[1] for k in liq_strength_metoh],label="MetOH-MetOH",linestyle='-')
# # for (i,k) in enumerate(liq_strength_metoh):
# #     plt.scatter(xL[i],k[1])
# #     plt.scatter(xL[i],k[3])
# plt.plot(xL,[k[3] for k in liq_strength_metoh],label="MetOH-OctOH",linestyle='-')


# plt.plot(xL,[k[3] for k in liq_strength_oct],label="OctOH-OctOH",linestyle='-')

# plt.legend()

#%%

plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):

    plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Methanol and (2) Octanol ECR 3B Liquid Phase",fontsize=10)
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$x_1$")
plt.legend()

plt.savefig("plots/metanol_octanol/X3B_liq.pdf",bbox_inches='tight') 
#%%

# for (i,s) in enumerate(sites):

#     plt.plot(xL,XL[:,i],label=s)


# plt.title("Methanol(2+ and 1-) and Octanol(1+ and 2-)  Liquid Phase")
# plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlabel(r"$x_1$") 
# plt.legend()
# plt.savefig("xassoc_plot/Methanol(2+ and 1-)_and_Octanol(1+ and 2-)_LiquidPhase.png",format='png')


#%%
plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):

    plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])

plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$y_1$") 

plt.title("(1) Methanol and (2) Octanol ECR 3B Vapor Phase",fontsize=10) 
plt.legend()
plt.savefig("plots/metanol_octanol/X3B_vap.pdf",bbox_inches='tight') 

#%%
#%%
# for (i,s) in enumerate(sites):

#     plt.plot(xV,XV[:,i],label=s)

# plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlabel(r"$x_1$") 

# plt.title("Methanol(2+ and 1-) and Octanol(1+ and 2-) Vapor Phase") 
# plt.legend()
# plt.savefig("xassoc_plot/Methanol(2+ and 1-)_and_Octanol(1+ and 2-)_VaporPhase.png",format='png')


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
bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1.0,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=360,
   y_sup=420,
   text=f"{P/1e3}kPa",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
    save_fig=yes_or_no,

   title="Propanoic Acid 1A(1) and Heptane(2)",
   exp_data=exp_data)

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
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=0.0,
   y_sup=30,
   text=f"{T}K",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   save_fig=True,
   plot_dir="plots/metanol_acetic",

   title="(1) Methanol 2B and (2) AcOH 1A ECR",
   exp_data=exp_data)


#%%

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-MetOH","+MetOH","+- AcOH"]
linestyle=['-','--','-.']

#%%

plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Methanol 2B and (2) AcOH 1A Liquid Phase",fontsize=10)
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$x_1$")
plt.legend()

plt.savefig("plots/metanol_acetic/X_liq.pdf",bbox_inches='tight') 


# %%
plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Methanol 2B and (2) AcOH 1A Vapor Phase",fontsize=10)
plt.ylabel("Fraction of Non-Bonded Sites")

plt.xlabel(r"$y_1$")
plt.legend()
plt.savefig("plots/metanol_acetic/X_vap.pdf",bbox_inches='tight') 
#%%

p=CPAParameters.from_records(
    cubic=[c_w,c_acoh],
    assoc=[a_w,a_acoh])

p.set_cubic_binary(0,1,0.0,-0.222)
p.set_assoc_binary(0,1,"ecr")
eos=EquationOfState.cpa(p)
# print(pWATER_ACETIC.as_string())

# p,y,vx=vle_diagram(T,peq,factor=1e3)

antoine=np.array([water_antoine,acoh_antoine])
xorv,porv=water_acoh["orv"]
xbol,pbol=water_acoh["bol"]
exp_data=[xorv,porv,xbol,pbol]
T=313.15



VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1e5,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=0.0,
   y_sup=0.1,
   text=f"{T}K",
   y_label="P/bar",
   x_label=r"$x_1,y_1$",
   save_fig=True,
   plot_dir="plots/water_acetic",
   title="(1) Water 4C and (2) AcOH 1A ECR",
   exp_data=exp_data)

#%%

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-H2O","+H2O","+- AcOH"]
linestyle=['-','--','-.']

#%%

plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Water 4C and (2) AcOH 1A  Liquid Phase",fontsize=10)
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$x_1$")
plt.legend()

plt.savefig("plots/water_acetic/X_liq.pdf",bbox_inches='tight') 


# %%
plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):
    
    
    plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Water 4C and (2) AcOH 1A Vapor Phase",fontsize=10)
plt.ylabel("Fraction of Non-Bonded Sites")

plt.xlabel(r"$y_1$")
plt.legend()
plt.savefig("plots/water_acetic/X_vap.pdf",bbox_inches='tight') 

#%%

p=CPAParameters.from_records(
    cubic=[c_w,c_mea],
    assoc=[a_w,a_mea])

p.set_cubic_binary(0,1,0.0,-0.165)
# p.set_assoc_binary(0,1,"ecr")
eos=EquationOfState.cpa(p)


antoine=np.array([water_antoine,mea_antoine])
P=101.33e3

xorv,porv=water_mea_1atm["orv"]
xbol,pbol=water_mea_1atm["bol"]
exp_data=[xorv,porv,xbol,pbol]



VAR,LIQUID,VAPOR=linspace_bubble_t(eos,P,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1.0,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=340,
   y_sup=480,
   text=f"{P/1e3}kPa",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
   save_fig=True,
   plot_dir="plots/water_mea",
   title="(1) Water and (2) MEA 4C CR1 ",
   exp_data=exp_data)

#%%

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-H2O","+H2O","-MEA","+MEA"]
linestyle=['-','--','-.',':']

#%%

plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Water and (2) MEA Liquid Phase",fontsize=10)


plt.ylabel("Fraction of Non-Bonded Sites",fontsize=10)
plt.xlabel(r"$x_1$")
plt.legend()

plt.savefig("plots/water_mea/X_liq.pdf",bbox_inches='tight')


#%%
plt.figure(figsize=(xsize,ysize)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Water and (2) MEA Vapor Phase")
plt.ylabel("Fraction of Non-Bonded Sites",fontsize=10)
plt.xlabel(r"$y_1$")
plt.legend()

plt.savefig("plots/water_mea/X_vap.pdf",bbox_inches='tight')

#%%
p=CPAParameters.from_records(
    cubic=[c_ethanol_3b,c_w],
    assoc=[a_ethanol_3b,a_w])
# 
# p=CPAParameters.from_records(
#     cubic=[c_w,c_ethanol_3b],
#     assoc=[a_w,a_ethanol_3b])
p.set_cubic_binary(0,1,0.0,0.0)
p.set_assoc_binary(0,1,"ecr")
eos=EquationOfState.cpa(p)
ethanol_antoine=np.array([1.24677,1598.673,-46.424])


antoine=np.array([ethanol_antoine,water_antoine])
# P=101.33e3
T=298.14

x,pexp=water_ethanol["orv"]
y,_=water_ethanol["bol"]

# xbol,pbol=water_mea_1atm["bol"]
exp_data=[x,pexp,y,pexp]




VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)


#%%
bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=2.0,
   y_sup=9.0,
   text=f"{T}K",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   exp_data=exp_data,
   save_fig=True,
   plot_dir="plots/water_ethanol",
   title="(1) Ethanol 3B and (2) Water 4C ECR",
   )

#%%

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-EtOH","+EtOH","-H2O","+H2O"]
linestyle=['-','--','-.',':']

#%%
ethanol=np.zeros_like(LIQUID,dtype=object)
water=np.zeros_like(LIQUID,dtype=object)

for (i,s) in enumerate(LIQUID):

    water[i]=LIQUID[i].association_strength()[2]
    ethanol[i]=LIQUID[i].association_strength()[0]
#%%
plt.scatter(xL,VAR/1e5)
#%%

# plt.figure(figsize=(5,4.0)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Ethanol and (2) Water Liquid Phase")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$x_1$")
plt.legend()

plt.savefig("plots/water_ethanol/X_liq.pdf",bbox_inches='tight')

#%%

plt.ylabel("$\Delta_{asc}$")
plt.xlabel(r"$x_{agua}$")
plt.plot(xL,[k[2] for k in water],label="Agua-Agua",linestyle='-')
plt.plot(xL,[k[0] for k in water],label="Agua-Ethanol",linestyle='-')
plt.plot(xL,[k[0] for k in ethanol],label="Ethanol-Ethanol",linestyle='-')

plt.legend()

# plt.plot(xL,[k[3] for k in liq_strength_oct],label="OctOH-OctOH",linestyle='-')
plt.savefig("plots/water_ethanol/Delta_assoc_liq.pdf",bbox_inches='tight')

#%%
plt.figure(figsize=(5,4.0)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Ethanol and (2) Water Vapor Phase")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$y_1$")
plt.legend()

plt.savefig("plots/water_ethanol/X_vap.pdf",bbox_inches='tight')
# %%
#%%

p=CPAParameters.from_records(
    cubic=[c_w,c_mea],
    assoc=[a_w,a_mea])

p.set_cubic_binary(0,1,0.0,-0.165)
# p.set_assoc_binary(0,1,"ecr")
eos=EquationOfState.cpa(p)


antoine=np.array([water_antoine,mea_antoine])
T=400.15

xorv,porv=water_mea_1atm["orv"]
xbol,pbol=water_mea_1atm["bol"]
# exp_data=[xorv,porv,xbol,pbol]



VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=ysize,
   x_figsize=xsize,
#    y_inf=300,
#    y_sup=500,
   text=f"{T}K",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   save_fig=True,
   plot_dir="plots/water_mea",
   title="(1) Water and (2) MEA 4C CR1 ",)

#%%

xL=[LIQUID[i].composition()[0] for i in range(100)]
xV=[VAPOR[i].composition()[0] for i in range(100)]
XL= get_non_bondend_sites_from_states(LIQUID)
XV= get_non_bondend_sites_from_states(VAPOR)
sites=["-H2O","+H2O","-MEA","+MEA"]
linestyle=['-','--','-.',':']

#%%

plt.figure(figsize=(5,4.0)) 

for (i,s) in enumerate(sites):
    
    plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


plt.title("(1) Water and (2) MEA Liquid Phase")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlabel(r"$x_1$")
plt.legend()

# plt.savefig("plots/water_mea/X_liq.pdf",bbox_inches='tight')

#%%

mea=np.zeros_like(LIQUID,dtype=object)
water=np.zeros_like(LIQUID,dtype=object)

for (i,s) in enumerate(LIQUID):

    water[i]=LIQUID[i].association_strength()[0]
    mea[i]=LIQUID[i].association_strength()[2]


#%%

plt.ylabel("$\Delta_{asc}$")
plt.xlabel(r"$x_{agua}$")
plt.plot(xL,[k[0] for k in water],label="Agua-Agua",linestyle='-')
plt.plot(xL,[k[2] for k in water],label="Agua-MEA",linestyle='-')
plt.plot(xL,[k[2] for k in mea],label="MEA-MEA",linestyle='-')

plt.legend()
# plt.plot(xL,[k[3] for k in liq_strength_oct],label="OctOH-OctOH",linestyle='-')
plt.savefig("plots/water_mea/Delta_assoc_liq.pdf",bbox_inches='tight')


#%%


# %%
