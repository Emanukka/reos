#%% 

from reos.eos import EquationOfState
from reos.state import State
from reos.cpa import CPAParameters

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

from phase_diagram.antoine import utilis
from phase_diagram.exp_data import *
from phase_diagram.vle_functions import *

xsize = 3.15
ysize = 3.15

save_fig = True
plot_dir = "plots/phase_diagrams"

# yes_or_no=False
# plt.rcParams.update({
#     "font.serif": ["Computer Modern"], 
#     "axes.labelsize": 10,
#     "font.size": 10,
#     "legend.fontsize": 10,
#     "xtick.labelsize": 10,
#     "ytick.labelsize": 10,
#     "figure.figsize": (xsize, ysize),  

# })

def run_bubble_t(name1, name2, P,
                 ppath = "../../parameters/cpa/kontogeorgis2010.json", 
                 bpath = "../../parameters/cpa/kontogeorgis2010_binary.json"):
    
    p = CPAParameters.from_json(
        [name1, name2],
        ppath,
        bpath)

    eos = EquationOfState.scpa(p)

    antoine = np.array([utilis[name1], utilis[name2]])


    T, LIQUID, VAPOR=linspace_bubble_t(eos, P, antoine, N = 100)
    
    return T, LIQUID, VAPOR

def run_bubble_p(name1, name2, T,
                 ppath = "../../parameters/cpa/kontogeorgis2010.json", 
                 bpath = "../../parameters/cpa/kontogeorgis2010_binary.json"):

    p = CPAParameters.from_json(
        [name1, name2],
        ppath,
        bpath)

    eos = EquationOfState.scpa(p)

    antoine = np.array([utilis[name1], utilis[name2]])


    P, LIQUID, VAPOR=linspace_bubble_p(eos, T, antoine, N = 100)
    
    return P, LIQUID, VAPOR


#%% 

name1 = "acetic acid"
name2 = "n-octane"

xorv,porv = acoh_octane["orv"]
xbol,pbol = acoh_octane["bol"]
exp_data = [xorv,porv,xbol,pbol]

T=343.2

# PRES, LIQUID, VAPOR = linspace_bubble_p(eos,T, antoine,N=100)
PRES, LIQUID, VAPOR = run_bubble_p(name1, name2, T)
PRES

#%%
bubble_diagram(
   PRES,
   LIQUID,
   VAPOR,
   plot_dir = plot_dir,
   factor = 1e3,
   y_figsize = ysize,
   x_figsize = xsize,
   y_inf=15,
   y_sup=30,
   text=f"{T}K",
   title="(1) Acetic Acid 1A and (2) Octane",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   save_fig = save_fig,
   

   exp_data=exp_data)

#%%

name1 = "methanol"
name2 = "1-octanol"

xorv, torv = metoh_otctanol["orv"]
xbol, tbol = metoh_otctanol["bol"]
exp_data = [xorv,torv,xbol,tbol]

P = 101.32e3

VAR,LIQUID,VAPOR = run_bubble_t(name1,name2, P, bpath=None)

#%%

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   plot_dir = plot_dir,
   title="(1) Methanol 2B and (2) Octanol 2B CR1 ",
   y_label="T/K",
   x_label=r"$x_1,y_1$",
   text=f"{P/1e3}kPa",

   factor=1.0,
   y_figsize=ysize,
   x_figsize=xsize,
   y_inf=300,y_sup=480,
   save_fig = save_fig,
   exp_data=exp_data)

#%%

name1 = "propionic acid"
name2 = "n-heptane"

propanoic_hep_teb
xorv,torv = propanoic_hep_teb["orv"]
xbol,tbol = propanoic_hep_teb["bol"]
exp_data=[xorv,torv,xbol,tbol]

P=101.33e3
exp_data = [xorv,torv,xbol,tbol]

VAR,LIQUID,VAPOR = run_bubble_t(name1,name2, P)

# VAR,LIQUID,VAPOR=linspace_bubble_t(PROPANOIC_HEPTANE,P,antoine,N=100)
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
   save_fig = save_fig,
   plot_dir = plot_dir,
   title="Propanoic Acid 1A(1) and Heptane(2)",
   exp_data=exp_data)


#%%

# pMETHANOL_ACETIC=CPAParameters.from_records(
#     cubic=[c_methanol_2b,c_acoh],
#     assoc=[a_methanol_2b,a_acoh])

# pMETHANOL_ACETIC.set_cubic_binary(0,1,0.0,-0.04)
# pMETHANOL_ACETIC.set_assoc_binary(0,1,"ecr")
# METHANOL_ACETIC=EquationOfState.cpa(pMETHANOL_ACETIC)


# antoine=np.array([metoh_antoine,acoh_antoine])
# T=308.15

# xorv,porv=metoh_acoh["orv"]
# xbol,pbol=metoh_acoh["bol"]
# exp_data=[xorv,porv,xbol,pbol]

# PRES,LIQUID,VAPOR=linspace_bubble_p(METHANOL_ACETIC,T,antoine,N=100)

# bubble_diagram(
#    PRES,
#    LIQUID,
#    VAPOR,
#    factor=1e3,
#    y_figsize=ysize,
#    x_figsize=xsize,
#    y_inf=0.0,
#    y_sup=30,
#    text=f"{T}K",
#    y_label="P/kPa",
#    x_label=r"$x_1,y_1$",
#    save_fig=True,
#    plot_dir="plots/metanol_acetic",

#    title="(1) Methanol 2B and (2) AcOH 1A ECR",
#    exp_data=exp_data)


#%%

name1 = "acetic acid"
name2 = "water"

xorv,porv=water_acoh["orv"]
xbol,pbol=water_acoh["bol"]
exp_data=[xorv,porv,xbol,pbol]

T = 313.15



VAR,LIQUID,VAPOR= run_bubble_p(name2, name1, T, 
                               ppath = "../../parameters/cpa/kontogeorgis2006.json",
                               bpath = "../../parameters/cpa/kontogeorgis2006_binary.json") 

#%%
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
save_fig = save_fig,
plot_dir = plot_dir,
title="(1) Water 4C and (2) AcOH 1A ECR",
exp_data=exp_data)

# #%%

# xL=[LIQUID[i].composition()[0] for i in range(100)]
# xV=[VAPOR[i].composition()[0] for i in range(100)]
# XL= get_non_bondend_sites_from_states(LIQUID)
# XV= get_non_bondend_sites_from_states(VAPOR)
# sites=["-H2O","+H2O","+- AcOH"]
# linestyle=['-','--','-.']


# plt.figure(figsize=(xsize,ysize)) 

# for (i,s) in enumerate(sites):
    
#     plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


# plt.title("(1) Water 4C and (2) AcOH 1A  Liquid Phase",fontsize=10)
# plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlabel(r"$x_1$")
# plt.legend()

# plt.savefig("plots/water_acetic/X_liq.pdf",bbox_inches='tight') 


# # %%
# plt.figure(figsize=(xsize,ysize)) 

# for (i,s) in enumerate(sites):
    
    
#     plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])


# plt.title("(1) Water 4C and (2) AcOH 1A Vapor Phase",fontsize=10)
# plt.ylabel("Fraction of Non-Bonded Sites")

# plt.xlabel(r"$y_1$")
# plt.legend()
# plt.savefig("plots/water_acetic/X_vap.pdf",bbox_inches='tight') 




# p=CPAParameters.from_records(
#     cubic=[c_ethanol_3b,c_w],
#     assoc=[a_ethanol_3b,a_w])
# # 
# # p=CPAParameters.from_records(
# #     cubic=[c_w,c_ethanol_3b],
# #     assoc=[a_w,a_ethanol_3b])
# p.set_cubic_binary(0,1,0.0,0.0)
# p.set_assoc_binary(0,1,"ecr")
# eos=EquationOfState.cpa(p)
# ethanol_antoine=np.array([1.24677,1598.673,-46.424])


# antoine=np.array([ethanol_antoine,water_antoine])
# # P=101.33e3
# T=298.14

# x,pexp=water_ethanol["orv"]
# y,_=water_ethanol["bol"]

# # xbol,pbol=water_mea_1atm["bol"]
# exp_data=[x,pexp,y,pexp]




# VAR,LIQUID,VAPOR = linspace_bubble_p(eos,T,antoine,N=100)


# #%%
# bubble_diagram(
#    VAR,
#    LIQUID,
#    VAPOR,
#    factor=1e3,
#    y_figsize=ysize,
#    x_figsize=xsize,
#    y_inf=2.0,
#    y_sup=9.0,
#    text=f"{T}K",
#    y_label="P/kPa",
#    x_label=r"$x_1,y_1$",
#    exp_data=exp_data,
#    save_fig=True,
#    plot_dir="plots/water_ethanol",
#    title="(1) Ethanol 3B and (2) Water 4C ECR",
#    )

#%%

# xL=[LIQUID[i].composition()[0] for i in range(100)]
# xV=[VAPOR[i].composition()[0] for i in range(100)]
# XL= get_non_bondend_sites_from_states(LIQUID)
# XV= get_non_bondend_sites_from_states(VAPOR)
# sites=["-EtOH","+EtOH","-H2O","+H2O"]
# linestyle=['-','--','-.',':']

# #%%
# ethanol=np.zeros_like(LIQUID,dtype=object)
# water=np.zeros_like(LIQUID,dtype=object)

# for (i,s) in enumerate(LIQUID):

#     water[i]=LIQUID[i].association_strength()[2]
#     ethanol[i]=LIQUID[i].association_strength()[0]
# plt.scatter(xL,VAR/1e5)
# #%%

# # plt.figure(figsize=(5,4.0)) 

# for (i,s) in enumerate(sites):
    
#     plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


# plt.title("(1) Ethanol and (2) Water Liquid Phase")
# plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlabel(r"$x_1$")
# plt.legend()

# plt.savefig("plots/water_ethanol/X_liq.pdf",bbox_inches='tight')

# #%%

# plt.ylabel("$\Delta_{asc}$")
# plt.xlabel(r"$x_{agua}$")
# plt.plot(xL,[k[2] for k in water],label="Agua-Agua",linestyle='-')
# plt.plot(xL,[k[0] for k in water],label="Agua-Ethanol",linestyle='-')
# plt.plot(xL,[k[0] for k in ethanol],label="Ethanol-Ethanol",linestyle='-')

# plt.legend()

# # plt.plot(xL,[k[3] for k in liq_strength_oct],label="OctOH-OctOH",linestyle='-')
# plt.savefig("plots/water_ethanol/Delta_assoc_liq.pdf",bbox_inches='tight')

# #%%
# plt.figure(figsize=(5,4.0)) 

# for (i,s) in enumerate(sites):
    
#     plt.plot(xV,XV[:,i],label=s,color="black",linestyle=linestyle[i])


# plt.title("(1) Ethanol and (2) Water Vapor Phase")
# plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlabel(r"$y_1$")
# plt.legend()

# plt.savefig("plots/water_ethanol/X_vap.pdf",bbox_inches='tight')
# %%

# p=CPAParameters.from_records(
#     cubic=[c_w,c_mea],
#     assoc=[a_w,a_mea])

# p.set_cubic_binary(0,1,0.0,-0.165)
# # p.set_assoc_binary(0,1,"ecr")
# eos=EquationOfState.cpa(p)


# antoine=np.array([water_antoine,mea_antoine])
# T=400.15

# xorv,porv=water_mea_1atm["orv"]
# xbol,pbol=water_mea_1atm["bol"]
# # exp_data=[xorv,porv,xbol,pbol]



# VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)

# bubble_diagram(
#    VAR,
#    LIQUID,
#    VAPOR,
#    factor=1e3,
#    y_figsize=ysize,
#    x_figsize=xsize,
# #    y_inf=300,
# #    y_sup=500,
#    text=f"{T}K",
#    y_label="P/kPa",
#    x_label=r"$x_1,y_1$",
#    save_fig=True,
#    plot_dir="plots/water_mea",
#    title="(1) Water and (2) MEA 4C CR1 ",)

# #%%

# xL=[LIQUID[i].composition()[0] for i in range(100)]
# xV=[VAPOR[i].composition()[0] for i in range(100)]
# XL= get_non_bondend_sites_from_states(LIQUID)
# XV= get_non_bondend_sites_from_states(VAPOR)
# sites=["-H2O","+H2O","-MEA","+MEA"]
# linestyle=['-','--','-.',':']

# #%%

# plt.figure(figsize=(5,4.0)) 

# for (i,s) in enumerate(sites):
    
#     plt.plot(xL,XL[:,i],label=s,color="black",linestyle=linestyle[i])


# plt.title("(1) Water and (2) MEA Liquid Phase")
# plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlabel(r"$x_1$")
# plt.legend()

# # plt.savefig("plots/water_mea/X_liq.pdf",bbox_inches='tight')

# #%%

# mea=np.zeros_like(LIQUID,dtype=object)
# water=np.zeros_like(LIQUID,dtype=object)

# for (i,s) in enumerate(LIQUID):

#     water[i]=LIQUID[i].association_strength()[0]
#     mea[i]=LIQUID[i].association_strength()[2]


# #%%

# plt.ylabel("$\Delta_{asc}$")
# plt.xlabel(r"$x_{agua}$")
# plt.plot(xL,[k[0] for k in water],label="Agua-Agua",linestyle='-')
# plt.plot(xL,[k[2] for k in water],label="Agua-MEA",linestyle='-')
# plt.plot(xL,[k[2] for k in mea],label="MEA-MEA",linestyle='-')

# plt.legend()
# # plt.plot(xL,[k[3] for k in liq_strength_oct],label="OctOH-OctOH",linestyle='-')
# plt.savefig("plots/water_mea/Delta_assoc_liq.pdf",bbox_inches='tight')

#%%
