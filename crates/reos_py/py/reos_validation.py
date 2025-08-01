
#%%


from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from aux.parameters import *
from aux.vle_functions import *
from aux.data import *



#%%
pACOH=CPAParameters.from_records(
    cubic=[c_acoh],
    assoc=[a_acoh])


ACOH=EquationOfState.cpa(pACOH)
# print(pWATER_ACETIC.as_string())

peq=PhaseEquilibrium(ACOH)

antoine=np.array([acoh_antoine])
# p,y,vx=vle_diagram(T,peq,factor=1e3)




#%%
pACOH_OCT=CPAParameters.from_records(
    cubic=[c_acoh,c_octane],
    assoc=[a_acoh,a_octane])

pACOH_OCT.set_cubic_binary(0,1,0.0, 0.064)

ACOH_OCT=EquationOfState.cpa(pACOH_OCT)
# print(pWATER_ACETIC.as_string())

peq=PhaseEquilibrium(ACOH_OCT)

antoine=np.array([acoh_antoine,octane_antoine])
# p,y,vx=vle_diagram(T,peq,factor=1e3)

xorv,porv=acoh_octane["orv"]
xbol,pbol=acoh_octane["bol"]
exp_data=[xorv,porv,xbol,pbol]

T=343.15

_=VLE_DIAGRAM(
    ("t",T),
    peq,antoine,
    y_label="P/kPa",
    x_label="x1,y1",
    y_lim=[15,31],
    x_figsize=5,
    factor=1e3,
    exp_data=exp_data,
    title="Acetic Acid(1) and Octane",
    save_fig=True,
    N_points=100)



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
expdata=[xorv,porv,xbol,pbol]
_,_,_=VLE_DIAGRAM(
    ("t",T),
    peq,antoine,
    y_label="P/kPa",
    x_label="x1,y1",
    factor=1e3,
    exp_data=expdata,
    title="Propanoic Acid(1) and Heptane",
    y_lim=[0.0,25.0],
    save_fig=True,
    N_points=100)


#%%
p=CPAParameters.from_records(
    cubic=[c_methanol_2b,c_octanol_2b],
    assoc=[a_methanol_2b,a_octanol_2b])

p.set_cubic_binary(0,1,0.0,0.0)
p.set_assoc_binary(0,1,"cr1")
eos_2b=EquationOfState.cpa(p)
peq_2b=PhaseEquilibrium(eos_2b)



antoine=np.array([metoh_antoine,octanol_antoine])
P=101.32e3

xorv,torv=metoh_otctanol["orv"]
xbol,tbol=metoh_otctanol["bol"]

exp_data=[xorv,torv,xbol,tbol]
_,_,_=VLE_DIAGRAM(
    ("p",P),
    peq_2b,
    antoine,
    y_label="T/K",
    exp_data=exp_data,
    x_label="x1,y1",
    y_lim=[300,480],
    x_figsize=5,
    title="Methanol(1) and Octanol CR1 2B",
    save_fig=True,
    N_points=100)

#%%
p=CPAParameters.from_records(
    cubic=[c_methanol_3b,c_octanol_3b],
    assoc=[a_methanol_3b,a_octanol_3b])

p.set_cubic_binary(0,1,0.0,-0.025)
p.set_assoc_binary(0,1,"ecr")
eos_3b=EquationOfState.cpa(p)
# print(pWATER_ACETIC.as_string())
peq_3b=PhaseEquilibrium(eos_3b)
_,_,_=VLE_DIAGRAM(
    ("p",P),
    peq_3b,
    antoine,
    y_label="T/K",
    exp_data=exp_data,
    x_label="x1,y1",
    y_lim=[300,480],
    x_figsize=5,
    title="Methanol(1) and Octanol ECR 3B",
    save_fig=True,
    N_points=100)


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
P=101.32e3

_,_,_=VLE_DIAGRAM(
    ("p",P),
    peq,
    antoine,
    y_label="T/K",
    y_figsize=5,
    x_figsize=5,
    exp_data=exp_data,
    x_label="x1,y1",
    y_lim=[360,420],
    title="Propanoic Acid 1A(1) and Heptane(2)",
    save_fig=True,
    N_points=100)


#%%

pMETHANOL_ACETIC=CPAParameters.from_records(
    cubic=[c_methanol_2b,c_acoh],
    assoc=[a_methanol_2b,a_acoh])

pMETHANOL_ACETIC.set_cubic_binary(0,1,0.0,-0.04)
# pMETHANOL_ACETIC.set_assoc_binary(0,1,"ecr")
pMETHANOL_ACETIC.set_assoc_binary(0,1,"ecr")
METHANOL_ACETIC=EquationOfState.cpa(pMETHANOL_ACETIC)
# print(pWATER_ACETIC.as_string())
peq=PhaseEquilibrium(METHANOL_ACETIC)

antoine=np.array([metoh_antoine,acoh_antoine])
T=308.15

xorv,porv=metoh_acoh["orv"]
xbol,pbol=metoh_acoh["bol"]
exp_data=[xorv,porv,xbol,pbol]

_,_,_=VLE_DIAGRAM(
    ("t",T),
    peq,
    antoine,
    y_label="P/kPa",
    y_figsize=5,
    x_figsize=5,
    factor=1e3,
    exp_data=exp_data,
    x_label="x1,y1",
    y_lim=[0.0,30.0],
    title="Methanol 2B(1) and AcOH 1A(2) (ECR)",
    save_fig=True,
    N_points=100)

#%%