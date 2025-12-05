#%% main

from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt

from auxiliary_functions.parameters import *
from auxiliary_functions.vle_functions import *
# from python_examples.auxiliary_functions.exp_data import *
from auxiliary_functions.exp_data import *
from auxiliary_functions.association_functions import get_non_bondend_sites_from_states
#%%
c_acoh=CubicRecord(
    a0= 0.70592,
    b=0.0478e-3,
    c1= 0.8808,
    tc=594.8,
)

a_acoh=AssociationRecord.associative(
    eps=188.23e2,
    beta=0.1408,
    b=0.0478e-3,
    na=1,
    nb=1,
    nc=0)
p=CPAParameters.from_records(
    cubic=[c_w,c_acoh],
    assoc=[a_w,a_acoh])

p.set_cubic_binary(0,1,0.0,0)
# p.set_assoc_binary(0,1,"ecr")
eos=EquationOfState.cpa(p)
asc=eos.get_association()



#%%
ACID_DN=2
ACID_AC=3
WATER_DN=0
WATER_AC=1


#%%

# print(pWATER_ACETIC.as_string())

# p,y,vx=vle_diagram(T,peq,factor=1e3)

antoine=np.array([water_antoine,acoh_antoine])
xorv,porv=water_acoh["orv"]
xbol,pbol=water_acoh["bol"]
exp_data=[xorv,porv,xbol,pbol]

T=313.15
#%%
N=50
vz=np.linspace(0.00001,0.9999,N)
PRES=np.zeros_like(vz)
LIQUID=np.zeros_like(vz,dtype=object)
# print(LIQUID)
VAPOR=np.zeros_like(vz,dtype=object)
# print(VAPOR)
  
xbol_calc=np.zeros_like(vz)
xorv_calc=np.zeros_like(vz)

for (i,z1) in enumerate(vz):
  try:
      print(i)
      z=np.array([z1,1.0-z1])
      liquid_phase,vapor_phase=bubble_p(eos,T,z,antoine)
      PRES[i]=liquid_phase.pressure()
    #   LIQUID[i]=liquid_phase 
    #   VAPOR[i]=vapor_phase
      xorv_calc[i]=liquid_phase.composition()[0]

      xbol_calc[i]=vapor_phase.composition()[0]
  except Exception as e:
      print(e)

#%%
T=313.15

N=50
# LIQUID=np.zeros_like(vz,dtype=object)
# print(LIQUID)
# VAPOR=np.zeros_like(vz,dtype=object)
# print(VAPOR)
  

xorv[0]=1e-16
xorv[-1]=0.99999999999
PRES=np.zeros_like(xorv)
xbol_calc=np.zeros_like(xorv)
xorv_calc=np.zeros_like(xorv)

for (i,z1) in enumerate(xorv):

  try:
      print(i)
      z=np.array([z1,1.0-z1])
      liquid_phase,vapor_phase=bubble_p(eos,T,z,antoine)
      PRES[i]=liquid_phase.pressure()
    #   LIQUID[i]=liquid_phase 
    #   VAPOR[i]=vapor_phase
      xorv_calc[i]=liquid_phase.composition()[0]

      xbol_calc[i]=vapor_phase.composition()[0]
  except Exception as e:
      print(e)
# for (i,z1) in enumerate(xbol):
#   try:
#       print(i)
#       z=np.array([z1,1.0-z1])
#       liquid_phase,vapor_phase=orv_p(eos,T,z,antoine)
#       PRES[i]=liquid_phase.pressure()
#     #   LIQUID[i]=liquid_phase 
#     #   VAPOR[i]=vapor_phase
#       xorv_calc[i]=liquid_phase.composition()[0]

#       xbol_calc[i]=vapor_phase.composition()[0]
#   except Exception as e:
#       print(e)


#%%
PRES=np.zeros_like(xorv)
# xbol_calc=np.zeros_like(xorv)
xorv_calc=np.zeros_like(xorv)
eps_12=17739.0
eps_03=17739.0
beta_03=0.09870845961719796
beta_12=0.09870845961719796

eps_and_beta=np.array([eps_03,eps_12,beta_03,beta_12])
#%%
def fobj_porv(eps_and_beta):
    p=CPAParameters.from_records(
    cubic=[c_w,c_acoh],
    assoc=[a_w,a_acoh])
    eps0 =eps_and_beta[0]
    eps1 =eps_and_beta[1]
    beta0=eps_and_beta[2]
    beta1=eps_and_beta[3]
    p.set_cubic_binary(0,1,0.0,-0.2)
    p.change_sites_p(0,3,eps0,beta0)
    p.change_sites_p(1,2,eps1,beta1)

    # p.set_assoc_binary(0,1,"ecr")
    eos=EquationOfState.cpa(p)
    for (i,z1) in enumerate(xorv):
        z=np.array([z1,1.0-z1])
        liquid_phase,vapor_phase=bubble_p(eos,T,z,antoine)
        PRES[i]=liquid_phase.pressure()/1e5
        #  LIQUID[i]=liquid_phase 
        #  VAPOR[i]=vapor_phase
        # xorv_calc[i]=liquid_phase.composition()[0]
        # xbol_calc[i]=vapor_phase.composition()[0]
        # print(liquid_phase.pressure()/1e5)
    err_porv=np.linalg.norm(PRES-porv)
    # err_xbol=np.linalg.norm(xbol_calc-xbol)
    print(eps_and_beta)
    print(err_porv)
    return err_porv 


#%%
import scipy.optimize as opt
#%%
opt.minimize(fobj_porv,eps_and_beta)

#%%
p=CPAParameters.from_records(
cubic=[c_w,c_acoh],
assoc=[a_w,a_acoh])
eps0 =eps_and_beta[0]
eps1 =eps_and_beta[1]
beta0=eps_and_beta[2]
beta1=eps_and_beta[3]
p.set_cubic_binary(0,1,0.0,0)
# p.change_sites_p(0,3,9,9)
p.change_sites_p(1,2,9,9)
print(p.as_string())
#%%
plt.scatter(xorv_calc,PRES/1e5,marker="s")
# plt.scatter(xbol_calc,PRES/1e5,marker="s")
# plt.scatter(xorv,porv)
# plt.scatter(xbol,pbol)

# VAR,LIQUID,VAPOR=linspace_bubble_p(eos,T,antoine,N=100)

# bubble_diagram(
#    VAR,
#    LIQUID,
#    VAPOR,
#    factor=1e5,
#    # y_figsize=ysize,
#    # x_figsize=xsize,
#    y_inf=0.0,
#    y_sup=0.1,
#    text=f"{T}K",
#    y_label="P/bar",
#    x_label=r"$x_1,y_1$",
#    save_fig=False,
#    plot_dir="plots/water_acetic",
#    title="(1) Water 4C and (2) AcOH 1A ECR",
#    exp_data=exp_data)
#%%

