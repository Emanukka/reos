#%%


from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from auxiliary_functions.parameters import *
# from auxiliary_functions.vle_functions import *
# from auxiliary_functions.data import *
from auxiliary_functions.wsat_data import *
from auxiliary_functions.water_sat import *
from auxiliary_functions.association_functions import *
import os
yes_or_no=True
plt.rcParams.update({
    "text.usetex": True,               # Usa LaTeX
    "font.family": "serif",            # Usa fonte serifada (como em artigos)
    "font.serif": ["Computer Modern"], # Fonte do LaTeX padrão
    "axes.labelsize": 12,
    "font.size": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})


#%%



#%%
pWATER_CO2=CPAParameters.from_records(
    cubic=[c_w,c_co2],
    assoc=[a_w,a_co2])

pure_water=CPAParameters.from_records(
   [c_w],[a_w]

)

eos_pure_water=EquationOfState.cpa(pure_water)

pWATER_CO2.set_cubic_binary(0,1,kij_a=0.000877,kij_b=-0.15508)

pWATER_CO2.set_assoc_binary(0,1,"mcr1",beta=0.1836)


eos=EquationOfState.cpa(pWATER_CO2)
#%%calc
T=np.array([298.15,308.15,323.15])
# T=np.array([298.15])

ydg=np.array([0.0,0.5])
pResultBAR=np.zeros_like(T,dtype=object)
yResultPPM=np.zeros_like(T,dtype=object)
stateResult=np.zeros_like(T,dtype=object)

for (i,temp) in enumerate(T):

  try:
     
    p,y,s=linspace_wsat(eos,eos_pure_water,temp,ydg)
    pResultBAR[i]=p*1
    yResultPPM[i]=y*1
    stateResult[i]=s

  except Exception as e:
    print(e)
    pass
    continue
#%%
markers = ['o', 's', '^', 'D', 'P', 'X'] 
lines= ['-','--',':']
plt.figure(figsize=(5, 5))

for (i,t) in enumerate(T):


  pe,ye=wsat_data[t]

  # plt.plot(pResultBAR[i],yResultPPM[i],)
  plt.plot(pResultBAR[i],
           yResultPPM[i],
           color="Black",
           linestyle=lines[i],
           label=f"{t}K")
  
  plt.scatter(pe,
              ye*1e6,
              marker=markers[i],
              facecolors='none', 
              edgecolors='black')

  plt.ylim(0,20*1e3)
  plt.ylabel("Water Mole Fraction (ppm)")
  plt.xlabel("P/bar")

  # plt.xlim(20,100)

plt.title(f"H2O 4C and CO2 1ea")
plt.legend()

os.makedirs("water_sat_plot", exist_ok=True)
filepath = os.path.join("water_sat_plot", "H2O_CO2")
plt.savefig(filepath)
# plt.text(500,17500,"H2O 4C, CO2 1ea")



#%%

vX=np.zeros_like(T,dtype=object)

for i in range(len(T)):

  vX[i]=get_non_bondend_sites_from_states(stateResult[i])

#%%

i=0
AGUA_NEGATIVO=np.zeros_like(stateResult[i],dtype=object)
AGUA_POSITIVO=np.zeros_like(stateResult[i],dtype=object)
CO2_POSITIVO=np.zeros_like(stateResult[i],dtype=object)

for j,state in enumerate(stateResult[i]):

    # idx 0 - agua negativo
    # idx 1 - agua positivo 
    # idx 2 - co2 positivo
    AGUA_NEGATIVO[j],AGUA_POSITIVO[j],CO2_POSITIVO[j]=vX[i][j]
    
    agua_negativo=AGUA_NEGATIVO[j]
    agua_positivo=AGUA_POSITIVO[j]
    co2_positivo=AGUA_POSITIVO[j]

    yw,yco2=state.composition()

#%%

plt.plot(pResultBAR[i],AGUA_NEGATIVO,color="black",linestyle=lines[0],label="H2O -")
plt.plot(pResultBAR[i],AGUA_POSITIVO,color="black",linestyle=lines[1],label="H2O +")
plt.plot(pResultBAR[i],CO2_POSITIVO, color="black",linestyle=lines[2],label="CO2 +")
plt.xlabel("P/bar")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlim(0,200)
plt.text(150,0.99,"298.15K",fontsize=12)
plt.legend()

os.makedirs("water_sat_plot", exist_ok=True)
filepath = os.path.join("water_sat_plot", "H2O_CO2_non_bonded_sites_298")
plt.savefig(filepath)

#%%
#%%



pWATER_CO2=CPAParameters.from_records(
    cubic=[c_w,c_co2],
    assoc=[a_w,a_co2])



pWATER_CO2.set_cubic_binary(0,1,kij_a=0.0,kij_b=0.1141)

pWATER_CO2.set_assoc_binary(0,1,"exp",eps=142e2,beta=0.0162)


eos=EquationOfState.cpa(pWATER_CO2)

T=np.array([298.15])

ydg=np.array([0.0,0.5])
pResultBAR=np.zeros_like(T,dtype=object)
yResultPPM=np.zeros_like(T,dtype=object)
stateResult=np.zeros_like(T,dtype=object)

for (i,temp) in enumerate(T):

  try:
     
    p,y,s=linspace_wsat(eos,temp,ydg)
    pResultBAR[i]=p*1
    yResultPPM[i]=y*1
    stateResult[i]=s

  except Exception as e:
    print(e)
    pass
    continue

#%%

# i=0


# pe,ye=wsat_data[T[i]]
# plt.plot(pResultBAR[i],yResultPPM[i],)
# plt.scatter(pe,ye*1e6)
# plt.ylim(0,20*1e3)
# # plt.xlim(20,100)
# plt.title(f"T={T[i]}K,")

i=1
pe,ye=wsat_data[T[i]]
plt.plot(pResultBAR[i],yResultPPM[i],)
plt.scatter(pe,ye*1e6)
plt.ylim(0,20*1e3)
# plt.xlim(20,100)
plt.title(f"T={T[i]}K,")