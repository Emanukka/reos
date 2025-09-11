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
yes_or_no=False
xsize=3.15
ysize=3.15
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

markers = ['o', 's', '^', 'D', 'P', 'X'] 
lines= ['-','--',':']
eos=EquationOfState.cpa(pWATER_CO2)
#%%calc
T=np.array([323.15])
# T=np.array([298.15])

ydg=np.array([0.0,0.5])
pResultBAR=np.zeros_like(T,dtype=object)
yResultPPM=np.zeros_like(T,dtype=object)
stateResult_rich=np.zeros_like(T,dtype=object)
stateResult_poor=np.zeros_like(T,dtype=object)

for (i,temp) in enumerate(T):

  try:
    
    #s1:developed phase
    #s2: incipient phase
    #s2 vai ser sempre fase rica em agua (co2 na agua)
    #s1 vai ser sempre fase pobre em agua (agua no co2)

    p,y,s1,s2=linspace_wsat(eos,eos_pure_water,temp,ydg,pf=400,N=100)
    pResultBAR[i]=p*1
    yResultPPM[i]=y*1
    stateResult_rich[i]=s2
    stateResult_poor[i]=s1

  except Exception as e:
    print(e)
    pass
    continue

#%%

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

  plt.ylim(0,15*1e3)
  plt.ylabel("Water Mole Fraction (ppm)")
  plt.xlabel("P/bar")

  plt.xlim(0,400)

plt.title(f"H2O 4C and CO2 1ea")
plt.legend()


plt.savefig("plots/water_co2/wsat_co2_323_vap_phase.pdf",bbox_inches='tight')
# plt.text(500,17500,"H2O 4C, CO2 1ea")



#%%

vX=np.zeros_like(T,dtype=object)

for i in range(len(T)):

  vX[i]=get_non_bondend_sites_from_states(stateResult_poor[i])

#%%

i=0
AGUA_NEGATIVO=np.zeros_like(stateResult_poor[i],dtype=object)
AGUA_POSITIVO=np.zeros_like(stateResult_poor[i],dtype=object)
CO2_POSITIVO=np.zeros_like(stateResult_poor[i],dtype=object)
yW=np.zeros_like(AGUA_NEGATIVO)
yCO2=np.zeros_like(AGUA_NEGATIVO)

for j,state in enumerate(stateResult_poor[i]):

    # idx 0 - agua negativo
    # idx 1 - agua positivo 
    # idx 2 - co2 positivo
    AGUA_NEGATIVO[j],AGUA_POSITIVO[j],CO2_POSITIVO[j]=vX[i][j]
    
    yW[j],yCO2[j]=state.composition()

#%%
AGUA_NEGATIVO=AGUA_NEGATIVO
AGUA_POSITIVO=AGUA_POSITIVO
CO2_POSITIVO=CO2_POSITIVO


#%%

plt.plot(pResultBAR[i],AGUA_NEGATIVO,color="black",linestyle=lines[0],label="H2O -")
plt.plot(pResultBAR[i],AGUA_POSITIVO,color="black",linestyle=lines[1],label="H2O +")
plt.plot(pResultBAR[i],CO2_POSITIVO, color="black",linestyle=lines[2],label="CO2 +")
plt.xlabel("P/bar")
plt.ylabel("Fraction of Non-Bonded Sites")
plt.xlim(0,400)
plt.text(150,0.8,"323.15K",fontsize=12)
plt.legend()
plt.title("Water 4C and CO2 1ea")
plt.savefig("plots/water_co2/X_323_vap.pdf",bbox_inches='tight')

#%% liq
vco2_per_gram_water=np.array([
    9.71,
    17.25,
    22.53,
    25.63,
    26.77,
    27.64,
    29.14,
    31.34,
    33.29
])

x=((18*1e5)/(8.314*273.15))*vco2_per_gram_water*1e-6

xco2=x/(1+x)

p_exp=np.array([
    25,
    50,
    75,
    100,
    125,
    150,
    200,
    300,
    400
])

for (i,t) in enumerate(T):


  pe,_=wsat_data[t]

  xco2_result=[ x.composition()[1] for x in stateResult_rich[i] ]

  # plt.plot(pResultBAR[i],yResultPPM[i],)
  plt.plot(pResultBAR[i],
           xco2_result,
           color="Black",
           linestyle=lines[i],
           label=f"{t}K")
  
  plt.scatter(p_exp,
              xco2,
              marker=markers[i],
              facecolors='none', 
              edgecolors='black') #COLOCAR PONTOS XCO2

  plt.ylim(0,0.03
           )
  # plt.ylabel("Water Mole Fraction (ppm)")
  plt.xlabel("P/bar")

  plt.xlim(0,410)

plt.title(f"H2O 4C and CO2 1ea")
plt.legend()


plt.savefig("plots/water_co2/wsat_co2_323_liq_phase.pdf",bbox_inches='tight')
# plt.text(500,17500,"H2O 4C, CO2 1ea")

#%%
vX=np.zeros_like(T,dtype=object)

for i in range(len(T)):

  vX[i]=get_non_bondend_sites_from_states(stateResult_rich[i])

#%%

i=0
AGUA_NEGATIVO=np.zeros_like(stateResult_rich[i],dtype=object)
AGUA_POSITIVO=np.zeros_like(stateResult_rich[i],dtype=object)
CO2_POSITIVO=np.zeros_like(stateResult_rich[i],dtype=object)
yW=np.zeros_like(AGUA_NEGATIVO)
yCO2=np.zeros_like(AGUA_NEGATIVO)

for j,state in enumerate(stateResult_rich[i]):

    # idx 0 - agua negativo
    # idx 1 - agua positivo 
    # idx 2 - co2 positivo
    AGUA_NEGATIVO[j],AGUA_POSITIVO[j],CO2_POSITIVO[j]=vX[i][j]
    
    yW[j],yCO2[j]=state.composition()

# plt.plot(xResult)

#%%
plt.plot(pResultBAR[i],AGUA_NEGATIVO,color="black",linestyle=lines[0],label="H2O -")
plt.plot(pResultBAR[i],AGUA_POSITIVO,color="black",linestyle=lines[1],label="H2O +")
plt.plot(pResultBAR[i],CO2_POSITIVO, color="black",linestyle=lines[2],label="CO2 +")
plt.xlabel("P/bar")
plt.ylabel("Fraction of Non-Bonded Sites")
# plt.xlim(0,400)
# plt.text(150,0.8,"323.15K",fontsize=12)
plt.legend()
plt.title("Water 4C and CO2 1ea (Liquid Phase)")
plt.savefig("plots/water_co2/X_323_liq_phase.pdf",bbox_inches='tight')

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