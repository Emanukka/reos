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



#%%



#%%
p=CPAParameters.from_records(
    cubic=[c_w,c_h2s],
    assoc=[a_w,a_h2s])

pure_water=CPAParameters.from_records(
   [c_w],[a_w]

)

eos_pure_water=EquationOfState.cpa(pure_water)

p.set_cubic_binary(0,1,kij_a=0.0,kij_b=0.1913)

p.set_assoc_binary(0,1,"exp",eps=108.78e2,beta=0.0624)


eos=EquationOfState.cpa(p)
#%%calc
T=np.array([310.9])
# T=np.array([298.15])

ydg=np.array([0.0,0.5])
pResultBAR=np.zeros_like(T,dtype=object)
yResultPPM=np.zeros_like(T,dtype=object)
stateResult=np.zeros_like(T,dtype=object)

for (i,temp) in enumerate(T):

  try:
     
    p,y,s=linspace_wsat(eos,eos_pure_water,temp,ydg,pi=5,pf=250)
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
plt.figure(figsize=(6, 5))

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

  plt.ylim(0,25*1e3)
  plt.ylabel("Water Mole Fraction (ppm)")
  plt.xlabel("P/bar")

  # plt.xlim(20,100)

# plt.title(f"CO2 1ea")

plt.legend()