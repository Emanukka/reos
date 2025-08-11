#%%
from reos.reos  import EquationOfState,State,CPAParameters,CubicParameters,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
from auxiliary_functions.parameters import *

plt.rcParams.update({
    "text.usetex": True,               
    "font.family": "serif",            
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 12,
    "font.size": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

#%%

pcpa=CPAParameters.from_records(cubic=[c_w],assoc=[a_w])
psrk=CubicParameters.from_records([c_w])
# print(p1.as_string())
srk=EquationOfState.srk(psrk)
cpa=EquationOfState.cpa(pcpa)

t=298.15
x=np.array([1.0])
vp=np.linspace(0.01,700,100)


z_srk=[State.tpx(srk,t,p*1e5,x,"vapor").compressibility() for p in vp]

z_cpa=[State.tpx(cpa,t,p*1e5,x,"vapor").compressibility() for p in vp]



plt.plot(vp,z_srk,linestyle="--",color="black")
plt.plot(vp,z_cpa,linestyle="-",color="black")