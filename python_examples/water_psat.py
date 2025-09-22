#%%
from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
from auxiliary_functions.parameters import *
from auxiliary_functions.association_functions import *

import matplotlib.pyplot as plt

parameters = CPAParameters.from_records(
    cubic = [c_w],
    assoc = [a_w])


eos = EquationOfState.cpa(parameters)



#%%

assoc=eos.get_association()

# array with the all sites in mixture
# 0=Site of Type A
# 1=Site of Type B
# 2=Site of Type C

F=assoc.get_sites_map()
T=assoc.get_tmat()

#%%

T=np.linspace(298.15,600,50)


def fun(p):

    s1=State.tpx(eos,t,p,np.array([1.0]),'vapor')
    s2=State.tpx(eos,t,p,np.array([1.0]),'liquid')

    phiV=np.exp(s1.ln_phi())
    phiL=np.exp(s2.ln_phi())

    d=abs(phiV-phiL)

    return d

#%%

t=310.15

P=np.linspace(1e5,1e6)
res=np.zeros_like(P)

for (i,p) in enumerate(P):

    res[i]=fun(p)


#%%
plt.plot(res,)
#%%
for (i,t) in enumerate(T):

    
