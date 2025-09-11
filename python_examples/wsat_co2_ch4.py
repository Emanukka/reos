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
p=CPAParameters.from_records(
    cubic=[c_w,c_co2,c_ch4],
    assoc=[a_w,a_co2,a_ch4])

pure_water=CPAParameters.from_records(
   [c_w],[a_w]

)

eos_pure_water=EquationOfState.cpa(pure_water)

p.set_cubic_binary(0,1,kij_a=0.000877,kij_b=-0.15508)

p.set_assoc_binary(0,1,"mcr1",beta=0.1836)

# markers = ['o', 's', '^', 'D', 'P', 'X'] 
# lines= ['-','--',':']
eos=EquationOfState.cpa(p)

# press={}
# press[50]=[2,25,]
# T=np.array([298.15])




#%%
import io
import numpy as np
tabela="""
50	25	2
50	30	2
50	35	2
50	40	2
50	45	25
55	25	25
55	30	25
55	35	25
55	40	40
55	45	40
60	25	40
60	30	40
60	35	60
60	40	60
60	45	60
70	25	2
70	30	2
70	35	2
70	40	2
70	45	2
75	25	25
75	30	25
75	35	25
75	40	25
75	45	25
80	25	40
80	30	40
80	35	40
80	40	40
80	45	40
85	25	60
85	30	60
85	35	60
85	40	60
85	45	60
120	35	6
120	40	6
120	45	6
150	35	6
150	40	6
150	45	6
180	35	6
180	40	6
180	45	6
180	20	3
180	10	3
180	4	3
200	20	3
200	10	3
200	4	3
220	20	3
220	10	3
220	4	3
250	20	3
250	10	3
250	4	3
550	20	25
550	10	60
550	4	80

"""

dados = np.loadtxt(io.StringIO(tabela), delimiter='\t')
#%%
NL,_=dados.shape
# print(dados.shape)
# print(dados[2,:])
    

yw_result=np.zeros(NL)
#%%
for i in range(NL):

    dado_i=dados[i,:]
    P_BAR,T_CELSIUS,yco2=dado_i

    P_PA=1e5*P_BAR
    T=T_CELSIUS+273.15
    yco2=yco2*0.01
    ydg=np.array([0.0,yco2,1.0-yco2])

    result,_=tpd_root_given_tp(eos,eos_pure_water,T,P_PA,ydg)

    yw_result[i]=result*1e6
    # print(P_PA,T,yco2,dado_i)

    # print(result*1e6)

#%%

P_PA=1e5*325
T=3+273.15
yco2=3*0.01
ydg=np.array([0.0,yco2,1.0-yco2])
result,_=tpd_root_given_tp(eos,eos_pure_water,T,P_PA,ydg)
result*1e6
#%%
print(yw_result)
# # %%
# ydg=np.array([0.0,0.02,0.98])
# T=[298.15]
# pResultBAR=np.zeros_like(T,dtype=object)
# yResultPPM=np.zeros_like(T,dtype=object)
# stateResult=np.zeros_like(T,dtype=object)

# for (i,temp) in enumerate(T):

#   try:
     
#     p,y,s=linspace_wsat(eos,eos_pure_water,temp,ydg)
#     pResultBAR[i]=p*1
#     yResultPPM[i]=y*1
#     stateResult[i]=s

#   except Exception as e:
#     print(e)
#     pass
#     continue
  
#%%

# P_PA=1e5*325
# T=3.0+273.15
# yco2=yco2*0.01
# ydg=np.array([0.0,yco2,1.0-yco2])
# result,_=tpd_root_given_tp(eos,eos_pure_water,T,P_PA,ydg)
    