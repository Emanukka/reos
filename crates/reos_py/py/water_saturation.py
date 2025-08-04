
#%%


from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from auxiliary_functions.parameters import *
from auxiliary_functions.vle_functions import *
from auxiliary_functions.data import *

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

tol_tpd=1e-6
def tpd_root_given_tp(T,P,y_dry_gas):

  def RES(T,P,v):
    
    
    yw=np.exp(v[0])
    # print(yw)

    z=(y_dry_gas/y_dry_gas.sum())*(1.-yw)
    z[0]=yw
    # print(sum(y))
    # print(yw)
    
    VAPOR=State.tpx(eos,T,P,z,"vapor")
    LIQUID=State.tpx(eos,T,P,z,"liquid")

    muL=LIQUID.ln_phi()+np.log(z)
    muV=VAPOR.ln_phi()+np.log(z)

    GL=np.sum(muL*z)
    GV=np.sum(muV*z)

    #Com certeza s2 é liquid, pois a pressao é muito alta
    s2=State.tpx(eos_pure_water,T,P,np.array([1.0]),"liquid")

    if GL<GV:
        s1=LIQUID
        return ( np.log(yw)+s1.ln_phi()[0] - s2.ln_phi()[0] )

    else:
        s1=VAPOR
       

        return ( np.log(yw)+s1.ln_phi()[0] - s2.ln_phi()[0] )

  xw=0.9999
  N=len(y_dry_gas)

  xguess=np.zeros(N)
  xguess[:]=(1.0-xw)/(N-1)
  xguess[0]=xw
  fRES = lambda X: RES(T,P,X)

  # print("aqui")
  # print(sum(incipient_phase_guess))
  yguess_log=np.log(2000*1e-6)
  yw_guess_log=opt.root(fRES,yguess_log).x[0]

  # zguess=np.array([yw_guess,1-yw_guess])

  # print("guess=",yw_guess)
#   print(np.exp(yw_guess_log))

  def F(v):

    # print(v)
    yw=np.exp(v[0])
    Z=(y_dry_gas/y_dry_gas.sum())*(1.-yw)

    Z[0]=yw
    # print("Z=",Z)
    LIQUID=State.tpx(eos,T,P,Z,"liquid")
    VAPOR=State.tpx(eos,T,P,Z,"vapor")

    muL=LIQUID.ln_phi()+np.log(Z)
    muV=VAPOR.ln_phi()+np.log(Z)

    GL=np.sum(muL*Z)
    GV=np.sum(muV*Z)

    if GL<GV:

        # mudar chute, poix xguess é pra quando é 
        return LIQUID.min_tpd("vapor",xguess,tol=tol_tpd)
    else:
       
        return VAPOR.min_tpd("liquid",xguess,tol=tol_tpd)
    
      
    # dg,states=state_z2.min_tpd("liquid",xguess)
    # dg,states=state_z.min_tpd("liquid",xguess,)



  # print(zguess)
  fun=lambda X: F(X)[0]

  yw_log_root=opt.root(fun,[yw_guess_log] ).x

  # print(pbol)

  return np.exp(yw_log_root)

#%%
T=298.15
#
PRESS_BAR=np.linspace(1,600,100)
yw=np.zeros_like(PRESS_BAR)
# press_orv=np.zeros_like(water_PPM)
# press_bolha=np.zeros_like(water_PPM)
y_dry_gas=np.array([0.0,0.5])
# y_dry_gas=np.array([0.0,1.0,1e-20])

for (i,P_BAR) in enumerate(PRESS_BAR):

  try:
    # xw=xw*1e-6
    # z=np.array([xw,1-xw])

    P_PA=P_BAR*1e5

    result=tpd_root_given_tp(T,P_PA,y_dry_gas)

    yw[i]=result[0]

  except Exception as e:
    print(e)
    pass
    continue
#%%
DataKonteo=dict()

# 298K
# DataKonteo[298]=(

# [1.01   ,25.33	 ,50.66	 ,101.33 ,151.99 ,202.65 ,455.96 ,481.29 ,506.63 ,5.04 ,10.07,14.96,24.83  ,34.91  ,16.66  ,32.58  ,63.91  ]

# ,

# [30245,1736,1368,3511,3800,3984,4234,4221,4195,6500,3800,2800,1400,1300,2680,1600,1940]

# )

plt.scatter(PRESS_BAR,yw*1e6)
# vx

plt.ylim(0,20*1e3)
plt.xlim(20,100)
plt.title(f"T={T}K,")

# pexp,yexp=DataKonteo[298][0],DataKonteo[298][1]
# plt.scatter(pexp,yexp,label=f"{T}K;Konteo. ref")

#%%
#%%
DataKonteo=dict()

# 298K
DataKonteo[298]=(

[1.01   ,25.33	 ,50.66	 ,101.33 ,151.99 ,202.65 ,455.96 ,481.29 ,506.63 ,5.04 ,10.07,14.96,24.83  ,34.91  ,16.66  ,32.58  ,63.91  ]

,

[30245,1736,1368,3511,3800,3984,4234,4221,4195,6500,3800,2800,1400,1300,2680,1600,1940]

)

plt.plot(PRESS_BAR,yw*1e6)
# vx

plt.ylim(0,20*1e3)

# plt.xlim(0,125)
# plt.xticks([50,75,100])
# plt.title(f"T={T}K, CO2 puro")

pexp,yexp=DataKonteo[298][0],DataKonteo[298][1]
plt.scatter(pexp,yexp,label=f"{T}K;Konteo. ref")