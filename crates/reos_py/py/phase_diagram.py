#%%

from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from auxiliary_functions.parameters import *
from auxiliary_functions.vle_functions import *
from auxiliary_functions.data import *


#%%
pACOH_OCT=CPAParameters.from_records(
    cubic=[c_acoh,c_octane],
    assoc=[a_acoh,a_octane])

pACOH_OCT.set_cubic_binary(0,1,0.0, 0.064)

ACOH_OCT=EquationOfState.cpa(pACOH_OCT)# print(pWATER_ACETIC.as_string())

# peq=PhaseEquilibrium(ACOH_OCT)

antoine=np.array([acoh_antoine,octane_antoine])
# p,y,vx=vle_diagram(T,peq,factor=1e3)

xorv,porv=acoh_octane["orv"]
xbol,pbol=acoh_octane["bol"]
exp_data=[xorv,porv,xbol,pbol]

T=343.2

PRES,LIQUID,VAPOR=linspace_bubble_p(ACOH_OCT,T,antoine,N=100)

bubble_diagram(
   PRES,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=15,
   y_sup=30,
   title="Acetic Acid 1A(1) and Octane(2)",
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   exp_data=exp_data)



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
exp_data=[xorv,porv,xbol,pbol]

VAR,LIQUID,VAPOR=linspace_bubble_p(PROPANOIC_HEPTANE,T,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=0,
   y_sup=25,
   y_label="P/kPa",
   x_label=r"$x_1,y_1$",
   title="Propanoic Acid(1) and Heptane",
   exp_data=exp_data)


#%%
p=CPAParameters.from_records(
    cubic=[c_methanol_2b,c_octanol_2b],
    assoc=[a_methanol_2b,a_octanol_2b])

p.set_cubic_binary(0,1,0.0,0.0)
p.set_assoc_binary(0,1,"cr1")
eos_2b=EquationOfState.cpa(p)


antoine=np.array([metoh_antoine,octanol_antoine])
P=101.32e3

xorv,torv=metoh_otctanol["orv"]
xbol,tbol=metoh_otctanol["bol"]

exp_data=[xorv,torv,xbol,tbol]


VAR,LIQUID,VAPOR=linspace_bubble_t(eos_2b,P,antoine,N=100)

bubble_diagram(VAR,LIQUID,VAPOR,factor=1.0,y_figsize=2.5,x_figsize=5,y_inf=300,y_sup=480,exp_data=exp_data)


#%%
p=CPAParameters.from_records(
    cubic=[c_methanol_3b,c_octanol_3b],
    assoc=[a_methanol_3b,a_octanol_3b])

p.set_cubic_binary(0,1,0.0,-0.025)
p.set_assoc_binary(0,1,"ecr")
eos_3b=EquationOfState.cpa(p)

VAR,LIQUID,VAPOR=linspace_bubble_t(eos_3b,P,antoine,N=100)

bubble_diagram(VAR,LIQUID,VAPOR,factor=1.0,y_figsize=2.5,x_figsize=5,y_inf=300,y_sup=480,exp_data=exp_data)

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



VAR,LIQUID,VAPOR=linspace_bubble_t(PROPANOIC_HEPTANE,P,antoine,N=100)

bubble_diagram(
   VAR,
   LIQUID,
   VAPOR,
   factor=1.0,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=360,
   y_sup=420,
   title="Propanoic Acid 1A(1) and Heptane(2)",
   exp_data=exp_data)


#%%

pMETHANOL_ACETIC=CPAParameters.from_records(
    cubic=[c_methanol_2b,c_acoh],
    assoc=[a_methanol_2b,a_acoh])

pMETHANOL_ACETIC.set_cubic_binary(0,1,0.0,-0.04)
pMETHANOL_ACETIC.set_assoc_binary(0,1,"ecr")
METHANOL_ACETIC=EquationOfState.cpa(pMETHANOL_ACETIC)


antoine=np.array([metoh_antoine,acoh_antoine])
T=308.15

xorv,porv=metoh_acoh["orv"]
xbol,pbol=metoh_acoh["bol"]
exp_data=[xorv,porv,xbol,pbol]

PRES,LIQUID,VAPOR=linspace_bubble_p(METHANOL_ACETIC,T,antoine,N=100)

bubble_diagram(
   PRES,
   LIQUID,
   VAPOR,
   factor=1e3,
   y_figsize=2.5,
   x_figsize=5,
   y_inf=0.0,
   y_sup=30,
   title="Methanol 2B(1) and AcOH 1A(2) (ECR)",
   exp_data=exp_data)


# BOL,_,linspaceZ,XASCL,XASCV=VLE_DIAGRAM(
#     ("t",T),
#     METHANOL_ACETIC,antoine,
#     y_label=r"P/kPa",
#     x_label=r"$x_1,y_1$",
#     y_lim=[0,30],
#     x_figsize=5,
#     y_figsize=5,
#     factor=1e3,
#     exp_data=exp_data,
#     title=r"Methanol 2B(1) and AcOH 1A(2) (ECR)",
#     save_fig=False,
#     N_points=100)

# #now testando modelos de grafico de X
# XASCV=np.stack(XASCV)
# XASCL=np.stack(XASCL)

# plt.figure()
# # plt.plot(BOL,XASC)
# # plt.plot(linspaceZ,XASCV[:,0],label="- MetOH VAP")
# # plt.plot(linspaceZ,XASCV[:,1],label="+ MetOH VAP")
# # plt.plot(linspaceZ,XASCV[:,2],label="+-AcOH  VAP")

# # plt.plot(linspaceZ,XASCL[:,0],label="-MetOH LIQ")
# # plt.plot(linspaceZ,XASCL[:,1],label="+MetOH LIQ")
# # plt.plot(linspaceZ,XASCL[:,2],label="+-AcOH LIQ")

# plt.plot(BOL,XASCV[:,0],label="- MetOH VAP")
# plt.plot(BOL,XASCV[:,1],label="+ MetOH VAP")
# plt.plot(BOL,XASCV[:,2],label="+-AcOH  VAP")
# plt.plot(BOL,XASCL[:,0],label="-MetOH LIQ")
# plt.plot(BOL,XASCL[:,1],label="+MetOH LIQ")
# plt.plot(BOL,XASCL[:,2],label="+-AcOH LIQ")

# plt.legend()
# # plt.xlim(-0.01,1.01)
# plt.ylim(-0.01,1.01)


#%%

co2=CPAParameters.from_records([c_co2],[a_co2])

eos=EquationOfState.cpa(co2)

x=np.array([1.0])
# state=State.tpx(eos,288,1e5,x,"vapor")
# phi=np.exp(state.ln_phi())

T=288.15
def F(P):
    state=State.tpx(eos,T,P,x,"liquid")
    lnphiL=state.ln_phi()*1

    state=State.tpx(eos,T,P,x,"vapor")
    lnphiV=state.ln_phi()*1

    print(lnphiL,lnphiV)
    return lnphiL-lnphiV
F(50e5)
from scipy import optimize as op
ans = opt.root(F,50e5)
print(ans)

# %%


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
