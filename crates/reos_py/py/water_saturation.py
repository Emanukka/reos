
#%%


from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
from auxiliary_functions.parameters import *
# from auxiliary_functions.vle_functions import *
# from auxiliary_functions.data import *
from auxiliary_functions.wsat_data import *


#%%


#%%
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


def linspace_wsat(eos,T,ydg,pi=1,pf=600,N=100):
   
  PRESS_BAR=np.linspace(1,600,100)
  yw_PPM=np.zeros_like(PRESS_BAR)

  for (i,P_BAR) in enumerate(PRESS_BAR):

    try:
      # xw=xw*1e-6
      # z=np.array([xw,1-xw])

      P_PA=P_BAR*1e5

      result=tpd_root_given_tp(T,P_PA,ydg)

      yw_PPM[i]=result[0]*1e6

    except Exception as e:
      print(e)
      pass
      continue
  
  return PRESS_BAR,yw_PPM


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
#%%
T=np.array([298.15,308.15])

ydg=np.array([0.0,0.5])
pResultBAR=np.zeros_like(T,dtype=object)
yResultPPM=np.zeros_like(T,dtype=object)
for (i,temp) in enumerate(T):

  try:
     
    p,y=linspace_wsat(eos,temp,ydg)
    pResultBAR[i]=p*1
    yResultPPM[i]=y*1

  except Exception as e:
    print(e)
    pass
    continue
#%%

i=0

pe,ye=wsat_data[T[i]]
# plt.plot(pResultBAR[i],yResultPPM[i],)
plt.scatter(pResultBAR[i],yResultPPM[i],)

# plt.scatter(pe,ye*1e6)
plt.ylim(0,20*1e3)
# plt.xlim(20,100)
plt.title(f"T={T[i]}K,")

yw_ppm=[32214.659481216913, 4858.681097055392, 2785.770691689941, 2035.0246242290398, 1656.0045230218059, 1434.8097516552368, 1296.9983919603833, 1210.7008752769873, 1161.3568388463696, 2942.3771211249136, 3221.951384091279, 1318.3409457948806, 3536.0774604583166, 3651.195424693426, 3751.3719547623723, 3840.664832834879, 3921.577690287392, 3995.7835604519437, 4064.4651707590065, 4128.495482652775, 4188.540018111133, 4245.1204339591895, 4298.654156307189, 4349.4818486102395, 4397.885287097065, 4444.101107442424, 4488.329880866098, 4530.743864274891, 4571.4919936692795, 4610.704056164973, 4648.494320826769, 4684.9634912550655, 4720.201226744937, 4754.287732013092, 4787.294802016899, 4819.287434223879, 4850.324428685772, 4880.459173750209, 4909.7405944818975, 4938.213266587783, 4965.918183730299, 4992.893089489673, 5019.172892723248, 5044.78970592505, 5069.773505311598, 5094.152095647827, 5117.951543371733, 5141.195636592322, 5163.907261270164, 5186.107532773146, 5207.816318287103, 5229.052287694755, 5249.833064750372, 5270.175103979498, 5290.09408739646, 5309.604784152106, 5328.721168139159, 5347.456476459282, 5365.8232154150855, 5383.833397268825, 5401.49830104004, 5418.828830113863, 5435.83498674815, 5452.526679301892, 5468.913236412998, 5485.003564720465, 5500.80613394937, 5516.329165074414, 5531.5804061204635, 5546.5673170889295, 5561.297047818402, 5575.776455592887, 5590.012158318885, 5604.010402260867, 5617.777297557343, 5631.318683654663, 5644.64017807249, 5657.7471897957475, 5670.644927018088, 5683.338376298984, 5695.832443077156, 5708.131762755063, 5720.240840901424, 5732.164027951713, 5743.905526440803, 5755.4693980918855, 5766.859569779216, 5778.079869104614, 5789.1339135080125, 5800.025290542297, 5810.757446253469, 5821.333719214388, 5831.757345077886, 5842.0314608873605, 5852.159083402551, 5862.143216119424, 5871.986697989711, 5881.692310340604, 5891.262754278138, 5900.700653821883]
PRES_BAR=[1.0, 7.05050505050505, 13.1010101010101, 19.15151515151515, 25.2020202020202, 31.252525252525253, 37.3030303030303, 43.35353535353535, 49.4040404040404, 55.45454545454545, 61.505050505050505, 67.55555555555556, 73.6060606060606, 79.65656565656565, 85.7070707070707, 91.75757575757575, 97.8080808080808, 103.85858585858585, 109.9090909090909, 115.95959595959596, 122.01010101010101, 128.06060606060606, 134.11111111111111, 140.16161616161617, 146.2121212121212, 152.26262626262624, 158.3131313131313, 164.36363636363635, 170.4141414141414, 176.46464646464645, 182.5151515151515, 188.56565656565655, 194.6161616161616, 200.66666666666666, 206.7171717171717, 212.76767676767676, 218.8181818181818, 224.86868686868686, 230.91919191919192, 236.96969696969697, 243.02020202020202, 249.07070707070704, 255.1212121212121, 261.17171717171715, 267.22222222222223, 273.27272727272725, 279.32323232323233, 285.37373737373736, 291.4242424242424, 297.47474747474746, 303.5252525252525, 309.57575757575756, 315.6262626262626, 321.67676767676767, 327.7272727272727, 333.77777777777777, 339.8282828282828, 345.8787878787879, 351.9292929292929, 357.979797979798, 364.030303030303, 370.0808080808081, 376.1313131313131, 382.1818181818182, 388.2323232323232, 394.28282828282823, 400.3333333333333, 406.38383838383834, 412.4343434343434, 418.48484848484844, 424.5353535353535, 430.58585858585855, 436.6363636363636, 442.68686868686865, 448.73737373737373, 454.78787878787875, 460.83838383838383, 466.88888888888886, 472.93939393939394, 478.98989898989896, 485.04040404040404, 491.09090909090907, 497.1414141414141, 503.19191919191917, 509.2424242424242, 515.2929292929292, 521.3434343434343, 527.3939393939394, 533.4444444444445, 539.4949494949494, 545.5454545454545, 551.5959595959596, 557.6464646464647, 563.6969696969696, 569.7474747474747, 575.7979797979798, 581.8484848484848, 587.8989898989898, 593.9494949494949, 600.0]

plt.scatter(PRES_BAR,yw_ppm)

# i=1
# pe,ye=wsat_data[T[i]]
# plt.plot(pResultBAR[i],yResultPPM[i],)
# plt.scatter(pe,ye*1e6)
# plt.ylim(0,20*1e3)
# # plt.xlim(20,100)
# plt.title(f"T={T[i]}K,")

#%%
#%%



pWATER_CO2=CPAParameters.from_records(
    cubic=[c_w,c_co2],
    assoc=[a_w,a_co2])



pWATER_CO2.set_cubic_binary(0,1,kij_a=0.0,kij_b=0.1141)

pWATER_CO2.set_assoc_binary(0,1,"exp",eps=142e2,beta=0.0162)


eos=EquationOfState.cpa(pWATER_CO2)

T=np.array([298.15,308.15])

ydg=np.array([0.0,0.5])
pResultBAR=np.zeros_like(T,dtype=object)
yResultPPM=np.zeros_like(T,dtype=object)
for (i,temp) in enumerate(T):

  try:
     
    p,y=linspace_wsat(eos,temp,ydg)
    pResultBAR[i]=p*1
    yResultPPM[i]=y*1

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