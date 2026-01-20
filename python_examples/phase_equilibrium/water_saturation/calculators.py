#%%
import numpy as np
from reos.cpa import CPAParameters
from reos.state import State
import scipy.optimize as opt
from si_units import BAR, KELVIN, PASCAL

#%%

def tpd_root_given_tp(eos, eos_water, T, P, y_dry_gas):

  def RES(T,P,v):
    
    
    yw = np.exp(v[0])
    z = y_dry_gas / y_dry_gas.sum() * (1.0 - yw)
    z[0] = yw
    
    s1 = State.tpx(eos, T, P, z) 
    # VAPOR=State.tpx(eos,T,P,z,"vapor")
    # LIQUID=State.tpx(eos,T,P,z,"liquid")

    # muL=LIQUID.ln_phi()+np.log(z)
    # muV=VAPOR.ln_phi()+np.log(z)
    # GL=np.sum(muL*z)
    # GV=np.sum(muV*z)

    s2 = State.tpx(eos_water, T, P, np.array([1.0]),"liquid")
    
    return np.log(yw) + s1.ln_phi()[0] - s2.ln_phi()[0]
    # if GL < GV:

    #     s1=LIQUID
    #     return ( np.log(yw) + s1.ln_phi()[0] - s2.ln_phi()[0] )

    # else:

    #     s1=VAPOR
    #     return ( np.log(yw) + s1.ln_phi()[0] - s2.ln_phi()[0] )

  xw = 0.9999
  N = len(y_dry_gas)

  xguess = np.zeros(N)
  xguess[:] = (1.0-xw)/(N-1)
  xguess[0] = xw
  fRES = lambda X: RES(T,P,X)

  yguess_log = np.log(2000*1e-6)
  yw_guess_log = opt.root(fRES,yguess_log).x[0]

  # zguess=np.array([yw_guess,1-yw_guess])

  # print("guess=",yw_guess)
#   print(np.exp(yw_guess_log))

  def F(v):

    # print(v)
    yw = np.exp(v[0])
    z = (y_dry_gas/y_dry_gas.sum())*(1.-yw)
    z[0] = yw
    # print("Z=",Z)

    developed = State.tpx(eos, T, P, z)
    # LIQUID = State.tpx(eos,T,P,Z,"liquid")
    # VAPOR = State.tpx(eos,T,P,Z,"vapor")

    # muL=LIQUID.ln_phi()+np.log(Z)
    # muV=VAPOR.ln_phi()+np.log(Z)

    # GL=np.sum(muL*Z)
    # GV=np.sum(muV*Z)

    dg = min_tpd.dg
    incipient = min_tpd.state
    
    return dg, np.array([developed, incipient])
  
    # if GL<GV: 

    #     return LIQUID.min_tpd("vapor",xguess,tol=tol_tpd),LIQUID
    # else:
       
    #     return VAPOR.min_tpd("liquid",xguess,tol=tol_tpd),VAPOR
    
      
    # dg,states=state_z2.min_tpd("liquid",xguess)
    # dg,states=state_z.min_tpd("liquid",xguess,)



  # print(zguess)
  fun = lambda X: F(X)[0]

  yw_log_root = opt.root(fun,[yw_guess_log] ).x

  # print(pbol)
  dg, states = F([yw_log_root])

  # xstate = states[1]
  # xstate = tpd_res[1]
  yw = np.exp(yw_log_root)

  return yw, states


def linspace_wsat(
    eos,
    eos_water,
    t,
    y_dry_gas, 
    pi = 1, 
    pf = 600, 
    N = 100):
   
  vpressure = np.linspace(pi, pf, N)
  yw_ppm = np.zeros(N)
  vdeveloped = np.zeros(N, dtype = object)
  vincipient = np.zeros(N, dtype = object)

  for (i, p) in enumerate(vpressure):

    try:
      # xw=xw*1e-6
      # z=np.array([xw,1-xw])

      p = p * 1e5

      result, states = tpd_root_given_tp(eos, eos_water, t, p, y_dry_gas)

      vdeveloped[i] = states[0]
      vincipient[i] = states[1]

      yw_ppm[i] = result[0] * 1e6

    except Exception as e:
      print(e)
      pass
      continue
  
  return vpressure, yw_ppm, vdeveloped, vincipient
# %%
