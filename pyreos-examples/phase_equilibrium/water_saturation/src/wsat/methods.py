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
    
    s1 = State.tpx(eos, T, P, z,"stable") 

    s2 = State.tpx(eos_water, T, P, np.array([1.0]), "liquid")

    lnphi1 = s1.lnphi()[0]
    lnphi2 = s2.lnphi()[0]

    return np.log(yw) + lnphi1  - lnphi2

  xw = 0.9999
  N = len(y_dry_gas)

  xguess = np.zeros(N)
  xguess[:] = (1.0-xw)/(N-1)
  xguess[0] = xw
  fRES = lambda X: RES(T,P,X)

  yguess_log = np.log(2000*1e-6)
  res = opt.root(fRES,yguess_log).x
  yw_guess_log = res[0]

  def F(v):

    yw = np.exp(v[0])
    z = (y_dry_gas/y_dry_gas.sum())*(1.-yw)
    z[0] = yw

    liq = State.tpx(eos, T, P, z, "liquid")
    vap = State.tpx(eos, T, P, z, "vapor")

    if liq.gibbs() > vap.gibbs():

      min_tpd = vap.min_tpd("liquid", xguess, 1e-9, 200)
      dg = min_tpd.dg
      incipient = min_tpd.state

      return dg, np.array([vap, incipient])
    
    else:

      min_tpd = liq.min_tpd("vapor", xguess, 1e-9, 200)
      dg = min_tpd.dg
      incipient = min_tpd.state
      
      return dg, np.array([liq, incipient])
    
  fun = lambda X: F(X)[0]
  yw_log_root = opt.root(fun,[yw_guess_log] ).x[0]
  _, states = F([yw_log_root])
  yw = np.exp(yw_log_root)

  return yw, states


def linspace_wsat(
    eos,
    eos_water,
    t,
    y_dry_gas, 
    p0, 
    pf, 
    N):
   
  vpressure = np.linspace(p0, pf, N)
  yw = np.zeros(N)
  vdeveloped = np.zeros(N, dtype = object)
  vincipient = np.zeros(N, dtype = object)

  for (i, p) in enumerate(vpressure):

    try:
      # xw=xw*1e-6
      # z=np.array([xw,1-xw])

      yw_frac, states = tpd_root_given_tp(eos, eos_water, t, p, y_dry_gas)

      vdeveloped[i] = states[0]
      vincipient[i] = states[1]

      yw[i] = yw_frac 

    except Exception as e:
      print(e)
      pass
      continue
  
  return vpressure, yw, vdeveloped, vincipient
