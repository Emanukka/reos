import numpy as np
from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import scipy.optimize as opt
tol_tpd=1e-6

def tpd_root_given_tp(eos,eos_pure_water,T,P,y_dry_gas):

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

        return LIQUID.min_tpd("vapor",xguess,tol=tol_tpd),LIQUID
    else:
       
        return VAPOR.min_tpd("liquid",xguess,tol=tol_tpd),VAPOR
    
      
    # dg,states=state_z2.min_tpd("liquid",xguess)
    # dg,states=state_z.min_tpd("liquid",xguess,)



  # print(zguess)
  fun=lambda X: F(X)[0][0]

  yw_log_root=opt.root(fun,[yw_guess_log] ).x

  # print(pbol)
  tpd_res,zstate=F([yw_log_root])

  xstate=tpd_res[1]
  yw=np.exp(yw_log_root)

  return yw,xstate,zstate


def linspace_wsat(eos,eos_pure_water,T,ydg,pi=1,pf=600,N=100):
   
  PRESS_BAR=np.linspace(pi,pf,N)
  yw_PPM=np.zeros_like(PRESS_BAR)
  zSTATE=np.zeros_like(PRESS_BAR,dtype=object)
  xSTATE=np.zeros_like(PRESS_BAR,dtype=object)

  for (i,P_BAR) in enumerate(PRESS_BAR):

    try:
      # xw=xw*1e-6
      # z=np.array([xw,1-xw])

      P_PA=P_BAR*1e5

      result,xstate,zstate=tpd_root_given_tp(eos,eos_pure_water,T,P_PA,ydg)

      zSTATE[i]=zstate
      xSTATE[i]=xstate


      yw_PPM[i]=result[0]*1e6

    except Exception as e:
      print(e)
      pass
      continue
  
  return PRESS_BAR,yw_PPM, zSTATE,xSTATE