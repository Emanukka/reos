from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord
import matplotlib.pyplot as plt
import scipy.optimize as opt
import os
import numpy as np
def bubble_diagram(

    VAR,
    LIQUID,
    VAPOR,
    plot_dir,
    x_figsize=5,
    y_figsize=5,
    y_inf=None,
    y_sup=None,
    x_inf=0.0,
    x_sup=1.0,
    x_label="",
    y_label="",
    title="",
    text="",
    bol_linestyle="-",
    orv_linestyle="-",
    exp_data=[None,None,None,None],
    save_fig=None,
    factor=1.0):

  """
  Calculate VAR,LIQUID,VAPOR using linspace_bubble_(t or p), and use the values inside this function

  antoine: np.ndarray[[A1,B1,C1],[A2,B2,C2],...[AN,BN,CN]]

  """

  VAR=VAR/factor
  N=len(VAR)
  vz=np.zeros_like(VAR)
  vy=np.zeros_like(VAR)

  for i,vapor in enumerate(VAPOR):
    
    y1=vapor.composition()[0]
    vy[i]=y1
  for i,liquid in enumerate(LIQUID):
    z1=liquid.composition()[0]

    vz[i]=z1

  if (y_sup==None) and (y_inf==None):

    y_inf=np.min(VAR)*0.5
    y_sup=np.max(VAR)*1.5
  # print(y_inf,y_sup)

  xorv,orv,xbol,bol=exp_data
  # Tamanho e limites
  plt.figure(figsize=(x_figsize, y_figsize))
  plt.xlim(-0.01,1.01)
  plt.ylim(y_inf,y_sup)

  # Bolha
  plt.plot(vy,VAR,label=text,linestyle=bol_linestyle,color='black')
  plt.scatter(xorv,orv,marker='o',facecolors='none',edgecolors='black',)
  # Orvalho 
  plt.plot(vz, VAR,linestyle=orv_linestyle,color='black')

  plt.scatter(xbol,bol,marker='o',facecolors='none',edgecolors='black')
  plt.xlabel(x_label)
  plt.ylabel(y_label)
  plt.title(title)
  plt.legend()
  plt.grid(True)
  
  if save_fig:
    os.makedirs(plot_dir, exist_ok=True)
    filename = f"{title}.pdf"
    filepath = os.path.join(plot_dir, filename)
    plt.savefig(filepath,bbox_inches='tight')
  plt.show()


def linspace_bubble_p(eos,t,antoine,N=100):

  """
    return PRES,LIQUID,VAPOR

  """

  # Composiçaõ do comp1
  vz=np.linspace(0.00001,0.9999,N)
  PRES=np.zeros_like(vz)

  LIQUID=np.zeros_like(vz,dtype=object)
  # print(LIQUID)
  VAPOR=np.zeros_like(vz,dtype=object)
  # print(VAPOR)
  

  for (i,z1) in enumerate(vz):

    print(i)
    z=np.array([z1,1.0-z1])
    liquid_phase,vapor_phase=bubble_p(eos,t,z,antoine)
    PRES[i]=liquid_phase.pressure()
    LIQUID[i]=liquid_phase 
    VAPOR[i]=vapor_phase 


  return PRES,LIQUID,VAPOR

def linspace_bubble_t(eos,p,antoine,N=100,tol=1e-8,it_max=100):

  vz=np.linspace(0.00001,0.9999,N)
  TEMP=np.zeros_like(vz)
  LIQUID=np.zeros_like(vz,dtype=object)
  VAPOR=np.zeros_like(vz,dtype=object)

  for (i,z1) in enumerate(vz):

    z=np.array([z1,1.0-z1])

    liquid_phase,vapor_phase=bubble_t(eos,p,z,antoine,tol,it_max)
    TEMP[i]=liquid_phase.temperature()
    LIQUID[i]=liquid_phase 
    VAPOR[i]=vapor_phase 


  return TEMP,LIQUID,VAPOR





def bubble_p(eos,t,z,antoine,tol=1e-8,it_max=100):
  
  """
  Return np.ndarray[State(T,Pb,z),State(T,Pb,y)]
  """
  psat=psat_antoine(t,antoine)
  P0=z.dot(psat)
  yguess=(z*psat)/P0
  def F(P):
    
    # print(P)
    state_z= State.tpx(eos,t,P,z,"liquid")
    dg,state_x=state_z.min_tpd("vapor",yguess,tol,it_max)
    return dg,np.array([state_z,state_x])

    # F=lambda P: peq.tpd(T,P,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[0]
    # F_return_x=lambda X: peq.tpd(T,X,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[1]
  
  f=lambda x: F(x)[0]

  RESULT=opt.root(f,P0)
  Pbubble=RESULT.x
  _,states=F(Pbubble)
  return states

def orv_p(eos,t,z,antoine,tol=1e-8,it_max=100):
  
  """
  Return np.ndarray[State(T,Pb,z),State(T,Pb,y)]
  """
  psat=psat_antoine(t,antoine)
  P0=z.dot(psat)
  yguess=(z*psat)/P0
  def F(P):
    
    # print(P)
    state_z= State.tpx(eos,t,P,z,"vapor")
    dg,state_x=state_z.min_tpd("liquid",yguess,tol,it_max)
    return dg,np.array([state_z,state_x])

    # F=lambda P: peq.tpd(T,P,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[0]
    # F_return_x=lambda X: peq.tpd(T,X,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[1]
  
  
  f=lambda x: F(x)[0]

  RESULT=opt.root(f,P0)
  Porv=RESULT.x
  _,states=F(Porv)
  return states
def bubble_t(eos,p,z,antoine,tol=1e-8,it_max=100):
  
  """

  z: composition of developd phase
  antoine: np.ndarray[[A1,B1,C1],[A2,B2,C2],...[AN,BN,CN]]
  Return np.ndarray[State(Tb,P,z),State(Tb,P,y)]
  """
  tsat=tsat_antoine(p,antoine)
  T0=z.dot(tsat)
  yguess=(z*tsat)/T0

  def F(T):
    state_z= State.tpx(eos,T,p,z,"liquid")
    dg,state_x=state_z.min_tpd("vapor",yguess,tol,it_max)
    return dg,np.array([state_z,state_x])

    # F=lambda P: peq.tpd(T,P,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[0]
    # F_return_x=lambda X: peq.tpd(T,X,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[1]
  f=lambda x: F(x)[0]
  
  RESULT=opt.root(f,T0)
  Tbubble=RESULT.x
  _,states=F(Tbubble)
  return states

def tpd_root_(tp,
              var,
              z,
              incipient_phase,
              peq:PhaseEquilibrium,
              antoine,tol=1e-6):


  #Se fixo T, então é P-VLE
  if tp=="t":
    T=var
    psat=psat_antoine(T,antoine)

    if incipient_phase=="vapor":

      x0=z.dot(psat)
      incipient_phase_guess=(z*psat)/x0


    else:
      x0=(z.dot(1/psat))**-1
      incipient_phase_guess=x0*(z/psat)


    F=lambda X: peq.tpd(T,X,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[0]
    F_return_x=lambda X: peq.tpd(T,X,z,incipient_phase,incipient_phase_guess,tol=tol,it_max=100)[1]


  #Se fixo P, então é T-VLE

  elif tp=="p":
    P=var
    tsat=tsat_antoine(P,antoine)

    # print("P=",P)
    # print(tsat)
    if incipient_phase=="vapor":

      x0=z.dot(tsat)
      incipient_phase_guess=(z*tsat)/x0



    else:
      x0=(z.dot(1/tsat))**-1
      incipient_phase_guess=x0*(z/tsat)

    # print(x0,incipient_phase_guess)

    # !!!
    # incipient_phase_guess=np.array([0.5,0.5])

    F=lambda X: peq.tpd(X,P,z,incipient_phase,incipient_phase_guess,tol=1e-10,it_max=200)[0]
    F_return_x=lambda X: peq.tpd(X,P,z,incipient_phase,incipient_phase_guess,tol=1e-10,it_max=200)[1]


  else:
    return "t_ou_p"



  # print(psat)
  # print(x0,sum(incipient_phase_guess))

  # print(P,x0,incipient_phase_guess,z)
  # F=lambda X: peq.tpd(T,X,z,incipient_phase,incipient_phase_guess,tol=1e-6,it_max=100)
  # print('xguess=',incipient_phase_guess,'x0=',x0)
  # print("x0_guess")

  #!!!
  # x0=500
  # print('z=',z,'T0=',x0,"T=",result,"xGUESS=",incipient_phase_guess,"xResult=",x)

  result=opt.root(F,x0,method='lm' ).x
  x=F_return_x(result)

  return (result[0],x[0])
  # return (F,result[0],x[0],x0,incipient_phase_guess)

  # return F
def psat_antoine(T,v):

  psats=np.zeros(len(v))

  for (i,vec) in enumerate(v):
    a,b,c=vec[0],vec[1],vec[2]

    log=a-b/(T+c)

    psats[i]=(10**log)*1e5

  return psats

def tsat_antoine(P,v):

  tsats=np.zeros(len(v))

  for (i,vec) in enumerate(v):

    a,b,c=vec[0],vec[1],vec[2]

    tsats[i]=b/(a-np.log10(P*1e-5))-c

  return tsats