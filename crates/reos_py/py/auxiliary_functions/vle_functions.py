
from auxiliary_functions import *
import matplotlib.pyplot as plt
import scipy.optimize as opt
import os
def vle_diagram(T,phase_eq,N=100,factor=1e5):
  vx=np.linspace(0.00001,0.9999,N)
  calc_p=np.zeros_like(vx)
  calc_vy=np.zeros_like(vx)

  for i in range(len(vx)):

      try:
        x=np.array([vx[i],1-vx[i]])
        # p,y,*_=BubblePy(eos,T,x,guess_p,guess_y) #calculado em Pa
        p,y=phase_eq.bbpy(T,x,tol_p=1e-6,tol_y=1e-6)
        # print(p)
        calc_p[i]=p/factor #armazena em bar
        calc_vy[i]=y[0] #composição de acido acetico

      except Exception as e:
        print(e)
        continue

  return [calc_p,calc_vy,vx]

def get_bbp(EOS,T,antoine,factor=1e3,N_points=100):

  peq=PhaseEquilibrium(EOS)
  linspaceZ=np.linspace(0.00001,0.9999,N_points)
  PBOL=np.zeros_like(linspaceZ)
  XPHASE=np.zeros_like(linspaceZ)
  # XASCL=np.zeros_like(linspaceZ,dtype=object)
  # XASCV=np.zeros_like(linspaceZ,dtype=object)
  MAE=np.zeros_like(linspaceZ,dtype=object)
  FILHA=np.zeros_like(linspaceZ,dtype=object)

  tol=1e-6
  for (i,z1) in enumerate(linspaceZ):

    try:
      z=np.array([z1,1-z1])
      # print(z)

      # print(z)
      pcalc,w1=tpd_root_("t",T,z,"vapor",peq ,antoine,tol)

      MAE[i]=State.tpx(EOS,T,pcalc,z,"liquid")

      y=np.array( [w1,1-w1] )
      FILHA[i]=State.tpx(EOS,T,pcalc,y,"vapor")

      PBOL[i]=pcalc/factor
      XPHASE[i]=w1

    except Exception as e:
      print(e)

  return PBOL,linspaceZ,XPHASE,MAE,FILHA

def get_bbt(EOS,P,antoine,N_points=100):

  peq=PhaseEquilibrium(EOS)
  Z=np.linspace(0.00001,0.9999,N_points)

  TBOL=np.zeros_like(Z)
  XPHASE=np.zeros_like(Z)
  # XASCL=np.zeros_like(linspaceZ,dtype=object)
  # XASCV=np.zeros_like(linspaceZ,dtype=object)
  MAE=np.zeros_like(Z,dtype=object)
  FILHA=np.zeros_like(Z,dtype=object)

  tol=1e-6
  for (i,z1) in enumerate(Z):

    try:
      z=np.array([z1,1-z1])
      # print(z)

      # print(z)
      tcalc,w1=tpd_root_("p",P,z,"vapor",peq ,antoine,tol)

      state_mae=State.tpx(EOS,tcalc,P,z,"liquid")
      MAE[i]= 1.0/state_mae.volume()

      w=np.array( [w1,1-w1] )

      state_filha=State.tpx(EOS,tcalc,P,w,"vapor")
      FILHA[i]=1.0/state_filha.volume()

      TBOL[i]=tcalc
      XPHASE[i]=w1

    except Exception as e:
      print(e)
      continue

  return TBOL,Z,XPHASE,MAE,FILHA

def get_orvt(EOS,P,antoine,N_points=100):

  peq=PhaseEquilibrium(EOS)
  Z=np.linspace(0.00001,0.9999,N_points)

  TORV=np.zeros_like(Z)
  XPHASE=np.zeros_like(Z)
  # XASCL=np.zeros_like(linspaceZ,dtype=object)
  # XASCV=np.zeros_like(linspaceZ,dtype=object)
  MAE=np.zeros_like(Z,dtype=object)
  FILHA=np.zeros_like(Z,dtype=object)

  tol=1e-6
  for (i,z1) in enumerate(Z):

    try:
      z=np.array([z1,1-z1])
      # print(z)

      # print(z)
      tcalc,w1=tpd_root_("p",P,z,"liquid",peq ,antoine,tol)

      state_mae=State.tpx(EOS,tcalc,P,z,"vapor")
      rhoV=1.0/state_mae.volume()
      FILHA[i]=rhoV

      w=np.array( [w1,1-w1] )

      state_filha=State.tpx(EOS,tcalc,P,w,"liquid")
      rhoL=1.0/state_filha.volume()
      FILHA[i]=rhoL
      TORV[i]=tcalc
      XPHASE[i]=w1

    except Exception as e:
      print(e)
      continue
      
  return TORV,Z,XPHASE,MAE,FILHA

def calc_X(EOS,T,RHO,X1):


  XASC=np.zeros_like(RHO,dtype=object)
  association_residual=EOS.get_association()

  for (i,rho) in enumerate(RHO):

    try:


      z=np.array([X1[i],1-X1[i]])

      X=association_residual.non_bonded_sites(T[i],rho,z)
      XASC[i]=X*1   

    except Exception as e:
      print(e)

  return XASC

def VLE_DIAGRAM(p_or_t,
                eos,antoine,
                y_label,
                x_label,
                title,
                exp_data=None,
                factor=1.0,
                y_lim=None,
                x_figsize=5,
                y_figsize=5,
                save_fig=False,
                plot_dir="plots",
                tol=1e-6,
                N_points=100):

  '''
  P_or_T: tupla ; ex: ("t",300) em Kelvin; ("p",500) em bar
  xorv,orv,xbol,bol=exp_data
  '''

  var_str,var=p_or_t
  EOS=eos
  peq=PhaseEquilibrium(EOS)
  # print(var_str,var)
  linspaceZ=np.linspace(0.00001,0.9999,N_points)
  BOL=np.zeros_like(linspaceZ)
  XPHASE=np.zeros_like(linspaceZ)
  XASCL=np.zeros_like(linspaceZ,dtype=object)
  XASCV=np.zeros_like(linspaceZ,dtype=object)
  
  for (i,z1) in enumerate(linspaceZ):

    try:
      z=np.array([z1,1-z1])
      # print(z)

      # print(z)
      calc_bol,w1=tpd_root_(var_str,var,z,"vapor",peq ,antoine,tol)

      #novo calculo
      #dado T,P,x, calcular XASC
      t=var
      p=calc_bol
      state_mae=State.tpx(EOS,t,p,z,"liquid")
      rho=1.0/state_mae.volume()
      association_residual=EOS.get_association()
      X=association_residual.non_bonded_sites(t,rho,z)
      XASCL[i]=X*1

      y=np.array( [w1,1-w1] )
      state_fil=State.tpx(EOS,t,p,y,"vapor")
      rho=1.0/state_fil.volume()
      association_residual=EOS.get_association()
      X=association_residual.non_bonded_sites(t,rho,y)
      XASCV[i]=X*1      

      # calc_bol,x1=tpd_root_(var_str,var,z,"vapor",peq ,antoine)
      # calc_orv=tpd_root_(var_str,var,z,"liquid",peq,antoine)
      # print('z=',z,'var=',var)
      # print('Tbol=',calc_bol)
      # print('x1=',x1)
      BOL[i]=calc_bol/factor
      XPHASE[i]=w1

    except Exception as e:
      print(e)
    # ORV[i]=calc_orv/factor

  if y_lim==None:
    y_sup=np.max(BOL)*1.1
    y_inf=np.min(BOL)*0.9
  else:
    y_sup=y_lim[1]
    y_inf=y_lim[0]

  # ORV
  # plt.scatter(XPHASE,BOL)
  if exp_data!=None:
    plt.figure(figsize=(x_figsize, y_figsize))
    plt.xlim(-0.01,1.01)
    plt.ylim(y_inf,y_sup)
    xorv,orv,xbol,bol=exp_data

    plt.plot(XPHASE,BOL,color='black')
    plt.scatter(xorv,orv,color='black')

    plt.plot(linspaceZ, BOL,color='black')
    plt.scatter(xbol,bol,color='black')

    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    if save_fig:
      os.makedirs(plot_dir, exist_ok=True)

      filename = f"{title}.png"
      filepath = os.path.join(plot_dir, filename)
      plt.savefig(filepath)

    plt.show()
  else:
    plt.figure(figsize=(x_figsize, y_figsize))
    plt.xlim(-0.01,1.01)
    plt.ylim(y_inf,y_sup)
    plt.plot(XPHASE,BOL,color='black')
    # BOL
    # plt.scatter(linspaceZ, BOL)
    plt.plot(linspaceZ, BOL,color='black')


    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.legend()
    plt.grid(True)
    if save_fig:
      os.makedirs(plot_dir, exist_ok=True)

      filename = f"{title}.png"
      filepath = os.path.join(plot_dir, filename)
      plt.savefig(filepath)
      
  
  

  return BOL,XPHASE,linspaceZ,XASCL,XASCV

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