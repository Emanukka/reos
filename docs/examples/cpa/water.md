---
hide:
  - toc
---

# Water

<!-- The correspondent python code can be found at [pure_water.py](../../../python_examples/phase_equilibrium/pure_water.py). -->

``` py title = "Initializing SCPA EOS"
import numpy as np

from reos.state import State
from reos.consts import Consts
from reos.cpa import CPAParameters
from reos.eos import EquationOfState


parameters = CPAParameters.from_json(["water"], "/parameters/cpa/kontogeorgis2006.json")
eos = EquationOfState.scpa(parameters)

print(parameters) 
```

``` sh title = "CPA Parameters"
CPAParameters(
  CubicParameters(
  a0=[0.12277],
  b=[0.000014515],
  c1=[0.67359],
  tc=[647.29]),
  AssociativeParameters(
  na=1, nb=1, nc=0,
  sites=[
     Site(type=A,owner=0,idx=0,mul=2,eps=16655,kappa=0.0692),
     Site(type=B,owner=0,idx=1,mul=2,eps=16655,kappa=0.0692)],
  interactions=[
     SiteInteraction(j=0,l=1,epsilon=16655,kappa=0.0692, rule='cr1')])
)
```

``` py title = "Psat algorithm "
def calc_psat(t, p0):
    
    e = 1
    it = 0

    while abs(e) > 1e-8 and it < 100:

        s1 = State.tpx(eos, t, p0, np.array([1.0]), 'vapor')
        s2 = State.tpx(eos, t, p0, np.array([1.0]), 'liquid')
        phiv = np.exp(s1.lnphi()[0])
        phil = np.exp(s2.lnphi()[0])
        r = phil / phiv
        e = 1.0 - r
        p0 = p0 * r
        it += 1

    return p0, s1, s2
```

``` py title = "Computing properties"
N = 100
T=np.linspace(250.0, 650, N)

PRES = np.zeros_like(T)
VAP = np.zeros_like(T,dtype=object)
LIQ = np.zeros_like(T,dtype=object)

XV = np.zeros_like(VAP,dtype=object)
XL = np.zeros_like(LIQ,dtype=object)

rhoV = np.zeros_like(VAP)
rhoL = np.zeros_like(LIQ)

entropyV = np.zeros_like(VAP)
entropyL = np.zeros_like(LIQ)

for (i,t) in enumerate(T):
    
    PRES[i], VAP[i], LIQ[i] = calc_psat(t, PRES[i-1] if i>0 else 1e5)
    
    liq = LIQ[i]
    vap = VAP[i]

    entropyL[i] = liq.entropy() 
    entropyV[i] = vap.entropy() 

    rhoV[i] = vap.density
    rhoL[i] = liq.density

DS = entropyV - entropyL
DH = T * DS
```

<figure markdown="span">
  ![ ](pdee.png){ width="800" }
  <figcaption></figcaption>
</figure>

``` py title = "Non-associated sites fraction"
for (i,t) in enumerate(T):
    
    ...

    XL[i] = eos.unbonded_sites_fraction(t, rhoL[i], np.array([1.0]))[0]
    XV[i] = eos.unbonded_sites_fraction(t, rhoV[i], np.array([1.0]))[0]
```

<figure markdown="span">
  ![ ](xsat.png){ width="400" }
  <figcaption></figcaption>
</figure>
<!-- ![alt text](./pdes.png) -->