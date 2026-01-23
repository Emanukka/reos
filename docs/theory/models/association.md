# Association

The implementation of association model in ℜ*eos* is based on the work of [Michelsen 2001](https://doi.org/10.1016/S0378-3812(01)00344-2) and [Michelsen 2006](https://doi.org/10.1021/ie060029x). The dimensionless association contribution to the Helmholtz free energy is given by:

$$
a^{asc} = \sum_j^S m_j \left( \ln X_j - \frac{X_j}{2} + \frac{1}{2} \right) M_j
$$

where $S$ is the number of distinguishable sites in the mixture, $m_j$ is the molar fraction of site $j$ and $M_j$ is the multiplicity of site $j$.

The fraction of non-bonded sites $X_j$ is obtained by solving implicitily:

$$
X_j = \frac{m_j}{m_j +  \sum_l^S K_{jl} X_l M_l}
$$
where $K_{jl}$ is the association constant of sites $j$ and $l$: 
$$
K_{jl} = \rho m_j m_l \Delta_{jl}
$$
and $\Delta^{jl}$ is the association strength between sites $j$ and $l$. Because the expressions for $\Delta^{jl}$ depend on the specific association model, we can use the general expression of the dimensionless association strength:
$$
\delta^{jl} = \kappa^{jl} \left[ \exp\left( \frac{\epsilon^{jl}}{RT} \right) - 1 \right]
$$
where $\kappa^{jl}$ is the association volume and $\epsilon^{jl}$ is the association energy between sites $j$ and $l$.
Then,

$$
\Delta^{jl} = f_v ^{ik} \delta^{jl}
$$
where $f_v^{ik}$ is a volumectric factor that depends on the specific model and the the combining rule choosen for the components hosts $i$ and $k$ of sites $j$ and $l$, respectively.

<div align="center">

<table>
  <thead>
    <tr>
      <th>Model</th>
      <th>Volumetric factor \(f_v^{ik}\)</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td align="center">PC-SAFT</td>
      <td align="center">\(d_{ik}^3\, g_{ik}^{seg}\)</td>
    </tr>
    <tr>
      <td align="center">CPA</td>
      <td align="center">\(b_{ik}\, g(\rho)\)</td>
    </tr>
  </tbody>
</table>

</div>

## Association sites

Each site $s_j$ has a owner $i$, a type $\alpha$, an index $j$, a multiplicity $M_j$ and the parameters $\epsilon^{j}$ and $\kappa^{j}$.
The type of each site defines the allowed interactions with other sites.

Then, to satisfy all association schemes proposed by [Huang and Radosz 1990](https://pubs.acs.org/doi/10.1021/ie00107a014), we propose 3 types of sites:
$A, B$ and $C$. $A$ and $B$ sites can interact with all types except with sites of the same type, while $C$ sites can interact with all types, including themselves.
<!-- ![alt text](grid.png) -->

<figure markdown="span">
  ![ ](grid.png){ width="600" }
  <figcaption></figcaption>
</figure>

``` py title="Example: Association parameters with pure water 4C in Reos"
AssociativeParameters(
na=1, nb=1, nc=0,
sites=[
   Site(type=A,owner=0,idx=0,mul=2,eps=16655,kappa=0.0692),
   Site(type=B,owner=0,idx=1,mul=2,eps=16655,kappa=0.0692)],
interactions=[
   SiteInteraction(j=0,l=1,epsilon=16655,kappa=0.0692, rule='cr1')
   ])
```

## Calculation of unbonded site fractions

To calculate the fractions of unbonded sites we can use the successive substitution method or the second-order method proposed by [Michelsen 2006](https://doi.org/10.1021/ie060029x):

$$
g_j = \frac{\partial Q}{\partial X_j}
= m_j \left( \frac{1}{X_j} - 1 \right)
- \sum_{l=1}^{S} K_{jl} X_l M_l
$$

$$
\hat{H}_{jl} = -\delta_{jl}\left( \frac{m_{j}}{X_{j}^2} \right) - K_{jl}
$$

Otherwise, in simple cases, the fractions of unbonded sites can be calculated analytically.

If the mixture contains only 1 site of type $A$ and 1 site of type $B$, the fractions of unbonded sites can be calculated analytically:

$$
X_2 = \frac{2 m_1 m_2 }{m_{1} m_2 + k\left(m_1 M_1 - m_2 M_2\right) + \sqrt{4km_2m_1^2M_1 + \left[m_1m_2 + k \left(m_2 M_2 m_1 M_1 \right) \right]^2}} 
$$

$$
X_1 = \frac{m_1}{m_1 + k X_2 M_2}
$$
where $k = K_{12} = K_{21}$.

If the mixture contains only 1 site of type $C$, the fraction of unbonded fraction is given by:

$$
X_1 = \frac{2 m_1}{m_1 + \sqrt{4 k M_1 m_1 + m_1^2}}
$$

where $k = K_{11}$.

<!-- ## CPA expressions  -->
<!-- 
$$
\frac{P^{asc}}{RT} = -\rho \left(1 + \rho \frac{\partial \ln g}{\partial \rho} \right)
$$

$$
\frac{\mu_i^{asc}}{RT} = \sum_{j \ \mid \ c\left(j\right) = i} \ln X_j M_j - \frac{1}{2}  \left( n \frac{\partial \ln g}{\partial n_i} \right) h
$$

$$
s^{asc}= - a^{asc} -  \frac{1}{2RT} \sum_{j}^S \sum_{l}^S K_{jl} \epsilon_{jl}X_j X_l M_j M_l
$$ -->
