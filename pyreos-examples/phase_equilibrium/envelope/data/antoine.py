#%%
import numpy as np
import pandas as pd
# from reos.reos  import EquationOfState,State,CPAParameters,PhaseEquilibrium,CubicRecord,AssociationRecord



# c_co2=CubicRecord(
#     a0=0.35079,
#     b=0.0272e-3,
#     c1=0.7602,
#     tc=304.12,
# )
# c_ch4=CubicRecord(
#     a0=0.23204,
#     b=0.0291e-3,
#     c1=0.447,
#     tc= 190.56,
# )
# c_w=CubicRecord(
#     a0= 0.12277,
#     b=0.0145e-3,
#     c1=0.6736,
#     tc=647.096,
# )
# c_w3b=CubicRecord(
#     a0= 3.005960e-1,
#     b=0.014969e-3,
#     c1=0.35928,
#     tc=647.096,
# )
# a_w3b=AssociationRecord.associative(
#     eps=207.97e2,
#     beta=21.3e-3,
#     b=0.014969e-3,
#     na=1,
#     nb=2,
#     nc=0)

# c_acoh=CubicRecord(
#     a0= 0.91196,
#     b=0.0468e-3,
#     c1=0.4644,
#     tc=594.8,
# )
# # antes
# c_propanoic=CubicRecord(
#     a0= 1.32676,
#     b=0.06406e-3,
#     c1= 0.6891,
#     tc= 607.0,
# )
# c_methanol_2b=CubicRecord(
#     a0= 0.40531,
#     b=0.0000309,
#     c1=0.4310,
#     tc=513.,
# )
# c_methanol_3b=CubicRecord(
#     a0= 4.5897e-1,
#     b= 0.0334e-3,
#     c1=1.0068,
#     tc=513.,
# )

# c_octanol_2b=CubicRecord(
#     a0= 4.15822,
#     b=0.1485e-3,
#     c1=1.1486,
#     tc=655.5,
# )
# c_octanol_3b=CubicRecord(
#     a0= 41.9005e-1,
#     b= 0.1489e-3,
#     c1= 1.0550,
#     tc=655.5,
# )
# c_octane=CubicRecord(
#     a0= 34.8750e-1,
#     b=0.1424e-3,
#     c1= 0.99415,
#     tc=568.7,
# )
# c_heptane=CubicRecord(
#     a0= 29.17800e-1,
#     b=0.125350e-3,
#     c1= 0.913700,
#     tc=540.0,
# )


# c_decane=CubicRecord(
#     a0= 4.7389,
#     b=0.17865e-3,
#     c1= 1.13243,
#     tc=617.7,
# )
# a_decane=AssociationRecord.inert(b=0.17865e-3)

# c_h2s=CubicRecord(
#     a0= 4.45050e-1,
#     b=0.0285e-3,
#     c1= 0.60265,
#     tc=373.3,
# )
# a_co2=AssociationRecord.solvate(
#     b=0.0272e-3,
#     na=0,
#     nb=1,
#     nc=0)

# a_w=AssociationRecord.associative(
#     eps=166.55e2,
#     beta=0.0692,
#     b=0.0145e-3,
#     na=2,
#     nb=2,
#     nc=0)

# a_acoh=AssociationRecord.associative(
#     eps=403.23e2,
#     beta=4.5e-3,
#     b=0.0468e-3,
#     na=0,
#     nb=0,
#     nc=1)

# a_methanol_2b=AssociationRecord.associative(
#     eps=24591.0,
#     beta=0.01610,
#     b=0.0000309,
#     na=1,
#     nb=1,
#     nc=0)

# a_methanol_3b=AssociationRecord.associative(
#     eps=160.70e2,
#     beta=34.4e-3,
#     b= 0.0334e-3,
#     na=2,
#     nb=1,
#     nc=0)

# a_propanoic=AssociationRecord.associative(
#     eps=399.75e2,
#     beta=0.0021,
#     b= 0.0641e-3,
#     na=0,
#     nb=0,
#     nc=1)

# a_octanol_2b=AssociationRecord.associative(
#     eps=267.59e2,
#     beta=0.14e-3,
#     b= 0.1485e-3,
#     na=1,
#     nb=1,
#     nc=0)
# a_octanol_3b=AssociationRecord.associative(
#     eps=250.00e2,
#     beta= 0.2e-3,
#     b= 0.1489e-3,
#     na=2,
#     nb=1,
#     nc=0)

# #  0.0500 8.5755 1.0564 150.00 17.3
# c_ethanol_3b=CubicRecord(
#     a0= 8.57550e-1,
#     b=0.0500e-3,
#     c1= 1.0564,
#     tc=514.0,
# )

# a_ethanol_3b=AssociationRecord.associative(
#     eps=150.00e2,
#     beta= 17.3e-3,
#     b= 0.0500e-3,
#     na=2,
#     nb=1,
#     nc=0)

# a_octane=AssociationRecord.inert(0.1424e-3)
# a_heptane=AssociationRecord.inert(0.125350e-3)
# a_ch4=AssociationRecord.inert(0.0291e-3)

# a_h2s=AssociationRecord.solvate(
#     b=0.0285e-3,
#     na=0,
#     nb=2,
#     nc=0)

# c_mea=CubicRecord(
#     a0= 14.112e-1,
#     b= 0.05656e-3,
#     c1= 0.7012,
#     tc=670.0,
# )
# a_mea=AssociationRecord.associative(
#     eps=181.77e2,
#     beta= 5.35e-3,
#     b= 0.05656e-3,
#     na=2,
#     nb=2,
#     nc=0)



# tudo do NIST
# log10, P/BAR,T/Kelvin

water_antoine = np.array([6.20963,2354.731,  7.559])
acoh_antoine = np.array([4.68206, 1642.54,    -39.764 ])
co2_antoine = np.array([6.81228,  1301.679,   -3.494])
octane_antoine = np.array([4.04867    ,1355.126   ,-63.633    ])
propanoic_antoine = np.array([4.74558 ,1679.869,  -59.832])

heptane_antoine = np.array([4.02832   ,1268.636   ,-56.199])

octanol_antoine = np.array([6.47682   ,2603.359   ,-48.799    ])
metoh_antoine = np.array([5.20409,    1581.341,   -33.50])

h2s_antoine = np.array([4.52887	,958.587	,-0.539	])

mea_antoine = np.array([4.29252	,1408.873	,-116.093])

ethanol_antoine = np.array([5.24677,1598.673,-46.424])


utilis = {
    "water":water_antoine,
    "acetic acid":acoh_antoine,
    "carbon dioxide":co2_antoine,
    "n-octane":octane_antoine,
    "propionic acid":propanoic_antoine,
    "n-heptane":heptane_antoine,
    "1-octanol":octanol_antoine,
    "methanol":metoh_antoine,
    "hydrogen sulfide":h2s_antoine,
    "mea":mea_antoine,
    "ethanol":ethanol_antoine
}


# __all__=[c_co2]

# %%

df = pd.DataFrame(utilis)

df.index = ["A","B","C"]# df.columns = ["A","B","C"]

with pd.ExcelWriter('data.xlsx',mode='a',if_sheet_exists='replace') as writer:  # doctest: +SKIP
    
    df.to_excel(writer, sheet_name='antoine')

# df.to_excel("data.xlsx",sheet_name="antoine")
#%%
# df.rename(columns={"0":"A", "1":"B", "2":"C"})
