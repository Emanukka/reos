#%%
import numpy as np
import pandas as pd
from math import isclose


from reos.cubic import CubicPureRecord,CubicBinaryRecord, CubicParameters
from reos.eos import EquationOfState
# from reos.consts import Consts
# from reos.cubic import CubicBinaryRecord
# R = Consts.ideal_gas_const()
substances = pd.read_json("../../substances.json")

#%%
csv_jaubert2016 = pd.read_csv("./jaubert2016.csv")

#%% jaubert2016 

records = []

for index, row in csv_jaubert2016.iterrows():

    name = row.Name.lower()
    # print(name)
    try:
        molar_weight = substances[name]["molar_weight"]
    except KeyError:
        print(f"Substance {name} not found in substances.json")
        molar_weight = 0.0
    
    if row.c =="N.A.":
        print(f"Substance {name} has c = N.A., skipping")
        continue

    record = CubicPureRecord.new(
        name=name, 
        molar_weight=molar_weight, 
        tc = row.tc,
        pc = row.pc * 1e5,
        l = row.L,
        m = row.M,
        n = row.N,
        c = float(row.c) * 1e-6
        )
    
    records.append(record)

#%%
# CubicPureRecord.new(
        # name=name, 
        # tc = row.tc,
        # pc = row.pc * 1e5,
        # w = row.w,
        # l = row.L,
        # m = row.M,
        # n = row.N,
        # c = float(row.c)
        # )
json = CubicParameters.to_json_vec("../jaubert2016", records, build = True)
 
# %%
