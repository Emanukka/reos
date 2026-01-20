#%%
import numpy as np
import pandas as pd
import tabula
from math import isclose


from reos.cpa import CPAPureRecord,CPABinaryRecord, CPAParameters
from reos.eos import EquationOfState
# from reos.consts import Consts
# from reos.cubic import CubicBinaryRecord
# R = Consts.ideal_gas_const()
substances = pd.read_json("../../substances.json")

#%%
csv_kontogeorgis2010 = pd.read_csv("./kontogeorgis2010.csv")
csv_kontogeorgis2006 = pd.read_csv("./kontogeorgis2006.csv")
csv_tsvintzelis2011 = pd.read_csv("./tsivintzelis2011.csv")

#%% kontogeorgis2010 

records = []

for index, row in csv_kontogeorgis2010.iterrows():

    name = row.compound.lower()
    record = CPAPureRecord.new(
        name, 
        substances[name]["molar_weight"], 
        a0 = row.a0 / 10, 
        b = row.b / 1000, 
        c1 = row.c1, 
        tc = row.tc, 
        epsilon = row.epsilon * 100, 
        kappa = row.kappa / 1000, 
        na = row.na, 
        nb = row.nb,
        nc = row.nc)
    records.append(record)


json = CPAPureRecord.to_json_vec("../kontogeorgis2010", records, build = True)
 
#%% kontogeorgis2006

records = []

for index, row in csv_kontogeorgis2006.iterrows():

    name = row.compound.lower()
    record = CPAPureRecord.new(
        name, 
        substances[name]["molar_weight"], 
        a0 = row.a0, 
        b = row.b, 
        c1 = row.c1, 
        tc = row.tc, 
        epsilon = row.epsilon,
        kappa = row.kappa, 
        na = row.na, 
        nb = row.nb,
        nc = row.nc)
    records.append(record)


json = CPAPureRecord.to_json_vec("../kontogeorgis2006", records, build = True)

# todo: from_multiple_jsons
# from_multiple_jsons(names = [names1[str],names2[str],...], ppaths = [json1,json2], bpath = json)

#%% tsivintzelis2011

records = []

for index, row in csv_tsvintzelis2011.iterrows():

    name = row.compound.lower()
    record = CPAPureRecord.new(
        name, 
        substances[name]["molar_weight"], 
        a0 = row.a0 / 10, 
        b = row.b / 1000, 
        c1 = row.c1, 
        tc = row.tc, 
        epsilon = row.epsilon * 100, 
        kappa = row.kappa / 1000, 
        na = row.na, 
        nb = row.nb,
        nc = row.nc)
    records.append(record)


json = CPAPureRecord.to_json_vec("../tsivintzelis2011", records, build = True)
# %%
