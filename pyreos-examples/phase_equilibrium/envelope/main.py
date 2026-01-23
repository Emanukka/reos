#%%

import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from reos.eos import EquationOfState
from reos.state import State
from reos.cpa import CPAParameters
from si_units import KELVIN, BAR, PASCAL, KILO
from vle_functions import *

xsize = 3.15
ysize = 3.15

SAVE = True
PLOT = True
PLTDIR = "plots"
DATAPATH = "data/data.xlsx"
HEAD = f"{'case':<{20}}, {'t (ms)':<10}, {'v (points/s)':<15}"

plt.rcParams.update({
    "text.usetex": True,
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "figure.figsize": (xsize, ysize),  

})

df_antoine = pd.read_excel(DATAPATH,sheet_name="antoine")

def run_bubble_t(name1, name2, P,
                 ppath = "../../../parameters/cpa/kontogeorgis2010.json", 
                 bpath = "../../../parameters/cpa/kontogeorgis2010_binary.json"):
    
    p = CPAParameters.from_json(
        [name1, name2],
        ppath,
        bpath)

    eos = EquationOfState.scpa(p)

    antoine1 = df_antoine[name1].to_numpy()
    antoine2 = df_antoine[name2].to_numpy()
    antoine = np.array([antoine1, antoine2])


    T, LIQUID, VAPOR = linspace_bubble_t(eos, P, antoine, N = 100)
    
    return T, LIQUID, VAPOR


def timing(f):
    
    def wrap(*args, **kwargs):

        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        dt = (time2-time1)
        v = 100 / dt
        case_name = ret[-1]
        info = f"{case_name:<20}, {dt*1000:<10.4f}, {v:<15.4f}"

        # print(HEAD)
        # print(info)
        return ret, info
    
    return wrap

@timing
def compute(case:dict):
    
    ppath = case.get("ppath", "../../../parameters/cpa/kontogeorgis2010.json")
    bpath = case.get("bpath", "../../../parameters/cpa/kontogeorgis2010_binary.json")

    name1 = case["id1"]
    name2 = case["id2"]
    
    p = CPAParameters.from_json(
        [name1, name2],
        ppath,
        bpath)

    eos = EquationOfState.scpa(p)

    var = case.get("var")
    
    if var is None : raise ValueError("suply 'var' ")

    antoine1 = df_antoine[name1].to_numpy()
    antoine2 = df_antoine[name2].to_numpy()
    antoine = np.array([antoine1, antoine2])

    try:
        if var.has_unit(PASCAL):
            
            p = var / PASCAL
            vtarget, vliq, vvap = linspace_bubble_t(eos, p, antoine, N=100)
        
        elif var.has_unit(KELVIN):
            t = var / KELVIN
            # print(t)
            vtarget, vliq, vvap = linspace_bubble_p(eos, t, antoine, N=100)
        
        case_name = "/".join([name1,name2])

        return vtarget, vliq, vvap, case_name 
    
    except Exception as e:
        print(e)
        return e


def unit_from_str(s:str):

    s = s.lower()

    if s == str(BAR).lower(): # '100 kpa'
        return BAR
    elif s == str(KILO * PASCAL).lower(): # '1 kpa'
        return KILO * PASCAL
    
    elif s == str(KELVIN).lower(): # '1 k'
        return KELVIN
    else:
        raise ValueError("unit didnt match")

def plot(
        case:dict, 
        vtarget,
        vliq,
        vvap, 
        save = False, 
        pltdir = "plots",
        xsize = 3.15,
        ysize = 3.15):
    

    var = case["var"]
    name1 = case["id1"]
    name2 = case["id2"]

    if var.has_unit(PASCAL):
        target_str = "t"
        unit_str = "K"
        cf = 1.0

    elif var.has_unit(KELVIN):
        target_str = "p"    
        unit_str = "bar"
        cf = 1e-5

    n = len(vtarget)
    vz = np.zeros(n)
    vy = np.zeros(n)

    for i,vapor in enumerate(vvap):
        
        y1 = vapor.composition[0]
        vy[i] = y1

    for i,liquid in enumerate(vliq):

        z1 = liquid.composition[0]
        vz[i] = z1

    plt.figure(figsize = (xsize, ysize))

    try:

        df = pd.read_excel(DATAPATH,sheet_name=",".join([name1,name2]))

        target_unit = unit_from_str(df["unit"][0])

        if target_unit.has_unit(PASCAL): 

            target_str = "p"
            cf = PASCAL / target_unit 

        elif target_unit.has_unit(KELVIN): 

            target_str = "t"
            cf = 1.0

        unit_str = str(target_unit)

        if unit_str.split()[0] == "1":
            unit_str = unit_str.split()[1]
        else:
            unit_str = unit_str

        exp_data = [df["x"],df[target_str + "x"],df["y"],df[target_str + "y"]]
        x, target_x, y, target_y = exp_data

        plt.scatter(x, target_x, marker='o', facecolors = 'none', edgecolors = 'black')
        plt.scatter(y, target_y, marker='o', facecolors = 'none', edgecolors = 'black')
    
    except Exception as e:
        print(e)

    plt.xlim(-0.01,1.01)

    vtarget *= cf

    plt.plot(vy, vtarget, color = 'black')
    plt.plot(vz, vtarget, color = 'black')

    plt.xlabel(r"$x_1,y_1$")


    plt.ylabel(fr"${target_str.upper()}/{unit_str}$")
    plt.grid(True)
    
    if save:
        name1 = "_".join(name1.split())
        name2 = "_".join(name2.split())
        prefx = name1 + "@" + name2 + "@" + "".join(str(var).split())
        os.makedirs(pltdir, exist_ok = True)
        filename = f"{prefx}.png"
        filepath = os.path.join(pltdir, filename)
        plt.savefig(filepath, bbox_inches='tight', dpi=300)

    # plt.show()
    return True


#%%Cases

cases = [
    {"id1":"propionic acid", "id2":"n-heptane", "var":101.33e3 * PASCAL,},

    {"id1":"water",          "id2":"acetic acid","var":313.15 * KELVIN,
        "ppath": "../../../parameters/cpa/kontogeorgis2006.json",
        "bpath": "../../../parameters/cpa/kontogeorgis2006_binary.json"},

    {"id1":"methanol", "id2":"1-octanol", "var":101.32e3 * PASCAL},

    {"id1":"acetic acid",    "id2":"n-octane",   "var":343.2 * KELVIN},

    {"id1":"ethanol",    "id2":"water",   "var":298.14 * KELVIN,
        "ppath": "../../../parameters/cpa/kontogeorgis2006.json",
        "bpath": "../../../parameters/cpa/kontogeorgis2006_binary.json"},
    
    ]

infos = [HEAD]
for case in cases:

    result, info = compute(case)
    vtarget, vliq, vvap, _ = result
    infos.append(info)
    plot(case, vtarget, vliq, vvap, save = False)

data = '\n'.join(infos)

print(data)

# #%%Test
# test_case =   {"id1":"methane",    "id2":"carbon dioxide",   "var": 240 * KELVIN,
#         "ppath": "../../../parameters/cpa/tsivintzelis2011.json",
#         "bpath": "../../../parameters/cpa/tsivintzelis2011_binary.json"}

# result, info = compute(test_case)
# vtarget, vliq, vvap, _ = result
# # infos.append(info)
# plot(case, vtarget, vliq, vvap, save = False)

#%%
try:

    from cpuinfo import get_cpu_info

    cpuinfo = get_cpu_info()
    filename = "bench/" + "_".join(cpuinfo["brand_raw"].split()) + ".csv"

    with open(filename, "x") as f:
        f.write(data)

except Exception as e:
    print(e)


