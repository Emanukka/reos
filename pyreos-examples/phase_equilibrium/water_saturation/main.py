#%%
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from os import listdir
from reos.eos import EquationOfState
from reos.state import State
from reos.cpa import CPAParameters
from si_units import KELVIN, BAR, PASCAL, CELSIUS
from calculators import *

xsize = 3.15
ysize = 3.15

SAVE = True
PLOT = True
PLTDIR = "./plots/"
DATAPATH = "./data/"
HEAD = f"{'case':<{20}}, {'t (ms)':<10}, {'v (points/s)':<15}"


plt.rcParams.update({
    "text.usetex": True,
    "font.serif": ["Computer Modern"], 
    "axes.labelsize": 14,
    "font.size": 14,
    "legend.fontsize": 14,
    "xtick.labelsize": 13.3,
    "ytick.labelsize": 13.3,
    "figure.figsize": (xsize * 1.5, ysize),  

})

pw = CPAParameters.from_json(["water"],"../../../parameters/cpa/tsivintzelis2011.json")
eosw = EquationOfState.scpa(pw)

def timing(f):
    
    def wrap(*args, **kwargs):

        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        dt = (time2-time1)
        v = 100 / dt
        case_name = ret[-1]
        info = f"{case_name:<20}, {dt*1000:<10.4f}, {v:<15.4f}"

        print(HEAD)
        print(info)
        # return ret, info
        return ret
    
    return wrap

@timing
def compute(
        case, 
        ppath = None,  
        bpath = None,
        p0 = 1 * BAR,
        pf = 600 * BAR, 
        N = 100):
    
    if ppath is None: ppath = "../../../parameters/cpa/tsivintzelis2011.json"
    if bpath is None: bpath = "../../../parameters/cpa/tsivintzelis2011_binary.json"

    case_name = "/".join(case["names"])
    names = ["ethane"] + case["names"]
    y_dry_gas = [0.0] + case["y_dry_gas"]
    t = case["t"] / KELVIN
    y_dry_gas = np.array(y_dry_gas)


    p = CPAParameters.from_json(names, ppath, bpath)

    eos = EquationOfState.scpa(p)
    vpressure, yw, vdeveloped, vincipient = linspace_wsat(
        eos, 
        eosw, 
        t, 
        y_dry_gas, 
        p0 / PASCAL, 
        pf / PASCAL, 
        N)
    
    return vpressure, yw, vdeveloped, vincipient, case_name


#%%

def nf(filename):

    casename = "".join(filename.split(".")[0])
    content = casename.split("@")

    old_fracs = content[1].split("#")
    old_names = content[0].split("#")
    
    assert(len(old_fracs) == len(old_names))

    new_names = []
    new_fracs = []
    
    n = len(old_fracs) 

    for i in range(n):

        old_name = old_names[i]
        old_frac = old_fracs[i]

        new_names.append(" ".join(old_name.split("_")))
        new_fracs.append(float(".".join(old_frac.split("_"))))

    return new_names, new_fracs, casename




#%%


filenames = listdir(DATAPATH)
casenames = [] 
pressure = []
water_content = []
nfiles = 0

for filename in filenames:

    n += 1
    excel = pd.ExcelFile(DATAPATH + filename)

    sheet_names = excel.sheet_names

    names, fracs, casename = nf(filename)

    for j, t in enumerate(sheet_names):

        casenames.append(casename + '@' + '_'.join(t.split('.') ) )      

        t = float(t) * KELVIN
        case = {"names": names, "y_dry_gas": fracs, "t": t}
        vpressure, yw, vdeveloped, vincipient, _ = compute(case, N = 100)

        pressure.append(vpressure * PASCAL / BAR)
        water_content.append(yw * 1e6)


#%%

mk = ["o","s","^",">"]

for (i, t) in enumerate(sheet_names):

    plt.figure(figsize=(xsize * 1.5, ysize))
    df = pd.read_excel(excel, sheet_name = t)

    plt.scatter(df["p"], df["y"] * 1e6, marker=mk[i], label = f"{t} K", facecolors="None", edgecolors="black")
    plt.plot(pressure[i], water_content[i], "black")

    plt.xlabel(r"$ \mathrm{P / bar}$")
    plt.ylabel(r"$ \mathrm{Water \ content / ppm}$")
    plt.ylim(0, 15_000)
    plt.xlim(0, 600)
    plt.tight_layout()
    plt.legend()

    if SAVE:

        casename = casenames[i]
        plt.savefig(PLTDIR + casename + '.png', dpi = 300)


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


