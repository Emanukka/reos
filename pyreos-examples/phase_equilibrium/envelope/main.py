#%%
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import listdir

from reos.eos import EquationOfState
from reos.state import State
from reos.cpa import CPAParameters
from si_units import KELVIN, BAR, PASCAL, KILO
from vle_functions import *

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
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "figure.figsize": (xsize, ysize),  

})

antoine_parameters = pd.read_csv("./antoine.csv")

#%%
def timing(f):
    
    def wrap(*args, **kwargs):

        time1 = time.time()
        ret = f(*args, **kwargs)
        time2 = time.time()
        dt = (time2-time1)

        # info = f"{case_name:<20}, {dt*1000:<10.4f}, {v:<15.4f}"
        # print(info)
        return ret, dt
    
    return wrap

@timing
def compute(case:dict):
    
    ppath = case.get("ppath", "../../../parameters/cpa/kontogeorgis2010.json")
    bpath = case.get("bpath", "../../../parameters/cpa/kontogeorgis2010_binary.json")

    name1 = case["id1"]
    name2 = case["id2"]
    
    parameters = CPAParameters.from_json([name1, name2], ppath, bpath)

    eos = EquationOfState.scpa(parameters)

    fix = case["fix"]
    
    antoine1 = antoine_parameters[name1].to_numpy()
    antoine2 = antoine_parameters[name2].to_numpy()
    antoine = np.array([antoine1, antoine2])

    try:
        if fix.has_unit(PASCAL):
            
            p = fix / PASCAL
            vvar, vliq, vvap = linspace_bubble_t(eos, p, antoine, N=100)
        
        elif fix.has_unit(KELVIN):

            t = fix / KELVIN
            vvar, vliq, vvap = linspace_bubble_p(eos, t, antoine, N=100)
        
        return vvar, vliq, vvap
    
    except Exception as e:
        print(e)
        return e

#%%

alias = {
    "k": KELVIN,
    "pa": PASCAL,
    "kpa": KILO * PASCAL,
    "bar": BAR,
}

def parse(s:str):
    
    rvalue, runit = s.split(" ")
    value = float(rvalue)
    unit = runit.lower()
    si_unit = alias[unit]

    return value * si_unit

assert(parse("298.15 k") == 298.15 * KELVIN)
assert(parse("298.15 K") == 298.15 * KELVIN)
assert(parse("101.32 kPa") == 101.32* KILO * PASCAL)
assert(parse("15 bar") == 15 * BAR)
assert(parse("1 " + "bar") == BAR)


#%%
# name1@name2@fixed var + unit@target's unit
# acetic acid@n-octane@343.15 K@kPa

def deserialize(sheetname):

    names, rfixed, rtunit = sheetname.split("@")
    names = names.split("#")
    rfixed = rfixed.split("#")
    
    rfixed = " ".join([".".join(rfixed[0].split("_")), rfixed[1]])
    name1 = " ".join(names[0].split("_"))
    name2 = " ".join(names[1].split("_"))

    fix = parse(rfixed)
    var = parse("1 " + rtunit)

    case = {"id1": name1, "id2": name2, "fix": fix, "var": var}

    return case


assert(
    {"id1": "acetic acid", "id2": "n-octane", "fix": 343.15 * KELVIN, "var": KILO * PASCAL }
    ==
    deserialize("acetic_acid#n-octane@343_15#K@kPa")
)
assert(
    {"id1": "acetic acid", "id2": "n-octane", "fix": 101.32 * KILO * PASCAL , "var": KELVIN }
    ==
    deserialize("acetic_acid#n-octane@101_32#kPa@K")
)

#%%

files = listdir("cases")
cases = []

for file in files:
    file = file.split(".")[0]
    case = deserialize(file)
    cases.append(case) 

# optpaths = 
        # "ppath": "../../../parameters/cpa/kontogeorgis2006.json",
        # "bpath": "../../../parameters/cpa/kontogeorgis2006_binary.json"}
#%%

vdt = []
var_list = []
vap_frac = []
liq_frac = []

for i, case in enumerate(cases):
    
    ret, dt = compute(case)
    var, liq, vap = ret
    
    vap_frac.append([state.composition[0] for state in vap])
    liq_frac.append([state.composition[0] for state in liq])
    var_list.append(var)
    vdt.append(dt)


#%%

if PLOT:

    for i, case in enumerate(cases):
        
        plt.figure()

        if case["var"].has_unit(PASCAL): 
            s = "p"
            cf = PASCAL / case["var"]

        elif case["var"].has_unit(KELVIN): 
            s = "t"
            cf = 1.0
        
        print(case["id1"],case["id2"])
        exp = pd.read_csv("./cases/" + files[i])
        x, targetx, y, targety =  exp["x"], exp[s + "x"], exp["y"], exp[s + "y"]
        plt.scatter(x, targetx, marker='o', facecolors = 'none', edgecolors = 'black')
        plt.scatter(y, targety, marker='o', facecolors = 'none', edgecolors = 'black')
    
        plt.plot(vap_frac[i], var_list[i] * cf, color = 'black')
        plt.plot(liq_frac[i], var_list[i] * cf, color = 'black')

        plt.xlabel(r"$x_1,y_1$")
        plt.ylabel(fr"${s.upper()}/{"".join(files[i].split(".")[0]).split("@")[-1]}$")
        plt.grid(True)
    
#     if save:
#         name1 = "_".join(name1.split())
#         name2 = "_".join(name2.split())
#         prefx = name1 + "@" + name2 + "@" + "".join(str(var).split())
#         os.makedirs(pltdir, exist_ok = True)
#         filename = f"{prefx}.png"
#         filepath = os.path.join(pltdir, filename)
#         plt.savefig(filepath, bbox_inches='tight', dpi=300)

    # vdt[dt] = 
#%%
idx = sorted(range(len(vdt)), key=lambda k: vdt[k])

# d = {key: val  for key, val in zip() }
info = {}
for i in idx:
    case = cases[i]
    case_name = "/".join([case["id1"], case["id2"]])
    info[case_name] = vdt[i] * 1000

#%%

info
#%%
#%%

    # names, fracs, casename = nf(filename)

    # for j, t in enumerate(sheet_names):

    #     casenames.append(casename + '@' + '_'.join(t.split('.') ) )      

    #     t = float(t) * KELVIN
    #     case = {"names": names, "y_dry_gas": fracs, "t": t}
    #     vpressure, yw, vdeveloped, vincipient, _ = compute(case, N = 100)

    #     pressure.append(vpressure * PASCAL / BAR)
        # water_content.append(yw * 1e6)
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


