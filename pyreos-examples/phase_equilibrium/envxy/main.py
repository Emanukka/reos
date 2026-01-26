#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from reos.eos import EquationOfState
from reos.cpa import CPAParameters
from si_units import KELVIN, PASCAL
from _tools import *

xsize = 3.15
ysize = 3.15

PLOT = False
SAVE = False
BENCH = True
PLTDIR = "./plots/"
DATAPATH = "./data/"

antoine_parameters = pd.read_csv("./antoine.csv")

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
            ret, dt = vle_functions.linspace_bubble_t(eos, p, antoine, N=100)
            vvar, vliq, vvap = ret
        
        elif fix.has_unit(KELVIN):

            t = fix / KELVIN
            ret,dt = vle_functions.linspace_bubble_p(eos, t, antoine, N=100)
            vvar, vliq, vvap = ret
        
        return vvar, vliq, vvap, dt
    
    except Exception as e:
        print(e)
        return e

#%%get files

files = os.listdir("cases")
cases = []

for file in files:
    file = file.split(".")[0]
    case = deserializer.deserialize(file)
    cases.append(case) 

#%%compute

vdt = []
var_list = []
vap_frac = []
liq_frac = []

for i, case in enumerate(cases):
    
    var, liq, vap, dt = compute(case)
    
    vap_frac.append([state.composition[0] for state in vap])
    liq_frac.append([state.composition[0] for state in liq])
    var_list.append(var)
    vdt.append(dt)


#%%plot

figures=[]

if PLOT:

    for i, case in enumerate(cases):

        fig = plt.figure()

        content = files[i].split(".")[0]

        if case["var"].has_unit(PASCAL): 
            s = "p"
            cf = PASCAL / case["var"]

        elif case["var"].has_unit(KELVIN): 
            s = "t"
            cf = 1.0
        
        exp = pd.read_csv("./cases/" + files[i])
        x, targetx, y, targety =  exp["x"], exp[s + "x"], exp["y"], exp[s + "y"]
        plt.scatter(x, targetx, marker='o', facecolors = 'none', edgecolors = 'black')
        plt.scatter(y, targety, marker='o', facecolors = 'none', edgecolors = 'black')
        plt.plot(vap_frac[i], var_list[i] * cf, color = 'black')
        plt.plot(liq_frac[i], var_list[i] * cf, color = 'black')
        plt.xlabel(r"$x_1,y_1$")
        plt.ylabel(fr"${s.upper()}/{"".join(content).split("@")[-1]}$")
        ax = fig.axes[0]
        ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
        figures.append(fig)
        plt.grid(True)

#%%save

if SAVE:
    for i,fig in enumerate(figures):

        content = files[i].split(".")[0]
        pltname = content
        os.makedirs("plots", exist_ok = True)
        filename = f"{pltname}.png"
        filepath = os.path.join(PLTDIR, filename)
        print(filepath)
        fig.savefig(filepath, bbox_inches = 'tight', dpi=300)
    # vdt[dt] = 
#%%bench

if BENCH:

    idxsort = sorted(range(len(vdt)), key=lambda k: vdt[k])

    cases_name = []
    vel = []
    t = []

    for j in idxsort:

        case = cases[j]
        case_name = "/".join([case["id1"], case["id2"]])

        cases_name.append(case_name)
        t.append(vdt[j] * 1000)
        vel.append(100 / vdt[j])

    data = pd.DataFrame({'cases':cases_name, 't (ms)':t, 'v (points/s)':vel})
    
    try:

        from cpuinfo import get_cpu_info

        cpuinfo = get_cpu_info()
        filename = "bench/" + "_".join(cpuinfo["brand_raw"].split()) + ".csv"
        data.to_csv(filename, index=False, float_format="%.2f")

    except Exception as e:
        print(e)


