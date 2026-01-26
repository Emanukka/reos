#%%
from si_units import KELVIN, BAR, PASCAL, KILO

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
