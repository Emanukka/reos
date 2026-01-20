#%%
import numpy as np
from math import isclose

from reos.cubic import CubicPureRecord,CubicBinaryRecord,CubicParameters
from reos.eos import EquationOfState
from reos.consts import Consts
# from reos.cubic import CubicBinaryRecord
R = Consts.ideal_gas_const()
#%%instance

name = "water"
molar_weight = 18.01528  # g/mol

set1 = CubicPureRecord.new(name, molar_weight, a0=0.0, b=0.0, c1=0.0, tc=0.0)
set2 = CubicPureRecord.new(name, molar_weight, tc=0.0, pc=0.0,w=0.0)
set3 = CubicPureRecord.new(name, molar_weight, b=0.0, c1=0.0, tc=0.0, pc=0.0,w=0.0)

print(set1)
print(set2)
print(set3)


#%%err
try:
    err = CubicPureRecord.new(name, molar_weight)
except TypeError as e:
    
    # using panic! macro
    # thread '<unnamed>' panicked at src/py_parameters/pyrecords.rs:99:9:
    # parameters empty!

    # returning a PyResult (doesnt panic!)
    print(f"Caught expected exception: {e}")


# %%

id1 = "water"
id2 = "co2"
set1 = CubicBinaryRecord.new(id1, id2, aij=0.1, bij=0.01)
set2 = CubicBinaryRecord.new(id1, id2, kij=0.2)

print(set1)
print(set2)

#%%err
try:
    _ = CubicBinaryRecord.new(id1, id2, aij=0.1, bij=0.01)
except TypeError as e:
    
    print(f"Caught expected exception: {e}")

#%%SRK
name = "water"
molar_weight = 18.01528 

pr1 = CubicPureRecord.new(name, molar_weight, a0=0.0, b=0.0, c1=0.0, tc=0.0)
name = "co2"
molar_weight = 44.009  
pr2 = CubicPureRecord.new(name, molar_weight, a0=0.0, b=0.0, c1=0.0, tc=0.0)

pure_records = [pr1, pr2]
binary_record = [CubicBinaryRecord.new("water", "co2", aij=0.1, bij=0.01)]

p = CubicParameters.from_records(pure_records, binary_record)

print(p)


#%% srk eos
name = "water"
molar_weight = 18.01528  # g/mol

pr1 = CubicPureRecord.new(name, molar_weight, a0=0.12277, b=0.0145e-3, c1=0.6736, tc=647.14)
p = CubicParameters.from_records([pr1], opt = "srk")

#%%

srk = EquationOfState.cubic(p)
t = 298.15
d = 1000.0
x = np.array([1.0])

P = srk.pressure(t, d, x)
Pres = P - srk.ideal_gas_pressure(t, d)
Pres_reduced = Pres / R / t
chem_pot_res_reduced = srk.lnphi(t, d, x)[0] + np.log(srk.compressibility(t, d, x))
helmholtz_res_reduced = srk.helmholtz(t, d, x) / R / t
entropy_res_reduced = srk.entropy(t, d, x) / R

assert(isclose(Pres_reduced, -57.5159551979349, rel_tol=1e-10))
assert(isclose(chem_pot_res_reduced, -0.115660251059, rel_tol=1e-10))
assert(isclose(helmholtz_res_reduced, -0.058144295861, rel_tol=1e-10))
assert(isclose(entropy_res_reduced, -0.041951593945, rel_tol=1e-10))



        # assert_relative_eq!(val[0], , epsilon = 1e-10)
        # assert_relative_eq!(val, -0.041951593945, epsilon = 1e-10)
        # assert_relative_eq!(val, ,epsilon = 1e-10)

# %%
