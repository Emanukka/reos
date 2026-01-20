import sys
import reos
mods = sys.modules

assert(mods.keys().__contains__("reos"))
assert(mods.keys().__contains__("reos.consts"))
assert(mods.keys().__contains__("reos.state"))
assert(mods.keys().__contains__("reos.cubic"))
assert(mods.keys().__contains__("reos.cpa"))
assert(mods.keys().__contains__("reos.eos"))

