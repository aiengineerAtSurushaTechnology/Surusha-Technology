import rdkit
from rdkit import Chem
from chempy import Substance

Chem.MolFromSmiles('CCO')

from rdkit.Chem import AllChem
AllChem.ReactionFromSmarts('C(=O)O>>C=O')
ferricyanide = Substance.from_formula('Fe(CN)6-3')
ferricyanide.composition == {0: -3, 26: 1, 6: 6, 7: 6}
print(ferricyanide.unicode_name)
print(ferricyanide.latex_name + ", " + ferricyanide.html_name)
print('%.3f' % ferricyanide.mass)

# Balancing stoichiometry of a chemical reaction
from chempy import balance_stoichiometry 
reac, prod = balance_stoichiometry({'NH4ClO4','Al'}, {'Al2O3','HCl','H2O','N2'})
from pprint import pprint
pprint(dict(reac))
pprint(dict(prod))

from chempy import mass_fractions
for fractions in map(mass_fractions, [reac,prod]):
    pprint({k: '{0:.3g} wt%'.format(v*100) for k, v in fractions.items()})

from chempy import Substance
substances = {s.name: s for s in [
    Substance('pancake', composition=dict(eggs=1, spoons_of_flour=2, cups_of_milk=1)),
    Substance('eggs_6pack', composition=dict(eggs=6)),
    Substance('milk_carton', composition=dict(cups_of_milk=4)),
    Substance('flour_bag', composition=dict(spoons_of_flour=60))
]}
pprint([dict(_) for _ in balance_stoichiometry({'eggs_6pack','milk_carton', 'flour_bag'},
                                              {'pancake'}, substances = substances)])


# Chempy can also handle reactions with linear dependencies e.g.,
pprint([dict(_) for _ in balance_stoichiometry({'C','O2'}, {'CO2', 'CO'})])

pprint([dict(_) for _ in balance_stoichiometry({'C','O2'}, {'CO2', 'CO'}, underdetermined=None)])

# Balancing Reactions
from chempy import Equilibrium
from sympy import symbols 
K1, K2, Kw = symbols('K1 K2 Kw')
e1 = Equilibrium({'MnO4-': 1, 'H+': 8, 'e-': 5}, {'Mn+2': 1, 'H2O': 4}, K1)
e2 = Equilibrium({'O2': 1, 'H2O': 2, 'e-': 4}, {'OH-': 4}, K2)
coeff = Equilibrium.eliminate([e1, e2], 'e-')
coeff

redox = e1*coeff[0] + e2*coeff[1]
print(redox)


autoprot = Equilibrium({'H2O':1}, {'H+': 1, 'OH-': 1}, Kw)
n = redox.cancel(autoprot)
n


redox2 = redox + n*autoprot
print(redox2)

r = Reaction.from_string("H2O -> H+ + OH-; 1e-4/s")
from chempy.units import to_unitless, default_units as u
to_unitless(r.param, 1/u.minute)

# Chemical equilibria
# if we want to predict pH of a bicarbonate solution we simply just need pKa and pKw values:
from collections import defaultdict
from chempy.equilibria import EqSystem
eqsys = EqSystem.from_string("""HCO3- = H+ + CO3-2; 10**-10.3
H2CO3 = H+ + HCO3-; 10**-6.3
H2O = H+ + OH-; 10**-14/55.4
""")
arr, info, sane = eqsys.root(defaultdict(float, {'H2O':55.4, 'HCO3-':1e-2}))
conc = dict(zip(eqsys.substances, arr))
from math import log10
print("pH: %.2f" % -log10(conc['H+']))

# here is another for ammonia:
from chempy import Equilibrium
from chempy.chemistry import Species
water_autop = Equilibrium({'H2O'}, {'H+', 'OH-'}, 10**-14)
ammonia_prot = Equilibrium({"NH4+"}, {"NH3", "H+"}, 10**-9.24)
substances = [Species.from_formula(f) for f in 'H2O OH- H+ NH3 NH4+'.split()]
eqsys = EqSystem([water_autop, ammonia_prot], substances)
print('\n'.join(map(str, eqsys.rxns)))

init_conc = defaultdict(float, {'H2O': 1, 'NH3': 0.1})
x, sol, sane = eqsys.root(init_conc)
assert sol['success'] and sane
print(', '.join('%.2g' % v for v in x))


from chempy.electrolytes import ionic_strength
ionic_strength({'Fe+3': 0.050, 'ClO4-': 0.150}) == .3

from chempy.henry import Henry
kH_O2 = Henry(1.2e-3, 1800, ref='carpenter_1966')
print('%.1e' % kH_O2(298.15))

from chempy import ReactionSystem  # The rate constants below are arbitrary
rsys = ReactionSystem.from_string("""2 Fe+2 + H2O2 -> 2 Fe+3 + 2 OH-; 42
    2 Fe+3 + H2O2 -> 2 Fe+2 + O2 + 2 H+; 17
    H+ + OH- -> H2O; 1e10
    H2O -> H+ + OH-; 1e-4""")  # "[H2O]" = 1.0 (actually 55.4 at RT)
from chempy.kinetics.ode import get_odesys
odesys, extra = get_odesys(rsys)
from collections import defaultdict
import numpy as np
tout = sorted(np.concatenate((np.linspace(0, 23), np.logspace(-8, 1))))
c0 = defaultdict(float, {'Fe+2': 0.05, 'H2O2': 0.1, 'H2O': 1.0, 'H+': 1e-2, 'OH-': 1e-12})
result = odesys.integrate(tout, c0, atol=1e-12, rtol=1e-14)
import matplotlib.pyplot as plt
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
for ax in axes:
    _ = result.plot(names=[k for k in rsys.substances if k != 'H2O'], ax=ax)
    _ = ax.legend(loc='best', prop={'size': 9})
    _ = ax.set_xlabel('Time')
    _ = ax.set_ylabel('Concentration')
_ = axes[1].set_ylim([1e-13, 1e-1])
_ = axes[1].set_xscale('log')
_ = axes[1].set_yscale('log')
_ = fig.tight_layout()

from chempy import Substance
from chempy.properties.water_density_tanaka_2001 import water_density as rho
from chempy.units import to_unitless, default_units as u
water = Substance.from_formula('H2O')
for T_C in (15, 25, 35):
    concentration_H2O = rho(T=(273.15 + T_C)*u.kelvin, units=u)/water.molar_mass(units=u)
    print('[H2O] = %.2f M (at %d °C)' % (to_unitless(concentration_H2O, u.molar), T_C))
