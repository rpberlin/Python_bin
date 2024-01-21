#!/Library/Frameworks/Python.framework/Versions/3.10/bin/python3

import cantera as ct

sodium = ct.Element('Na')
sodium.TP = 300,101325
water = ct.Water()
water.TP = 300, 101325
mix_phases = [(sodium, 1.0), (water, 0.0)]
mix = ct.Mixture(mix_phases)
mix.TP = 300,101325
mix.equilibrate('UV')


#nSpecies = len(species_to_track)
#Z_burnt = g1.mixture_fraction(fuel, air)
