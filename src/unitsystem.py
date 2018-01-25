#Units and physical constants

import scipy.constants

pc=scipy.constants.physical_constants

#Boltzmann constant
kB=pc['Boltzmann constant in eV/K'][0]

#Avogadro's number
avogadro=pc['Avogadro constant'][0]

#Unit conversion factors
#These are designed to go from the designated unit to model units
#This way, inputs are specfied as X*unit to convert them from "unit" to the corresponding model uit
#and outputs are converted to other units as Y/unit.

#Concentration
molar=avogadro/1e24
millimolar=molar*1e-3
micromolar=molar*1e-6

#Time
second=1e-9
millisecond=1e-6
microsecond=1e-3
nanosecond=1

