# Kinetics

Kinetics code allows building Matlab inputs to calculate:
- thermodynamic properties of the different speces: Cp, H, S and G
- reaction constants
- constant temperature kinetics
- variable temperature kinetics
- temperature programmed desorptions
- further analyses will be implemented soon

according to the process described in the mechanisms included and the reaction conditions specified.

# Usage

Upon having a plain input with few initial data (conditions, elementary reactions and species), kinetics.sh calls input2mk.py. The python script is built on [ASE](https://wiki.fysik.dtu.dk/ase/) and reads outputs files from VASP and FHAims to extract the key information for each system to write the actual kinetics2Matlab input. Kinetics.sh subsequentially calls Kinetics2Matlab.pl to generate different Matlab inputs.

# License
The software here described is licenced according to the LICENCE.md and it can be used upon recognition.
Please cite [https://doi.org/10.1039/D1NA00015B](https://doi.org/10.1039/D1NA00015B).

# Case Study
A detailed methodology and flow work is found on:
- [https://doi.org/10.1039/D1NA00015B](https://doi.org/10.1039/D1NA00015B)

