'''

	by Alberto Roldan:	04/2020

	Reads information from outputs to generate a input file for microkinetics

'''

import os
import re
import subprocess
import sys
import numpy as np
# from numbers import Number
from ase import Atoms
from ase.io import read, write

Parameters = ["SYSTEM", "SYSPATH","FREQPATH", "ISITES", "SSITES", "IPRESSURE", "RPRESSURE", "ICOVERAGE", "RCOVERAGE","END"]

#print("\n")


def Common_Properties(system, name):

# Exports the geometry in XYZ format
	if not os.path.exists("./XYZ"):
		os.makedirs("./XYZ")
	write(str("./XYZ/" + name + ".xyz"), Atoms(system))

# Get the ground state multiplicity
	try:
		mag_moment = system.get_magnetic_moment() + 1
	except:
		mag_moment = 0.00 + 1

	return str("./XYZ/" + name + ".xyz"), mag_moment

#..............................................................................................................

def System_Properties(File, system, ssites):
	del system.constraints
	symmetry_factor = None
	Area = None

	system_type = "Molecule"
# Decide if the system is a naked SURFACE through SSites
	for i in range(3):
		if min(Atoms(system).get_scaled_positions()[:, i]) < 0.1 and \
				max(Atoms(system).get_scaled_positions()[:, i]) > 0.9:
			if ssites is not None:
				system_type = "Surface"
# Decide if the system combines adsorbate + adsorbent (CATALYSTS)
			else:
				system_type = "Catalyst"

	if system_type == "Molecule":
		pos = [1 for i in system.get_scaled_positions() if i[0] == 0. and i[1] == 0. and i[2] == 0.]
		if len(pos) > 1:
			if ssites is not None:
				system_type = "Surface"
# Decide if the system combines adsorbate + adsorbent (CATALYSTS)
			else:
				system_type = "Catalyst"
#		answer = input("Is " + name +" an isolated molecule? ")
#		if answer.startswith('y') is True or answer.startswith('Y') is True:
#			system_type = "Molecule"

# Searches for the symmetry of the Molecule
	if system_type == "Molecule":

		output = open(File,"r")
# if File is from VASP (OUTCAR) ----------------------------------------------------------------- not completed
		try:
			out = [out.strip() for out in output if re.search("The static configuration has the point symmetry",out) is not None][-1]
			symmetry_factor = float(out[-3])    ##### CHECK!
		except:
			symmetry_factor = input("What's the symmetry factor for " + name +"? ")
# if File is from FHI-AIMS
		try:
			out = [out.strip() for out in output if re.search("Space group",out) is not None][-1]
			symmetry_factor = float(out[-3])    ##### CHECK!
		except:
			symmetry_factor = input("What's the symmetry factor for " + name +"? ")


# Decide if the system is a naked SURFACE through SSites

	if system_type == "Surface":
# find the area of the surface in m^2
		if system.get_pbc().all == True:
			Area = system.cell[0][0] * system.cell[1][1] * 1e-20
		else:
			lattice = [system.cell[i][i] for i in range(len(system.cell)) if system.get_pbc()[i] == True]
			Area = lattice[0] * lattice[1] * 1e-20


	return system_type, symmetry_factor, Area

#..............................................................................................................

def Frequencies(inputfile, system_type, chemical_symbols, software):
	ff = open(inputfile)
	natoms = len(chemical_symbols)
# Enters in the loop to read the frequencies file
	frequencies = []
	frequencies_2D = []

	if software == "VASP":
		lines = ff.readlines()
		nline = 0
		frequencies_2D = []
		while nline <= len(lines)-1:
			words = [i.strip() for i in lines[nline].split() if i]
			nline += 1
# Looking for frequencies from VASP
			if len(words) > 7 and "f" in words[1] and 'THz' in words[4] and 'cm-1' in words[8]:
				frequencies.append(float(words[7]))
			elif len(words) > 6 and "f/i" in words[1] and 'THz' in words[3] and 'cm-1' in words[7]:
				frequencies.append(-float(words[6]))
# Looking for frequency displacements and freq2D from VASP
				nline += 1
				positions = []
				new_positions = []
				if system_type == "Molecule":
					for i in range(len(chemical_symbols)):
						words = [ i.strip() for i in lines[nline].split() if i]
						nline += 1
						positions.append([float(words[j]) for j in range(3)])
						new_positions.append([float(words[j]) + float(words[j+3]) for j in range(3)])

					system = Atoms(chemical_symbols, positions = positions, pbc=[True,True,True])
					new_system = Atoms(chemical_symbols, positions = new_positions, pbc=[True,True,True])
					shift = new_system.get_center_of_mass() - system.get_center_of_mass()
					shift_module = np.sqrt(shift[0]**2 + shift[1]**2 + shift[2]**2)
					if shift_module > 0.1:
						frequencies_2D += [frequencies[-1]]

	elif software == "FHI-aims":
		word = [i.strip() for i in ff.readline().split(" ") if i][-1]
		try:
			unconstrained_natoms = int(word)
		except ValueError:
			unconstrained_natoms = 0
			print("   "+system_type+" is not a valid for Frequecies!")
		if unconstrained_natoms > 0:
			lines = ff.readlines()
			nline = 1
			frequencies_2D = []
			while nline <= len(lines)-1:
				words = [i.strip() for i in lines[nline].split() if i]
				nline += 1
# Looking for frequencies from jmol (FHI-Aims)
				if len(words) > 0 and words[0] == "Mode":
					if words[4].endswith('i') is True:
						words[4] = float(words[4].rstrip("i")) * -1
					frequencies += [float(words[4])]

# Looking for frequency displacements and freq2D from jmol (FHI-Aims)
					positions = []
					new_positions = []
					unconstrained_elements = []
					words = [i.strip() for i in lines[nline].split() if i]
					nline += 1
					if system_type == "Molecule" and len(words) == 7 and words[0] in chemical_symbols:
						unconstrained_elements.append(words[0])
						positions.append([float(words[j]) for j in range(1,4)])
						new_positions.append([float(words[j]) + float(words[j+3]) for j in range(1,4)])
						for i in range(unconstrained_natoms -1):
							words = [i.strip() for i in lines[nline].split() if i]
							nline += 1
							unconstrained_elements.append(words[0])
							positions.append([float(words[j]) for j in range(1,4)])
							new_positions.append([float(words[j]) + float(words[j+3]) for j in range(1,4)])

						system = Atoms(unconstrained_elements, positions = positions, pbc=[True,True,True])
						new_system = Atoms(unconstrained_elements, positions = new_positions, pbc=[True,True,True])
						shift = new_system.get_center_of_mass() - system.get_center_of_mass()
						shift_module = np.sqrt(shift[0]**2 + shift[1]**2 + shift[2]**2)
						if shift_module > 0.001:
							frequencies_2D += [abs(frequencies[-1])]

# try to get the Infrared Spectrum
#		if len(frequencies) > 0:
#				print("infrared spectrum for",File,"not available!")
#                               system.get_dipole_moment()
#                               ir = Infrared(system)
#                               print (ir.summary())


# Frequencies
	frequencies_2D.sort(reverse = True)
	frequencies.sort(reverse = True)

# Adjust the number of freq to 3N-6 for molecules
	if system_type == "Molecule":
		if natoms < 3:          ### linear
			frequencies = [frequencies[i] for i in range(len(frequencies)) if i <= 3 * natoms - 6]  # 6 because starts by 0
			frequencies_2D = [frequencies_2D[i] for i in range(len(frequencies_2D)) if i <= 2 * natoms - 4]
		else:
			frequencies = [frequencies[i] for i in range(len(frequencies)) if i <= 3 * natoms - 7]  # 7 because starts by 0
			frequencies_2D = [frequencies_2D[i] for i in range(len(frequencies_2D)) if i <= 2 * natoms - 5]
		if len(frequencies_2D) == 0:
			frequencies_2D = [0.0]
	else:
		frequencies_2D = []

	return frequencies,frequencies_2D

#############################################################################################################


# reads the initial file
try:
	InputFile = sys.argv[1]
	f = open(InputFile)
	lines = f.readlines()
except IOError:
	raise Exception("   "+InputFile+" is not a valid file!")

OutputFile = str(InputFile+".mk.in")

# wipes the variables
File = None
FreqFile = None
system = None
name = None
frequencies = None
frequencies_2D = None
Inertia_moments = None
ssites = None
sites = None
ipressure = None
rpressure = None
icoverage = None
rcoverage = None

ifile = open(OutputFile,'a+')

# Enters in the loop to read the initial file looking for particular tags
for line in lines:
	words = line.split(" ")
	words = [i.strip() for i in words if i]

# Looking for parameters in the InputFile
	if words[0] == "SYSTEM":
		name = words[2]
	elif words[0] == "SYSPATH":
		File = words[2]
	elif words[0] == "FREQPATH":
		FreqFile = words[2]
	elif words[0] == "ISITES":
		sites = words[2:]
	elif words[0] == "SSITES":
		ssites = words[2:]
	elif words[0] == "IPRESSURE":
		ipressure = words[2]
	elif words[0] == "RPRESSURE":
		rpressure = words[2:]
	elif words[0] == "ICOVERAGE":
		icoverage = words[2]
	elif words[0] == "RCOVERAGE":
		rcoverage = words[2:]


	if words[0] == "END" or words[0] == "end":
# seeks for common and system properties
		try:
			fd = open(File)
			system = read(File, index=-1)
			xyzPath, mag = Common_Properties(system, name)
			system_type, symmetry_factor, Area = System_Properties(File, system, ssites)
			Inertia_moments = system.get_moments_of_inertia()

			software = None
			while software == None:
				line = fd.readline().split(" ")
				words = [i.strip() for i in line if i]
				for word in words:
					if word == "FHI-aims":
						software = "FHI-aims"
					elif word.startswith('vasp') is True:
						software = "VASP"

		except IOError:
			raise Exception("   "+File+" is not a valid file!")

# Works out the 3D and 2D frequencies
		if FreqFile is None and software == "VASP":
			FreqFile = File
		try:
			frequencies, frequencies_2D = Frequencies(FreqFile, system_type, system.get_chemical_symbols(), software)
		except IOError:
			raise Exception("   "+FreqFile+" is not a valid file for FREQUENCIES!")

# writes the information required from each system
#		print("   SYSTEM = %s" %(name) )
		ifile.write ("\nSYSTEM = %s\n" %(name))
		ifile.write ("SYSPATH = %s\n" %(File))
		ifile.write ("XYZPATH = %s\n" %(str("./XYZ/" + name + ".xyz")))
		if FreqFile is not None:
			ifile.write ("FREQPATH = %s\n" %(FreqFile))
		ifile.write (" E0 = %f\n" %(system.get_total_energy()))
		ifile.write (" DEGENERATION = %d\n" %(np.abs(round(mag))))
		if frequencies:
			ifile.write (" FREQ =")
			for freq in frequencies:
				ifile.write (" %.1f" %(float(freq)))
			ifile.write ("\n")
#		else:
#			ifile.write (" FREQ = 0.0\n")
		if system_type == "Molecule":
			if frequencies_2D is not None:
				if len(frequencies_2D) < 1:
					ifile.write (" FREQ2D = 0.0")
				else:
					ifile.write (" FREQ2D =")
					for freq2D in frequencies_2D:
						ifile.write (" %.1f" %(float(freq2D)))
				ifile.write ("\n")
			ifile.write (" IMASS =")
			masses = []
			n_masses = [0]*len(system.get_masses())
			n=-1
			for mass in system.get_masses():
				if mass not in masses:
					masses.append(mass)
					n += 1
					n_masses[n] = 1
				else:
					n_masses[n] += 1
			for mass in masses:
				ifile.write (" %.3f" %(float(mass)))
			ifile.write ("   # amu\n")
			ifile.write (" INATOMS =")
			for n in n_masses:
				if n > 0:
					ifile.write (" %d " %(int(n)))
			ifile.write ("\n")
			ifile.write (" SYMFACTOR = %d\n" %(int(symmetry_factor)))
			ifile.write (" INERTIA =")
			for IM in Inertia_moments:
				ifile.write (" %.3g" %(IM*1.66053904E-47))
			ifile.write ("   # kg*m^2\n")
			ifile.write (" ISITES =")
			for s in sites:
				ifile.write (" %s" %(s))
			ifile.write ("\n")
			if ipressure is not None:
				ifile.write (" IPRESSURE = %.3f\n" %(float(ipressure)))
			elif rpressure is not None:
				ifile.write (" RPRESSURE = %.3f %.3f %.3f\n" %(float(rpressure[0]),
															  float(rpressure[1]),
															  float(rpressure[2])))       # initial, final and pressure step
			else:
				ifile.write (" IPRESSURE = 0.0\n")
		if system_type == "Surface":
			ifile.write (" ISITES =")
			for s in ssites:
				ifile.write (" %s" %(s))
			ifile.write ("\n")
			ifile.write (" IACAT = %.5g    # m^2\n" %(float(Area)))
		if system_type == "Catalyst":
			ifile.write (" ISITES =")
			for s in sites:
				ifile.write (" %s" %(s))
			ifile.write ("\n")
			if icoverage is not None:
				ifile.write (" ICOVERAGE = %.3f\n" %(float(icoverage)))
			elif rcoverage is not None:
				ifile.write (" RCOVERAGE = %.3f %.3f %.3f\n" %(float(rcoverage[0]),
															  float(rcoverage[1]),
															  float(rcoverage[2])))       # initial, final and coverage step
			else:
				ifile.write (" ICOVERAGE = 0.0\n")

		ifile.write ("END\n")

# wipes the variables at the end of each system
		system = None
		File = None
		FreqFile = None
		name = None
		frequencies = None
		frequencies_2D = None
		Inertia_moments = None
		ssites = None
		sites = None
		ipressure = None
		rpressure = None
		icoverage = None
		rcoverage = None
	elif words[0] not in Parameters:
# write lines not related with the data from the systems
	   ifile.write (line)
ifile.close()
