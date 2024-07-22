# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 13:17:02 2020

@author: ndron
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


def read_data(readname, writename, element, targ, A):
    print(readname, writename, element, targ, A)
    with open(readname, 'r') as data_file:
        saved_lines = []
        #runs through each line in the file
        for linenum, data_line in enumerate(data_file):
            if (targ == "be"):
                lower = 22
                upper = 67
            if (targ == "epoxy"):
                lower = 24
                upper = 96
            elif (targ == "fiber"):
                lower = 23
                upper = 68
            elif (targ == "CD2"):
                #100->100,000 keV
                lower = 23
                upper = 103
            elif (targ == "CH2"):
                #100->100,000 keV
                lower = 23
                upper = 103
            elif (targ == "C"):
                #100->100,000 keV
                lower = 22
                upper = 102
            elif (targ == "Au"):
                #100->100,000 keV
                lower = 22
                upper = 102
            elif (targ == "Si"):
                #100->100,000 keV
                lower = 22
                upper = 102
            else:
                sys.exit("need a type of targ")
            if (linenum > 22 and linenum < 111): #for H-C files
            #if (linenum > lower and linenum < upper): #for Be files
                #parse each line, save it in a list
                splitline = data_line.split()
                print(splitline)
                if (splitline[1] == 'keV'):
                    dEdx = str(round(float(splitline[2]), 4))
                    E = str(round(float(splitline[0])/1000/A, 4))
                    line_to_write = E + " " + dEdx + '\n'
                elif (splitline[1] == 'MeV'):
                    dEdx = str(round(float(splitline[2]), 4))
                    E = str(round(float(splitline[0])/A, 4))
                    line_to_write = E + " " + dEdx + '\n'
                elif (splitline[1] == 'GeV'):
                    dEdx = str(round(float(splitline[2]), 4))
                    E = str(round(float(splitline[0])*1000/A, 4))
                    line_to_write = E + " " + dEdx + '\n'
                else:
                    print("didn't understand " + splitline[1])
                saved_lines.append(line_to_write)


    with open(writename, 'w') as writer:
        if (targ == "be"):
            writer.write("energy loss of " + element + " in Beryllium MeV/A--Mev/mg/cm2" + "\n")
        if (targ == "epoxy"):
            writer.write("energy loss of " + element + " in Expoxy (PCB) MeV/A--Mev/mg/cm2" + "\n")
        elif (targ == "fiber"):
            writer.write("energy loss of " + element + " in Fiber Ary (C9H10, D=1.032g/cm3) MeV/A--Mev/mg/cm2" + "\n")
        elif (targ == "CD2"):
            writer.write("energy loss of " + element + " in Deuterated polyethylene (C2D4, D=1.06g/cm3) MeV/A--Mev/mg/cm2" + "\n")
        elif (targ == "CH2"):
            writer.write("energy loss of " + element + " in polypropylene (C3H6, D=0.9g/cm3) MeV/A--Mev/mg/cm2" + "\n")
        elif (targ == "C"):
            writer.write("energy loss of " + element + " in Carbon (C, D=1g/cm3) MeV/A--Mev/mg/cm2" + "\n")
        elif (targ == "Au"):
            writer.write("energy loss of " + element + " in Gold (Au) MeV/A--Mev/mg/cm2" + "\n")
        elif (targ == "Si"):
            writer.write("energy loss of " + element + " in Silicon (Si) MeV/A--Mev/mg/cm2" + "\n")

        writer.write(str(len(saved_lines)) + "\n")

        for line in saved_lines:
            writer.write(line)

#srim runs for CD2 target
#run srim from 100->100,00keV, density=16/14*0.93=1.06, increase the mass of hydrogen to 2.014
# read_data("Hydrogen in CD2.txt", "Hydrogen_CD2.loss", "Hydrogen", "CD2", 1.008)
# read_data("Helium in CD2.txt", "Helium_CD2.loss", "Helium", "CD2", 4.003)
# read_data("Lithium in CD2.txt", "Lithium_CD2.loss", "Lithium", "CD2", 7.016)
# read_data("Beryllium in CD2.txt", "Beryllium_CD2.loss", "Beryllium", "CD2", 9.012)
# read_data("Boron in CD2.txt", "Boron_CD2.loss", "Boron", "CD2", 11.009)
# read_data("Carbon in CD2.txt", "Carbon_CD2.loss", "Carbon", "CD2", 12.000)

#srim runs for polypropylene target
#run srim from 100->100,000keV, density=0.9
# read_data("Hydrogen in Polypropylene.txt", "Hydrogen_CH2.loss", "Hydrogen", "CH2", 1.008)
# read_data("Helium in Polypropylene.txt", "Helium_CH2.loss", "Helium", "CH2", 4.003)
# read_data("Lithium in Polypropylene.txt", "Lithium_CH2.loss", "Lithium", "CH2", 7.016)
# read_data("Beryllium in Polypropylene.txt", "Beryllium_CH2.loss", "Beryllium", "CH2", 9.012)
# read_data("Boron in Polypropylene.txt", "Boron_CH2.loss", "Boron", "CH2", 11.009)
# read_data("Carbon in Polypropylene.txt", "Carbon_CH2.loss", "Carbon", "CH2", 12.000)

#srim runs for Carbon target
#run srim from 100->100,000keV, density=2.253
# read_data("Hydrogen in Carbon.txt", "Hydrogen_C.loss", "Hydrogen", "C", 1.008)
read_data("Helium in Carbon.txt", "Helium_C.loss", "Helium", "C", 4.003)
# read_data("Lithium in Carbon.txt", "Lithium_C.loss", "Lithium", "C", 7.016)
# read_data("Beryllium in Carbon.txt", "Beryllium_C.loss", "Beryllium", "C", 9.012)
# read_data("Boron in Carbon.txt", "Boron_C.loss", "Boron", "C", 11.009)
# read_data("Carbon in Carbon.txt", "Carbon_C.loss", "Carbon", "C", 12.000)

#srim runs for Gold target
#run srim from 100->100,000keV, density=19.311
# read_data("Hydrogen in Gold.txt", "Hydrogen_Au.loss", "Hydrogen", "Au", 1.008)
# read_data("Helium in Gold.txt", "Helium_Au.loss", "Helium", "Au", 4.003)
# read_data("Lithium in Gold.txt", "Lithium_Au.loss", "Lithium", "Au", 7.016)
# read_data("Beryllium in Gold.txt", "Beryllium_Au.loss", "Beryllium", "Au", 9.012)
# read_data("Boron in Gold.txt", "Boron_Au.loss", "Boron", "Au", 11.009)
# read_data("Carbon in Gold.txt", "Carbon_Au.loss", "Carbon", "Au", 12.000)

#srim runs for Silicon detectors
#run srim from 100->100,000keV, density=2.32120
#read_data("Hydrogen in Silicon.txt", "Hydrogen_Si.loss", "Hydrogen", "Si", 1.008)
#read_data("Helium in Silicon.txt", "Helium_Si.loss", "Helium", "Si", 4.003)
#read_data("Lithium in Silicon.txt", "Lithium_Si.loss", "Lithium", "Si", 7.016)
#read_data("Helium in Epoxy.txt", "Helium_epoxy.loss", "Helium", "epoxy", 4.003)



#Srim runs needed for Ca36 experiment, loss in Be and scintillating fiber array
# read_data("Phosphorus in Beryllium.txt", "Phosphorus_Be.loss", "Phosphorus30", "be",30)
# read_data("Sulfur in Beryllium.txt", "Sulfur_Be.loss", "Sulfur31", "be", 31)
# read_data("Chlorine in Beryllium.txt", "Chlorine_Be.loss", "Chlorine32", "be", 32)
# read_data("Argon in Beryllium.txt", "Argon_Be.loss", "Argon34", "be", 34)
# read_data("Potassium in Beryllium.txt", "Potassium35_Be.loss", "Potassium35", "be", 35)
# read_data("Calcium in Beryllium.txt", "Calcium_Be.loss", "Calcium36", "be", 36)
# read_data("Scandium in Beryllium.txt", "Scandium_Be.loss", "Scandium37", "be", 37)
# read_data("Titanium in Beryllium.txt", "Titanium_Be.loss", "Titanium39", "be", 39)
# read_data("Phosphorus in  H- C.txt", "Phosphorus_fiber.loss", "Phosphorus30", "fiber",30)
# read_data("Sulfur in  H- C.txt", "Sulfur_fiber.loss", "Sulfur31", "fiber",31)
# read_data("Chlorine in  H- C.txt", "Chlorine_fiber.loss", "Chlorine32", "fiber",32)
# read_data("Argon in  H- C.txt", "Argon_fiber.loss", "Argon", "fiber",34)
# read_data("Potassium in  H- C.txt", "Potassium35_fiber.loss", "Potassium", "fiber", 35)
# read_data("Calcium in  H- C.txt", "Calcium_fiber.loss", "Calcium", "fiber",36)
# read_data("Scandium in  H- C.txt", "Scandium_fiber.loss", "Scandium", "fiber",37)
# read_data("Titanium in  H- C.txt", "Titanium_fiber.loss", "Titanium", "fiber",39)

