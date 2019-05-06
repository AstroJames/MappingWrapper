# Dependencies
#####################################################################################################

import numpy as np
import os
import re
from matplotlib import pyplot as plt
import argparse

# Command Line Arguments
#####################################################################################################

ap 			= argparse.ArgumentParser(description = 'Input arguments')
ap.add_argument('-spectra',required=False,default=True,help='the number of models that will be run.',type=bool)
args 		= vars(ap.parse_args())


os.system("ls slab00* > Data.txt")

f = open("Data.txt")


OIIflux     = []
OIIlambda   = []
NIIflux     = []
NIIlambda   = []
Temp        = []
Dens        = []

# For each of the data
for Data in f.readlines():

    # Open the data
    g           = open(Data.rstrip('\n'))
    lineCounter = 0
    fluxCounter = 0
    Lambda      = []
    Energy      = []
    Flux        = []
    Species     = []
    Kind        = []
    Accuracy    = []

    # for each of the lines
    for line in g.readlines():

        if lineCounter < 20:
            TempMatch = re.search(r"TEMP\.\s+:\s+(\d\D\d+\w\D\d+)",line)
            if TempMatch is not None:
                Temp.append(float(TempMatch.groups()[0]))

            DenMatch = re.search(r"DENSITY\s+:\s+(\d\D\d+\w\D\d+)",line)
            if DenMatch is not None:
                Dens.append(float(DenMatch.groups()[0]))

        if args['spectra']:
            if fluxCounter == 1:
                fluxCounter +=1
            elif fluxCounter == 2:
                lineSplit = line.split(',')
                try:
                    Lambda.append(float(lineSplit[0]))
                    Energy.append(float(lineSplit[1]))
                    Flux.append(float(lineSplit[2]))
                    Species.append(lineSplit[3].strip())
                    Kind.append(lineSplit[4].strip())
                    Accuracy.append(float(lineSplit[5]))
                except:
                    break

            FluxMatch = re.search(r"Lambda\(A\)",line)
            if FluxMatch is not None:
                fluxCounter += 1
            lineCounter += 1


    g.close()

    OIIindex = np.where(np.array(Species) == 'O  II')[0][0]
    NIIindex = np.where(np.array(Species) == 'N  II')[0][0]
    OIIlambda.append(Lambda[OIIindex])
    NIIlambda.append(Lambda[NIIindex])
    OIIflux.append(Flux[OIIindex])
    NIIflux.append(Flux[NIIindex])

f.close()
