# Dependencies
#####################################################################################################

import numpy as np
import os
import re
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import argparse
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']

# Command Line Arguments
#####################################################################################################

ap 			= argparse.ArgumentParser(description = 'Input arguments')
ap.add_argument('-spectra',required=False,default=True,help='the number of models that will be run.',type=bool)
ap.add_argument('-dir',required=False,default=True,help='the directory.',type=str)
ap.add_argument('-type',required=False,default="Dens",help='the directory.',type=str)
args 		= vars(ap.parse_args())

dir = args['dir']

os.system("ls {}/slab0* > Data.txt".format(dir))

f = open("Data.txt")


OII_01flux   = []
OII_02flux   = []
SII_1flux   = []
SII_2flux   = []
OIII_01flux  = []
OIII_02flux  = []
NII_01flux   = []
NII_02flux   = []
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

            DenMatch = re.search(r"N\sELECTRONS\s+:\s+(\d\D\d+\w\D\d+)",line)
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
    if args['type'] == "Dens":
        try:
            OII_01 = np.where(np.array(Lambda) == 3728.815)[0][0]
            OII_02 = np.where(np.array(Lambda) == 3726.032)[0][0]
            SII_01 = np.where(np.array(Lambda) == 6730.816)[0][0]
            SII_02 = np.where(np.array(Lambda) == 6716.440)[0][0]
            OII_01flux.append(Flux[OII_01])
            OII_02flux.append(Flux[OII_02])
            SII_1flux.append(Flux[SII_01])
            SII_2flux.append(Flux[SII_02])
        except:
            continue
    else:
        try:
            OIII_01 = np.where(np.array(Lambda) == 4363.209)[0][0]
            OIII_02 = np.where(np.array(Lambda) == 5006.843)[0][0]
            NII_01  = np.where(np.array(Lambda) == 5754.595)[0][0]
            NII_02  = np.where(np.array(Lambda) == 6583.454)[0][0]
            OIII_01flux.append(Flux[OIII_01])
            OIII_02flux.append(Flux[OIII_02])
            NII_01flux.append(Flux[NII_01])
            NII_02flux.append(Flux[NII_02])
        except:
            continue

f.close()

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

if args['type'] == "Dens":
    fig, ax = plt.subplots(2,1,dpi=200,sharex=True)
    ax[0].plot(np.log10(Dens[0:len(OII_01flux)]),np.array(OII_01flux)/np.array(OII_02flux),'blue')
    ax[0].set_ylabel(r"$O[II](3729) / O[II](3726)$",fontsize=12)
    #ax[0].annotate("[OII]3726",xy=(0.8,0.8),xycoords ='axes fraction',fontsize=16)
    ax[1].plot(np.log10(Dens[0:len(OII_01flux)]),np.array(SII_1flux)/np.array(SII_2flux),'blue')
    #ax[1].annotate("[NII]6548",xy=(0.8,0.8),xycoords ='axes fraction',fontsize=16)
    ax[1].set_xlabel(r"$n_e$",fontsize=16)
    ax[1].set_ylabel(r"$S[II](6731) / S[II](6717)$",fontsize=12)

    np.save("OII_01{}".format(dir),np.array(OII_01flux))
    np.save("OII_02{}".format(dir),np.array(OII_02flux))
    np.save("SII_01{}".format(dir),np.array(SII_1flux))
    np.save("SII_02{}".format(dir),np.array(SII_2flux))
    np.save("Dens_{}".format(dir),np.array(Dens[0:len(OII_01flux)]))
    plt.show()
    plt.close()

else:
    fig, ax = plt.subplots(2,1,dpi=200,sharex=True)
    ax[0].plot(np.log10(Temp[0:len(OIII_01flux)]),np.array(OIII_01flux)/np.array(OIII_02flux),'blue')
    ax[0].set_ylabel(r"$O[III](4363) / O[III](5007)$",fontsize=12)
    ax[1].plot(np.log10(Temp[0:len(OIII_01flux)]),np.array(NII_01flux)/np.array(NII_02flux),'blue')
    ax[1].set_xlabel(r"$T_e \, [\text{K}]$",fontsize=16)
    ax[1].set_ylabel(r"$N[II](5754) / N[II](6583)$",fontsize=12)

    np.save("OIII_01{}".format(dir),np.array(OIII_01flux))
    np.save("OIII_02{}".format(dir),np.array(OIII_02flux))
    np.save("NII_01{}".format(dir),np.array(NII_01flux))
    np.save("NII_02{}".format(dir),np.array(NII_02flux))
    np.save("Temp_{}".format(dir),np.array(Temp[0:len(OIII_01flux)]))
    plt.show()
    plt.close()
