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

#####################################################################################################

ap 			= argparse.ArgumentParser(description = 'Input arguments')
ap.add_argument('-type',required=False,default="Dens",help='the directory.',type=str)
args 		= vars(ap.parse_args())


OII_1_20   = np.load("OII_01T20e3.npy")
OII_1_5    = np.load("OII_01T5e3.npy")
OII_1_10   = np.load("OII_01T1e4.npy")
OII_2_20   = np.load("OII_02T20e3.npy")
OII_2_5    = np.load("OII_02T5e3.npy")
OII_2_10   = np.load("OII_02T1e4.npy")
SII_1_5    = np.load("SII_01T5e3.npy")
SII_1_10   = np.load("SII_01T1e4.npy")
SII_1_20   = np.load("SII_01T20e3.npy")
SII_2_5    = np.load("SII_02T5e3.npy")
SII_2_10   = np.load("SII_02T1e4.npy")
SII_2_20   = np.load("SII_02T20e3.npy")
Dens_20    = np.load("Dens_T20e3.npy")
Dens_5     = np.load("Dens_T5e3.npy")
Dens_10    = np.load("Dens_T1e4.npy")

OIII_1_100  = np.load("OIII_01N1e4.npy")
OIII_2_100  = np.load("OIII_02N1e4.npy")
OIII_1_1    = np.load("OIII_01N1e2.npy")
OIII_2_1    = np.load("OIII_02N1e2.npy")
NII_1_100   = np.load("NII_01N1e4.npy")
NII_2_100   = np.load("NII_02N1e4.npy")
NII_1_1     = np.load("NII_01N1e2.npy")
NII_2_1     = np.load("NII_02N1e2.npy")
Temp_1      = np.load("Temp_N1e2.npy")
Temp_100    = np.load("Temp_N1e4.npy")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

if args['type'] == "Dens":
    fig, ax = plt.subplots(1,2,dpi=200)
    ax[0].plot(np.log10(Dens_20),OII_1_20/OII_2_20,'blue',label="T = 20,000 [K]")
    ax[0].plot(np.log10(Dens_10),OII_1_10/OII_2_10,'red',label="T = 10,000 [K]")
    ax[0].plot(np.log10(Dens_5),OII_1_5/OII_2_5,'green',label="T = 5,000 [K]")
    ax[0].axvline(x=1.5,ls=':',color="k",linewidth=1)
    ax[0].axvline(x=4.5,ls=':',color="k",linewidth=1)
    ax[0].annotate("Density sensitive",xy=(1.6,0.4))
    ax[0].set_xlabel(r"$\log_{10}(n_e) \, [\text{cm}^{-3}]$",fontsize=16)
    ax[0].set_ylabel(r"$O[II](3729\, [\AA]) / O[II](3726 \, [\AA])$",fontsize=12)
    ax[0].legend()
    ax[1].plot(np.log10(Dens_20),SII_1_20/SII_2_20,'blue',label="T = 20,000 [K]")
    ax[1].plot(np.log10(Dens_10),SII_1_10/SII_2_10,'red',label="T = 10,000 [K]")
    ax[1].plot(np.log10(Dens_5),SII_1_5/SII_2_5,'green',label="T = 5,000 [K]")
    ax[1].annotate("Density sensitive",xy=(2.1,0.75))
    ax[1].axvline(x=2,ls=':',color="k",linewidth=1)
    ax[1].axvline(x=4.5,ls=':',color="k",linewidth=1)
    ax[1].set_xlabel(r"$\log_{10}(n_e) \, [\text{cm}^{-3}]$",fontsize=16)
    ax[1].set_ylabel(r"$S[II](6731 \, [\AA]) / S[II](6717 \, [\AA])$",fontsize=12)
    ax[1].legend()
    plt.show()
else:
    fig, ax = plt.subplots(1,2,dpi=200)
    ax[0].plot(np.log10(Temp_1),OIII_1_1/OIII_2_1,'blue',label=r"n$_e$ = 100 [\text{cm}$^{-3}$]")
    ax[0].plot(np.log10(Temp_100),OIII_1_100/OIII_2_100,'red',label=r"n$_e$ = 10,000 [\text{cm}$^{-3}$]")
    ax[0].set_xlabel(r"$\log_{10}(T) \, [\text{K}]$",fontsize=16)
    ax[0].set_ylabel(r"$O[III](4363 \, [\AA]) / O[III](5007 \, [\AA])$",fontsize=12)
    ax[0].legend()
    ax[1].plot(np.log10(Temp_1),NII_1_1/NII_2_1,'blue',label=r"n$_e$ = 100 [\text{cm}$^{-3}$]")
    ax[1].plot(np.log10(Temp_100),NII_1_100/NII_2_100,'red',label=r"n$_e$ = 10,000 [\text{cm}$^{-3}$]")
    ax[1].set_xlabel(r"$\log_{10}(T) \, [\text{K}]$",fontsize=16)
    ax[1].set_ylabel(r"$N[II](5754 \, [\AA]) / N[II](6583\, [\AA])$",fontsize=12)
    ax[1].legend()
    plt.show()
