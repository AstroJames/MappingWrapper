# Dependencies
#####################################################################################################

import numpy as np
import os
import re
from matplotlib import pyplot as plt
from matplotlib import rc
import matplotlib as mpl
import argparse
from scipy import signal
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
from scipy.optimize import least_squares
from scipy.optimize import curve_fit
from scipy.integrate import simps
from sklearn.metrics import mean_squared_error, r2_score

# Command Line Arguments
#####################################################################################################

ap 			= argparse.ArgumentParser(description = 'Input arguments')
ap.add_argument('-data',required=False,help='the directory.',type=str)
args 		= vars(ap.parse_args())

#####################################################################################################

def ReadSpectrum(data):
    dir     = "../ProvidedSpectrum/"
    data    = dir + data

    Wavelength_A        = []
    Flux_erg_s_cm2_A    = []
    lineCounter         = 0


    f = open(data)
    for line in f.readlines():
        if lineCounter == 0:
            lineCounter +=1
            continue

        line = line.rstrip('\n')
        try:
            Flux        = float(line.split(' ')[10])
            Flux_erg_s_cm2_A.append(Flux)
            Wavelength  = float(line.split(' ')[7])
            Wavelength_A.append(Wavelength)
        except:
            try:
                Flux        = float(line.split(' ')[9])
                Flux_erg_s_cm2_A.append(Flux)
                Wavelength  = float(line.split(' ')[7])
                Wavelength_A.append(Wavelength)
            except:
                print("did not work on line {}".format(lineCounter))

        lineCounter +=1


    f.close()

    peaks, _ = signal.find_peaks(Flux_erg_s_cm2_A,prominence=0.06*np.max(Flux_erg_s_cm2_A))

    return np.array(Wavelength_A), np.array(Flux_erg_s_cm2_A), peaks

WL_Spec1, Flux_Spec1, peak1 = ReadSpectrum("Spectrum1.txt")
WL_Spec2, Flux_Spec2, peak2 = ReadSpectrum("Spectrum2.txt")
WL_Spec3, Flux_Spec3, peak3 = ReadSpectrum("Spectrum3.txt")
WL_Spec4, Flux_Spec4, peak4 = ReadSpectrum("Spectrum4.txt")

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

y_formatter = mpl.ticker.ScalarFormatter(useOffset=False)
fs = 16

peaks, _ = signal.find_peaks(Flux_Spec1,prominence=(None, 0.5e-17))

def FluxCorrection(x,A,B,C,D):
    return A*x**3 + B*x**2 + C*x + D



Cor_Coefs, error    = curve_fit(FluxCorrection, xdata=WL_Spec1[peaks], ydata=Flux_Spec1[peaks])
WaveLength          = np.linspace(3500,7025,1000)

fig, ax = plt.subplots(1,dpi=200,sharex=True)
ax.plot(WL_Spec1,Flux_Spec1,color='blue',linewidth=1)
ax.plot(WL_Spec1[peaks],Flux_Spec1[peaks],'.r')
ax.plot(WaveLength,FluxCorrection(WaveLength,*Cor_Coefs),'yellow')
ax.annotate(r"$H_{\beta}$",xy=(4861,1.1e-15))
ax.annotate(r"$H_{\alpha}$",xy=(6564,2.6e-15))
ax.annotate(r"$OIII$",xy=(5007,4.4e-15))
ax.annotate(r"$OIII$",xy=(4958.92,2e-15))
ax.annotate(r"$OII$",xy=(3727,1.6e-15))
ax.annotate(r"$OI$",xy=(6300,1.5e-15))
ax.annotate(r"$NeIII$",xy=(3869,0.8e-15))
ax.annotate(r"$NII$",xy=(6548.03,2.0e-15))
ax.annotate(r"$NII$",xy=(6584,1.4e-15))
ax.annotate(r"Continuum emission",xy=(4500,0.2e-15))
ax.annotate(r"Spectrum 1",xy=(0.05,0.9),xycoords ='axes fraction',fontsize=fs-2)
ax.axhline(y=0,ls=':',color='black',linewidth=1)
ax.set_ylabel(r"Flux [erg s cm$^2$ $\AA$]",fontsize=fs)
ax.set_ylim(0,5e-15)
plt.close()

fig, ax = plt.subplots(2,1,dpi=200,sharex=True)
ax[0].plot(WL_Spec1,Flux_Spec1,color='blue',linewidth=1)
ax[0].plot(WL_Spec1[peak1],Flux_Spec1[peak1],'.r')
ax[0].plot(WaveLength,FluxCorrection(WaveLength,*Cor_Coefs),'purple',label=r"$\Phi(\lambda) = \beta_3 \lambda^3 + \beta_2 \lambda^2 + \beta_1 \lambda + \beta_0$")
ax[0].annotate(r"$H_{\beta}$",xy=(4861,1.1e-15))
ax[0].annotate(r"$H_{\alpha}$",xy=(6564,2.6e-15))
ax[0].annotate(r"$OIII$",xy=(5007,4.4e-15))
ax[0].annotate(r"$OIII$",xy=(4958.92,2e-15))
ax[0].annotate(r"$OII$",xy=(3727,1.6e-15))
ax[0].annotate(r"$OI$",xy=(6300,1.5e-15))
ax[0].annotate(r"$NeIII$",xy=(3869,0.8e-15))
ax[0].annotate(r"$NII$",xy=(6548.03,2.0e-15))
ax[0].annotate(r"$NII$",xy=(6584,1.4e-15))
#ax[0].annotate(r"Continuum emission",xy=(4500,0.2e-15))
ax[0].annotate(r"Spectrum 1",xy=(0.05,0.9),xycoords ='axes fraction',fontsize=fs-2)
ax[0].axhline(y=0,ls=':',color='black',linewidth=1)
ax[0].set_ylabel(r"$\Phi$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]",fontsize=fs)
ax[0].set_ylim(0,5e-15)
ax[0].legend(fontsize=12)
ax[1].plot(WL_Spec2,Flux_Spec2,color='blue',linewidth=1)
ax[1].axhline(y=np.mean(Flux_Spec2))
ax[1].plot(WL_Spec2[peak2],Flux_Spec2[peak2],'.r')
ax[1].annotate(r"$H_{\beta}$",xy=(4861,2.8e-10))
ax[1].annotate(r"$H_{\alpha}$",xy=(6564,7.8e-10))
ax[1].annotate(r"$OIII$",xy=(5007,5.4e-10))
ax[1].annotate(r"$OIII$",xy=(4958.92,1.9e-10))
ax[1].annotate(r"$OII$",xy=(3727,3.6e-10))
ax[1].annotate(r"$OII$",xy=(3729,1.1e-10))
ax[1].annotate(r"$NII$",xy=(6548.03,1.4e-10))
ax[1].annotate(r"$NII$",xy=(6584,4.3e-10))
ax[1].annotate(r"$H_{\gamma}$",xy=(4340,1.2e-10))
ax[1].annotate(r"Spectrum 2",xy=(0.05,0.9),xycoords ='axes fraction',fontsize=fs-2)
ax[1].set_ylabel(r"$\Phi$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]",fontsize=fs)
ax[1].set_xlabel(r"$\lambda$ [$\AA$]",fontsize=fs)
ax[1].set_ylim(0,8.5e-10)
plt.close()


fig, ax = plt.subplots(2,1,dpi=200,sharex=True)
ax[0].plot(WL_Spec3,Flux_Spec3,color='blue',linewidth=1)
ax[0].plot(WL_Spec3[peak3],Flux_Spec3[peak3],'.r')
ax[0].set_ylabel(r"$\Phi$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]",fontsize=fs)
ax[0].annotate(r"$H_{\beta}$",xy=(4861,5.1e-14))
ax[0].annotate(r"$H_{\alpha}$",xy=(6564,1.9e-13))
ax[0].annotate(r"$OIII$",xy=(5007,2.6e-14))
ax[0].annotate(r"$OIII$",xy=(4958.92,1.06e-14))
ax[0].annotate(r"$OII$",xy=(3727,3.9e-14))
ax[0].annotate(r"$NII$",xy=(6548.03,3.1e-14))
ax[0].annotate(r"$NII$",xy=(6584,8.7e-14))
ax[0].annotate(r"$H_{\gamma}$",xy=(4340,2.3e-14))
ax[0].annotate(r"Spectrum 3",xy=(0.05,0.9),xycoords ='axes fraction',fontsize=fs-2)
ax[0].set_ylim(0,2.15e-13)
ax[1].plot(WL_Spec4,Flux_Spec4,color='blue',linewidth=1)
ax[1].plot(WL_Spec4[peak4],Flux_Spec4[peak4],'.r')
ax[1].annotate(r"$H_{\beta}$",xy=(4861,3.8e-14))
ax[1].annotate(r"$H_{\alpha}$",xy=(6564,1.2e-13))
ax[1].annotate(r"$OIII$",xy=(5007,5.1e-14))
ax[1].annotate(r"$OIII$",xy=(4958.92,2.0e-14))
ax[1].annotate(r"$OII$",xy=(3727,5e-14))
ax[1].annotate(r"$NII$",xy=(6548.03,1.9e-14))
ax[1].annotate(r"$NII$",xy=(6584,4.8e-14))
ax[1].annotate(r"$H_{\gamma}$",xy=(4340,1.8e-14))
ax[1].annotate(r"Spectrum 4",xy=(0.05,0.9),xycoords ='axes fraction',fontsize=fs-2)
ax[1].set_ylim(0,1.3e-13)
ax[1].set_ylabel(r"$\Phi$ [erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$]",fontsize=fs)
ax[1].set_xlabel(r"$\lambda$ [$\AA$]",fontsize=fs)
plt.close()


def CorrectReddening(lambda1,lambda2,Hbeta,Halpha,C1,C2):

    Ratio = np.log10(lambda1/lambda2) + C1*np.log10(Halpha/Hbeta) - C2

    return Ratio

# First Spectra

def Redshift(obs,emit):
    return (obs - emit)/emit

# Spectrum 1
###############################################################################################

#5007/Hbeta
OIII5007    = Flux_Spec1[peak1[9]]  - FluxCorrection(WL_Spec1[peak1[9]],*Cor_Coefs)     # 5007
Redshift(WL_Spec1[peak1[9]],5007.0)
Hbeta       = Flux_Spec1[peak1[7]]  - FluxCorrection(WL_Spec1[peak1[7]],*Cor_Coefs)     # 4861
Redshift(WL_Spec1[peak1[7]],4861.0)
Halpha      = Flux_Spec1[peak1[18]] - FluxCorrection(WL_Spec1[peak1[18]],*Cor_Coefs)    # 6564
Redshift(WL_Spec1[peak1[18]],6564.0)
C1_01       = -0.11
C2_01       = -0.05
Line1Cor    = CorrectReddening(OIII5007,Hbeta,Hbeta,Halpha,C1_01,C2_01)

#3727/5007
OII3727     = Flux_Spec1[peak1[5]] - FluxCorrection(WL_Spec1[peak1[5]],*Cor_Coefs)        # 3727
Redshift(WL_Spec1[peak1[5]],3727.0)
C1_02       = 0.98
C2_02       = 0.45
Line2Cor_1  = CorrectReddening(OII3727,OIII5007,Hbeta,Halpha,C1_02,C2_02)
Brent_1     = CorrectReddening(OIII5007,OII3727,Hbeta,Halpha,C1_02,C2_02)

#6584/Halpha
NII6584     = Flux_Spec1[peak1[19]] - FluxCorrection(WL_Spec1[peak1[19]],*Cor_Coefs)      # 6584
Redshift(WL_Spec1[peak1[19]],6584.0)
C1_03       = -0.01
C2_03       = 0
Line3Cor_1  = CorrectReddening(NII6584,Halpha,Hbeta,Halpha,C1_03,C2_03)

#6300/Halpha
OI6300      = Flux_Spec1[peak1[16]] - FluxCorrection(WL_Spec1[peak1[16]],*Cor_Coefs)
Redshift(WL_Spec1[peak1[16]],6300.0)
C1_04       = 0.12
C2_04       = 0.05
Line4Cor_1    = CorrectReddening(OI6300,Halpha,Hbeta,Halpha,C1_04,C2_04)

# Delta Construction

#DeltaE (5007 / Hbeta)

DeltaE_5007_1 = Line1Cor + np.log10(0.32 + 10**Line2Cor_1) - 0.44
DeltaE_6584_1 = 0.5*(Line3Cor_1 - np.log10( 10**Line2Cor_1/(10**Line2Cor_1 + 1.93) ) + 0.37 )
DeltaE_6300_1 = 0.2*(Line4Cor_1 + 2.23)
DeltaE_1      = 0.33*(DeltaE_5007_1 + DeltaE_6584_1 + DeltaE_6300_1)

# Spectrum 2
###############################################################################################

#5007/Hbeta
OIII5007    = Flux_Spec2[peak2[6]]  - np.mean(Flux_Spec2)    # 5007
Redshift(WL_Spec2[peak2[6]],5007.0)
WL_Spec2[peak2[6]]
Hbeta       = Flux_Spec2[peak2[4]]  - np.mean(Flux_Spec2)    # 4861
Redshift(WL_Spec2[peak2[4]],4861.0)
WL_Spec2[peak2[4]]
Halpha      = Flux_Spec2[peak2[8]]  - np.mean(Flux_Spec2)    # 6564
Redshift(WL_Spec2[peak2[8]],6564.0)
WL_Spec2[peak2[8]]
C1_01       = -0.11
C2_01       = -0.05
Line1Cor    = CorrectReddening(OIII5007,Hbeta,Hbeta,Halpha,C1_01,C2_01)

#3727/5007
OII3727     = Flux_Spec2[peak2[1]] - np.mean(Flux_Spec2)        # 3727
Redshift(WL_Spec2[peak2[1]],3727.0)
WL_Spec2[peak2[1]]
C1_02       = 0.98
C2_02       = 0.45
Line2Cor_2  = CorrectReddening(OII3727,OIII5007,Hbeta,Halpha,C1_02,C2_02)
Brent_2     = CorrectReddening(OIII5007,OII3727,Hbeta,Halpha,C1_02,C2_02)

#6584/Halpha
NII6584     = Flux_Spec2[peak2[9]] - np.mean(Flux_Spec2)     # 6584
WL_Spec2[peak2[9]]
Redshift(WL_Spec2[peak2[9]],6584.0)
C1_03       = -0.01
C2_03       = 0
Line3Cor_2  = CorrectReddening(NII6584,Halpha,Hbeta,Halpha,C1_03,C2_03)


#6300/Halpha
OI6300      = 7.30e-12 - np.mean(Flux_Spec2)
Redshift(6302.05,6300.0)
C1_04       = 0.12
C2_04       = 0.05
Line4Cor_2    = CorrectReddening(OI6300,Halpha,Hbeta,Halpha,C1_04,C2_04)


DeltaE_5007_2 = Line1Cor + np.log10(0.32 + 10**Line2Cor_2) - 0.44
DeltaE_6584_2 = 0.5*(Line3Cor_2 - np.log10( 10**Line2Cor_2/(10**Line2Cor_2 + 1.93) ) + 0.37 )
DeltaE_6300_2 = 0.2*(Line4Cor_2 + 2.23)
DeltaE_2      = 0.33*(DeltaE_5007_2 + DeltaE_6584_2 + DeltaE_6300_2)
#
# Spectrum 3
###############################################################################################

#5007/Hbeta
OIII5007    = Flux_Spec3[peak3[3]]  - np.mean(Flux_Spec3)    # 5007
WL_Spec3[peak3[3]]
Redshift(WL_Spec3[peak3[3]],5007.0)
Hbeta       = Flux_Spec3[peak3[2]]  - np.mean(Flux_Spec3)    # 4861
WL_Spec3[peak3[2]]
Redshift(WL_Spec3[peak3[2]],4861.0)
Halpha      = Flux_Spec3[peak3[5]]  - np.mean(Flux_Spec3)    # 6564
Redshift(WL_Spec3[peak3[5]],6564.0)
WL_Spec3[peak3[5]]
C1_01       = -0.11
C2_01       = -0.05
Line1Cor    = CorrectReddening(OIII5007,Hbeta,Hbeta,Halpha,C1_01,C2_01)

#3727/5007
OII3727     = Flux_Spec3[peak3[0]] - np.mean(Flux_Spec3)        # 3727
Redshift(WL_Spec3[peak3[0]],3727.0)
WL_Spec3[peak3[0]]
C1_02       = 0.98
C2_02       = 0.45
Line2Cor_3  = CorrectReddening(OII3727,OIII5007,Hbeta,Halpha,C1_02,C2_02)
Brent_3     = CorrectReddening(OIII5007,OII3727,Hbeta,Halpha,C1_02,C2_02)

#6584/Halpha
NII6584     = Flux_Spec3[peak3[6]] - np.mean(Flux_Spec3)     # 6584
Redshift(WL_Spec3[peak3[6]],6584.0)
WL_Spec3[peak3[6]]
C1_03       = -0.01
C2_03       = 0
Line3Cor_3  = CorrectReddening(NII6584,Halpha,Hbeta,Halpha,C1_03,C2_03)

#6300/Halpha
OI6300      = 3.08e-15
#\lambda = 6323.97
Redshift(6323.97,6300.0)
C1_04       = 0.12
C2_04       = 0.05
Line4Cor_3  = CorrectReddening(OI6300,Halpha,Hbeta,Halpha,C1_04,C2_04)


DeltaE_5007_3 = Line1Cor + np.log10(0.32 + 10**Line2Cor_3) - 0.44
DeltaE_6584_3 = 0.5*(Line3Cor_3 - np.log10( 10**Line2Cor_3/(10**Line2Cor_3 + 1.93) ) + 0.37 )
DeltaE_6300_3 = 0.2*(Line4Cor_3 + 2.23)
DeltaE_3      = 0.33*(DeltaE_5007_3 + DeltaE_6584_3 + DeltaE_6300_3)
#
# Spectrum 4
###############################################################################################

#5007/Hbeta
OIII5007    = Flux_Spec4[peak4[5]]  - np.mean(Flux_Spec4)    # 5007
WL_Spec4[peak4[5]]
Redshift(WL_Spec4[peak4[5]],5007.0)
Hbeta       = Flux_Spec4[peak4[3]]  - np.mean(Flux_Spec4)    # 4861
Redshift(WL_Spec4[peak4[3]],4861.0)
WL_Spec4[peak4[3]]
Halpha      = Flux_Spec4[peak4[7]]  - np.mean(Flux_Spec4)    # 6564
WL_Spec4[peak4[7]]
Redshift(WL_Spec4[peak4[7]],6564.0)
C1_01       = -0.11
C2_01       = -0.05
Line1Cor    = CorrectReddening(OIII5007,Hbeta,Hbeta,Halpha,C1_01,C2_01)

#3727/5007
OII3727     = Flux_Spec4[peak4[0]] - np.mean(Flux_Spec4)        # 3727
WL_Spec4[peak4[0]]
Redshift(WL_Spec4[peak4[0]],3727.0)
C1_02       = 0.98
C2_02       = 0.45
Line2Cor_4  = CorrectReddening(OII3727,OIII5007,Hbeta,Halpha,C1_02,C2_02)
Brent_4     = CorrectReddening(OIII5007,OII3727,Hbeta,Halpha,C1_02,C2_02)

#6584/Halpha
NII6584     = Flux_Spec4[peak4[8]] - np.mean(Flux_Spec4)     # 6584
WL_Spec4[peak4[8]]
Redshift(WL_Spec4[peak4[8]],6584.0)
C1_03       = -0.01
C2_03       = 0
Line3Cor_4    = CorrectReddening(NII6584,Halpha,Hbeta,Halpha,C1_03,C2_03)

#6300/Halpha
OI6300      = 7.69e-15 - np.mean(Flux_Spec4)
#\lambda = 6366.8
Redshift(6366.8,6300.0)
C1_04       = 0.12
C2_04       = 0.05
Line4Cor_4  = CorrectReddening(OI6300,Halpha,Hbeta,Halpha,C1_04,C2_04)


DeltaE_5007_4 = Line1Cor + np.log10(0.32 + 10**Line2Cor_4) - 0.44
DeltaE_6584_4 = 0.5*(Line3Cor_4 - np.log10( 10**Line2Cor_4/(10**Line2Cor_4 + 1.93) ) + 0.37 )
DeltaE_6300_4 = 0.2*(Line4Cor_4 + 2.23)
DeltaE_4      = 0.33*(DeltaE_5007_4 + DeltaE_6584_4 + DeltaE_6300_4)


DeltaE      = np.array([DeltaE_1,DeltaE_2,DeltaE_3,DeltaE_4])
Line5007    = np.array([Line2Cor_1,Line2Cor_2,Line2Cor_3,Line2Cor_4])
Brent      = np.array([Brent_1,Brent_2,Brent_3,Brent_4])
Line6300    = np.array([Line4Cor_1,Line4Cor_2,Line4Cor_3,Line4Cor_4])
x = np.linspace(min(Line6300),max(Line6300),1000)
x2 = np.linspace(-1.23,max(Line6300),1000)
HIIcomp = -1.7*x - 2.63
LINERcomp = 1*x2 + 0.7


fig, ax = plt.subplots(1,2,dpi=200)
for counter in xrange(0,4):
    ax[0].plot(Line5007[counter],DeltaE[counter],'*',markersize=12)
    ax[0].annotate("S{}".format(counter+1),xy=(Line5007[counter],DeltaE[counter]+0.03))
ax[0].set_ylim(-0.1,0.55)
ax[0].annotate("Shock Heating",xy=(0.65,0.5),xycoords ='axes fraction',fontsize=12)
ax[0].annotate("Region",xy=(0.65,0.45),xycoords ='axes fraction',fontsize=12)
ax[0].annotate("Power-Law",xy=(0.28,0.52),xycoords ='axes fraction',fontsize=12)
ax[0].annotate("Photo-Ionisation",xy=(0.28,0.47),xycoords ='axes fraction',fontsize=12)
ax[0].annotate("Region",xy=(0.28,0.42),xycoords ='axes fraction',fontsize=12)
ax[0].annotate("Planetaries",xy=(0.05,0.95),xycoords ='axes fraction',fontsize=12,color="red")
ax[0].annotate("HII Regions",xy=(0.05,0.86),xycoords ='axes fraction',fontsize=12,color="red")
ax[0].annotate("Baldwin et al. 1981",xy=(0.01,0.01),xycoords ='axes fraction',fontsize=10)
ax[0].axhline(y=0.5,ls=':',color="red",linewidth=1)
ax[0].axvline(x=0,ls=':',color="black",linewidth=1)
ax[0].set_ylabel(r'$\left\langle \Delta E \right\rangle$',fontsize=16)
ax[0].set_xlabel('OII[3727]/OIII[5007]',fontsize=16)
for counter in xrange(0,4):
    ax[1].plot(Line6300[counter],Brent[counter],'*',markersize=12)
    ax[1].annotate("S{}".format(counter+1),xy=(Line6300[counter],Brent[counter]+0.08))
ax[1].plot(x,HIIcomp,color='black',ls=":",linewidth=1)
ax[1].plot(x2,LINERcomp,color='black',ls=':',linewidth=1)
ax[1].annotate("Seyferts",xy=(0.5,0.8),xycoords ='axes fraction',fontsize=12)
ax[1].annotate("HII regions",xy=(0.1,0.5),xycoords ='axes fraction',fontsize=12)
ax[1].annotate("LINERS",xy=(0.8,0.4),xycoords ='axes fraction',fontsize=12)
ax[1].annotate("Kewley et al. 2006",xy=(0.01,0.01),xycoords ='axes fraction',fontsize=10)
ax[1].set_ylabel('OIII[5007]/OII[3727]',fontsize=16)
ax[1].set_xlabel(r'OI[6300]/H$_{\alpha}$',fontsize=16)
plt.show()
