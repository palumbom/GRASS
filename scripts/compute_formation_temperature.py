"""
AUTHOR : Khaled Al Moulla
DATE   : 2023-03-16

Compute line formation temperature with PySME.
"""

#%% MODULES
import pdb
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pysme import sme as SME
from pysme.abund import Abund
from pysme.synthesize import Synthesizer
from pysme.linelist.vald import ValdFile
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d

#%% DATA
data = pd.read_csv("/Users/michael/Desktop/test_spec.csv")

# Wavelength and flux
wave  = data.wave.to_numpy()[0:1100] # add wavelength
flux  = data.flux.to_numpy()[0:1100] # add flux
Nwave = len(wave)

#%% SYNTHESIS | WITHOUT INSTRUMENTAL RESOLUTION & ROTATIONAL BROADENING

# Initiate SME
sme0 = SME.SME_Structure()
sme0.abund = Abund.solar()
sme0.linelist = ValdFile('Sun_VALD.lin')
sme0.atmo.source = 'marcs2012p_t1.0.sav'

# Parameters
sme0.teff , sme0.logg, sme0.monh = 5770, 4.00, 0.00

# Wavelength and flux
sme0.wave = wave
sme0.spec = flux

# Synthesize
synthesizer0 = Synthesizer()
sme0 = synthesizer0.synthesize_spectrum(sme0)

# Save
sme0.save('Sun_SME_0.npy')

#%% SYNTHESIS | WITH INSTRUMENTAL RESOLUTION & ROTATIONAL BROADENING

# Initiate SME
smeR = SME.SME_Structure()
smeR.abund = Abund.solar()
smeR.linelist = ValdFile('Sun_VALD.lin')
smeR.atmo.source = 'marcs2012p_t1.0.sav'

# Parameters
smeR.teff , smeR.logg, smeR.monh = 5777, 4.40, 0.00
smeR.vsini, smeR.vmic, smeR.vmac = 1.63, 0.85, 3.98

# Resolution
smeR.ipres, smeR.iptype = 700000, 'gauss'

# Wavelength and flux
smeR.wave = wave
smeR.spec = flux

# Solve
synthesizerR = Synthesizer()
smeR = synthesizerR.synthesize_spectrum(smeR)

# Save
smeR.save('Sun_SME_R.npy')

#%% FORMATION TEMPERATURE

# Temperature and density
temp = smeR.atmo.temp
rhox = smeR.atmo.rhox

# Empty array for T_1/2
T1o2_spec = np.empty(Nwave)
flux_spec = np.empty(Nwave)

# Loop wavelength points
for i in range(Nwave):
    
    # Opacity and source function
    lineop, contop, _, source, _ = synthesizerR.dll.GetLineOpacity(wave[i])
    totop = lineop + contop
    
    # Optical depth
    tau0 = rhox[0]*totop[0]/2
    tauw = np.append(tau0, tau0 + np.cumsum((rhox[1:]-rhox[:-1])*(totop[:-1]+totop[1:])/2))
    
    # Contribution function
    cont = np.cumsum(source*np.exp(-tauw))

    # get flux by integrating contribution function
    flux_spec[i] = trapezoid(cont)

    # normalize to max
    cont /= max(cont)
    
    # Interpolate temperature at 50%
    if (min(cont) < 1/2) & (max(cont) > 1/2):
        temp_func = interp1d(cont, temp)
        T1o2      = temp_func(1/2)
    else:
        T1o2      = np.nan
    
    # Store
    T1o2_spec[i] = T1o2

# Save
# np.savetxt('T_1/2_spec.txt', T1o2_spec)

pdb.set_trace()
