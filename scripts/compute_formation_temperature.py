"""
ORIGINAL AUTHOR : Khaled Al Moulla
DATE   : 2023-03-16

MODIFIED BY: Michael Palumbo

Compute line formation temperature with PySME.
"""

#%% MODULES
import os
import pdb
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pysme import sme as SME
from pysme.abund import Abund
from pysme.synthesize import Synthesizer
from pysme.linelist.vald import ValdFile
from scipy.integrate import trapezoid
from scipy.interpolate import interp1d
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import matplotlib.cm as cmx

#%% DATA
fpath, _ = os.path.split(__file__)
dpath = os.path.abspath(os.path.join(fpath, "..", "data"))
dfile = os.path.abspath(os.path.join(dpath, "line_info.csv"))
data = pd.read_csv(dfile)
line_names = data.name
airwavs = data.air_wavelength
avgtemp80 = np.zeros_like(airwavs)
avgtemp50 = np.zeros_like(airwavs)

spectra_files = glob.glob("/Users/michael/Desktop/mean_spec/*.csv")

# loop over wavelengths in lines
for i in range(len(airwavs)):
    # find the right data
    try:
        file_idx = np.where([(line_names[i] in file) for file in spectra_files])[0][0]
    except:
        continue
    the_file = spectra_files[file_idx]
    spec_data = pd.read_csv(the_file)

    # get size of resolution elements
    res = 7e5
    delta_wav = (airwavs[i]/res) / 2.0

    # Wavelength and flux
    buff  = 0.5
    wave  = spec_data.wavs.values # np.arange(airwavs[i] - buff, airwavs[i] + buff, step=delta_wav)
    flux  = spec_data.flux.values # np.zeros_like(wave)
    Nwave = len(wave)

    #%% SYNTHESIS | WITHOUT INSTRUMENTAL RESOLUTION & ROTATIONAL BROADENING

    # # Initiate SME
    # sme0 = SME.SME_Structure()
    # sme0.abund = Abund.solar()
    # sme0.linelist = ValdFile(os.path.join(dpath, "Sun_VALD.lin"))
    # sme0.atmo.source = 'marcs2012p_t1.0.sav'

    # # Parameters
    # sme0.teff , sme0.logg, sme0.monh = 5770, 4.00, 0.00

    # # Wavelength and flux
    # sme0.wave = wave
    # sme0.spec = flux

    # # Synthesize
    # synthesizer0 = Synthesizer()
    # sme0 = synthesizer0.synthesize_spectrum(sme0)

    # Save
    # sme0.save('Sun_SME_0.npy')

    #%% SYNTHESIS | WITH INSTRUMENTAL RESOLUTION & ROTATIONAL BROADENING

    # Initiate SME
    smeR = SME.SME_Structure()
    smeR.abund = Abund.solar()
    smeR.linelist = ValdFile(os.path.join(dpath, "Sun_VALD.lin"))
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
    # smeR.save('Sun_SME_R.npy')

    #%% FORMATION TEMPERATURE

    # Temperature and density
    temp = smeR.atmo.temp
    rhox = smeR.atmo.rhox

    # Empty array for T_1/2
    T1o2_spec = np.empty(Nwave)
    flux_spec = np.empty(Nwave)

    # Loop wavelength points
    for j in range(Nwave):

        # Opacity and source function
        lineop, contop, _, source, _ = synthesizerR.dll.GetLineOpacity(wave[j])
        totop = lineop + contop

        # Optical depth
        tau0 = rhox[0]*totop[0]/2
        tauw = np.append(tau0, tau0 + np.cumsum((rhox[1:]-rhox[:-1])*(totop[:-1]+totop[1:])/2))

        # Contribution function
        cont = np.cumsum(source*np.exp(-tauw))

        # get flux by integrating contribution function
        flux_spec[j] = trapezoid(cont)

        # normalize to max
        cont /= max(cont)

        # Interpolate temperature at 50%
        if (min(cont) < 1/2) & (max(cont) > 1/2):
            temp_func = interp1d(cont, temp)
            T1o2      = temp_func(1/2)
        else:
            T1o2      = np.nan

        # Store
        T1o2_spec[j] = T1o2

    # get normalized spectrum
    flux_spec_norm = flux_spec / np.max(flux_spec)

    # measure average formation temp
    wav_idx = np.searchsorted(wave, airwavs[i])
    min_idx = np.argmin(flux_spec_norm[wav_idx - 5:wav_idx+5]) + (wav_idx - 5)

    # plot it
    cmap = plt.get_cmap("inferno")
    cnorm = colors.Normalize(vmin=np.min(T1o2_spec), vmax=np.max(T1o2_spec))
    scalarMap = cmx.ScalarMappable(norm=cnorm, cmap=cmap)

    cols = np.linspace(0,1,len(wave))
    points = np.array([wave, flux_spec_norm]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    fig, ax1 = plt.subplots()
    lc = LineCollection(segments, cmap=cmap, norm=cnorm)
    lc.set_array(T1o2_spec)
    lc.set_linewidth(2)
    line = ax1.add_collection(lc)
    fig.colorbar(line,ax=ax1)
    ax1.axvline(wave[min_idx], c="k", ls="--")
    ax1.set_xlabel("Wavelength")
    ax1.set_ylabel("Normalized Flux")
    ax1.set_title(data.name[i])
    ax1.set_xlim(np.min(wave)-0.1, np.max(wave)+0.1)
    ax1.set_ylim(-0.1, 1.1)
    fig.savefig("/Users/michael/Desktop/" + data.name[i] + ".pdf")
    plt.clf(); plt.close()

    # get wings up to 20% depth
    line_depth = 1.0 - flux_spec_norm[min_idx]
    wing_depth80 = 0.2 * line_depth
    idxl80 = min_idx - np.searchsorted(flux_spec_norm[min_idx:-1:1], 1.0-wing_depth80)
    idxr80 = np.searchsorted(flux_spec_norm[min_idx:], 1.0-wing_depth80) + min_idx

    wing_depth50 = 0.5 * line_depth
    idxl50 = min_idx - np.searchsorted(flux_spec_norm[min_idx:-1:1], 1.0-wing_depth50)
    idxr50 = np.searchsorted(flux_spec_norm[min_idx:], 1.0-wing_depth50) + min_idx


    # get average temp for line
    avgtemp80[i] = np.mean(T1o2_spec[idxl80:idxr80])
    avgtemp50[i] = np.mean(T1o2_spec[idxl50:idxr50])

# write it out to the csv
data["avg_temp_80"] = avgtemp80
data["avg_temp_50"] = avgtemp50
data.to_csv(dfile)
