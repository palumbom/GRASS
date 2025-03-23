import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime
from astropy.time import Time
from barycorrpy import get_BC_vel

grass_data_combined = h5py.File("projected_SSD_gpu_2026_combined.jld2", "r")
GRASS_rv_combined  = grass_data_combined["RV_list_no_cb"][()]
time_extended  = grass_data_combined["time"][()]
phase_angle  = grass_data_combined["phase_angle"][()]

grass_data_CB = h5py.File("projected_SSD_gpu_2026_CB.jld2", "r")
GRASS_rv_CB  = grass_data_CB["RV_list_no_cb"][()]

grass_data_SH = h5py.File("projected_SSD_gpu_2026_SH.jld2", "r")
GRASS_rv_SH  = grass_data_SH["RV_list_no_cb"][()]

UTC_time_ex = []
time_julian_ext = []
for i in time_extended:
    dt = datetime.strptime(i.decode("utf-8"), "%Y-%m-%dT%H:%M:%S.%f")
    UTC_time_ex.append(dt)
    time_julian_ext.append((Time(dt)).jd)

vb_ext, warnings, flag = get_BC_vel(JDUTC=time_julian_ext, lat=31.9583 , longi=-111.5967, alt=209.7938, SolSystemTarget='Sun', predictive=False,zmeas=0.0)

def jld2_read(jld2_file, variable, v, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array + v)
    array -= array[-1]
    return array

for i in range(0, 1):
    
    GRASS_rv_array_SH = jld2_read(grass_data_SH, GRASS_rv_SH, vb_ext, i)
    GRASS_rv_array_CB = jld2_read(grass_data_CB, GRASS_rv_CB, vb_ext, i)
    GRASS_rv_array_combined = jld2_read(grass_data_combined, GRASS_rv_combined, vb_ext, i)
    phase_angle_array = grass_data_combined[phase_angle[i]][()]

    # rm curve 
    fig, ax = plt.subplots()
    ax.plot(UTC_time_ex, GRASS_rv_array_combined, color = 'r', linewidth = 2, label = "Combined")
    ax.plot(UTC_time_ex, GRASS_rv_array_SH, color = 'b', linewidth = 2, label = "SH")
    ax.plot(UTC_time_ex, GRASS_rv_array_CB, color = 'g', linewidth = 2, label = "CB")
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.set_xlabel("Time on 01/10/2026 (UTC)", fontsize=12)
    ax.set_ylabel("RV [m/s]", fontsize=12)
    ax.legend()
    plt.savefig("2026.png")
    plt.clf()

    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].plot(UTC_time_ex, GRASS_rv_array_combined, color = 'r', linewidth = 2, label = "Combined")
    axs[0].plot(UTC_time_ex, GRASS_rv_array_SH, color = 'b', linewidth = 2, label = "SH")
    axs[0].plot(UTC_time_ex, GRASS_rv_array_CB, color = 'g', linewidth = 2, label = "CB")
    axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    start_date = pd.Timestamp('2026-01-10T08:00:00')
    end_date = pd.Timestamp('2026-01-10T11:30:00')
    axs[0].axvspan(start_date, end_date, color='gray', alpha=0.3)
    axs[1].axvspan(start_date, end_date, color='gray', alpha=0.3)
    axs[0].set_xlabel("Time on 01/10/2026 (UTC)", fontsize=12)
    axs[0].set_ylabel("RV [m/s]", fontsize=12)
    axs[0].legend(fontsize=12)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    #residuals
    axs[1].scatter(UTC_time_ex, phase_angle_array, color = 'r', marker = "x", s = 3)  
    axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axs[1].set_xlabel("Time on 01/10/2026 (UTC)", fontsize=12)
    axs[1].set_ylabel("Phase Angle (degrees)", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("2026_phase_angle")
    plt.clf()