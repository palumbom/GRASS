import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
mpl = plt.matplotlib 
import matplotlib.dates as mdates
from datetime import datetime, timedelta
from astropy.time import Time
from barycorrpy import get_BC_vel, exposure_meter_BC_vel

# model
file = h5py.File("neid_all_lines_rv_off_SSD_gpu.jld2", "r")
RV_list_no_cb = file["rv"][()]

def jld2_read(jld2_file, variable, index):
    array = jld2_file[variable[index]][()]
    array = np.array(array)
    array -= array[-1]
    return array

def plot_line(model, model_label):
    fig, axs = plt.subplots(2, sharex=True, sharey=False, gridspec_kw={'hspace': 0, 'height_ratios': [3, 1]})
    axs[0].plot(range(0, len(model)), model, color = 'r', linewidth = 2, label = model_label)
    axs[0].scatter(range(0, len(model)), model, color = 'b', linewidth = 2, label = model_label)
    # axs[0].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[0].set_xlabel("Time (UTC)", fontsize=12)
    # axs[0].set_ylabel("RV [m/s]", fontsize=12)
    # axs[0].legend(fontsize=12)
    # rms_model_no_cb = round(np.sqrt((np.nansum((line_rv_array - model)**2))/len(line_rv_array - model)),2)
    # axs[0].text(UTC_time[-40], 800, "Model RMS {}".format(rms_model_no_cb))
    # rms_grass_no_cb = round(np.sqrt((np.nansum((line_rv_array - GRASS)**2))/len(line_rv_array - GRASS)),2)
    # axs[0].text(UTC_time[-40], 700, "GRASS RMS {}".format(rms_grass_no_cb))
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    #residuals
    # axs[1].scatter(UTC_time, line_rv_array - model, color = 'r', marker = "x", s = 3) 
    # axs[1].scatter(UTC_time, line_rv_array - GRASS, color = 'b', marker = "x", s = 3)  
    # axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    # axs[1].set_xlabel("Time (UTC)", fontsize=12)
    # axs[1].set_ylabel("Residuals", fontsize=12) 
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig("test_gpu_projected.png")
    plt.clf()


RV_list_no_cb_array = jld2_read(file, RV_list_no_cb, 0)

    # rm curve 
plot_line(RV_list_no_cb_array, "Weighted RVs - No CB")