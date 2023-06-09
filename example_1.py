#%% Recommended to run this script in Jupyter // VSCode (Free)
# Imports!
%matplotlib inline
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import math
from scipy.signal import argrelextrema
import crystalspacing #local llibrary 

#%%
#1 Import data.
filename = "example_1-tbig-ggg-bi2te3.asc"
path = os.getcwd() + "\\" + filename
xrd_file = open(path)

# %%
# 2 File IO
def xrrxrd_filereader(file):
    """This function reads ascii files (.asc) and converts the data to a readable set of keys."""
    params = {}
    data = []
    for line in file:
        if line[0] == "*":
            parts = line.split("=")
            if len(parts) > 0 and type(parts) == list:
                params[parts[0].strip()] = [a.strip() for a in parts[1:]]
        else:
            if len(line.strip()) > 0:
                data += [float(a.strip())
                         for a in line.replace("\n", "").split(",")]
    return np.array(data), params

xr_data, xr_params = xrrxrd_filereader(xrd_file)
xrd_file.close()


xr_START = float(xr_params["*START"][0])
xr_STOP = float(xr_params["*STOP"][0])
xr_COUNT = int(xr_params["*COUNT"][0])
xr_STEP = float(xr_params["*STEP"][0])

# calculate 1D array of angles
xr_theta = np.linspace(xr_START, xr_STOP, xr_COUNT)

#%%
#3 Generate XRD peaks of materials according to Bragg law.
# n * theta = 2 d sin(theta).
# lattice constants defined below. Need to use various crytal planes to identify peaks.
spacing_hex = crystalspacing.spacing.hex
spacing_cubic = crystalspacing.spacing.cubic


crystallography = {
    "bi2te3" : ([04.38, 04.38, 30.50], spacing_hex), #lattice angstroms, hexagonal spacing. #Space group: R -3 m H   (166)
    # "al2o3"  : ([4.7607, 4.7607, 12.9947], spacing_hex), #Space group: R-3 c H (167)
    # "y3fe5o12" : ([12.376, 12.376, 12.376 ], spacing_cubic),
    "tb3fe5o12": ([12.435, 12.435, 12.435 ], spacing_cubic),
    "gd3ga5o12": ([12.3772, 12.3772, 12.3772], spacing_cubic),
}


#%%
#Set default plotting rcParams
acsNano_Shared_Params = {
    'figure.dpi':300,
    "font.family":
        'arial',
    'font.size': 7, #4.5 Minimum!
    "axes.titlesize":8,
    "axes.labelsize":8,
    "xtick.labelsize":7,
    "ytick.labelsize":7,
    "legend.fontsize":7,
    "legend.framealpha":0,
    "axes.titlesize":7,
    "axes.linewidth":0.5,
}
acsNano_2Col_Params = acsNano_Shared_Params.copy()
acsNano_2Col_Params.update({
    "figure.figsize" : [7,3.5], #4.167in (300pt) to 7in (504pt) wide maximum, 9.167 in (660pt) high.
})
plt.rcParams.update(acsNano_2Col_Params)
#%%
# Plot the raw data
xr_labels = {"xlabel":r"$\theta$","ylabel":r"Counts (A.U.)","title":r"$\theta/2\theta$"}
fig = plt.semilogy(xr_theta, xr_data,label="Raw XRR Data")[0].get_figure()
fig.axes[0].set_xlabel(xr_labels["xlabel"])
fig.axes[0].set_ylabel(xr_labels["ylabel"])
fig.axes[0].set_title(xr_labels["title"])

#%%
# Identify Peaks in data

#smooth data using savgol
# xr_data_smooth = savgol_filter(xr_data, window_length=31, polyorder=4) #Note this gives negative values!!

#smooth data using gaussian convolution
# sigma = 0.055 #degrees
sigma = 0.055 #degrees
xx = np.arange(-3*sigma, 3*sigma, xr_STEP)
#NOTE THIS WONT WORK (Len(xr_theta_smooth)) IF LEN(xx) IS EVEN.

def gaussian(x):
    return np.exp(-(x/sigma)**2/2)
gx = gaussian(xx)

xr_data_smooth = np.convolve(xr_data, gx, mode="valid") / (2 * math.pi)
halfw = math.floor((len(xx))/2)
xr_theta_smooth = xr_theta[halfw-1:-halfw-1]

fig=plt.semilogy(xr_theta_smooth,xr_data_smooth)[0].get_figure()
#take maxima
smooth_order = 30
maxima=argrelextrema(xr_data_smooth, np.greater, order=smooth_order)[0]
outofbounds = (maxima > smooth_order) & (maxima < len(xr_data_smooth) - smooth_order)
maxima = maxima[outofbounds]

#Check that maxima can be seen locally in the data.
valid_maxima = []
req_flof = 0.95
req_flof_w = 40

for m in maxima:
    xr_local = xr_data_smooth[m-req_flof_w:m+req_flof_w]
    # xr_local
    s1 = np.average(np.log10(xr_local[:req_flof_w]))
    s2 = np.average(np.log10(xr_local[req_flof_w:]))
    # s1
    # s2
    # np.log10(xr_data_smooth[m])
    #Check peak value can be identified as more than 5% of log value.
    peak = np.log10(xr_data_smooth[m])
    if abs(s1) < req_flof*peak or abs(s2) < req_flof*peak:
        #Valid peak
        valid_maxima.append(m)
len(valid_maxima)

# xr_theta_smooth = np.convolve(xr_theta, gaussian(gx), mode="valid") / (2 * math.pi)
# xr_theta_smooth = np.convolve(xr_theta, gaussian(gx), mode="full") / (2 * math.pi) + (len(gx)/2 + 1.5) * xr_STEP

# Check that data before and after match
i1 = np.where(xr_data == max(xr_data))[0][0]
i2 = np.where(xr_data_smooth == max(xr_data_smooth))[0][0]
xr_theta[i1]
xr_theta_smooth[i2]


for m in valid_maxima:
    fig.axes[0].plot([xr_theta_smooth[m],xr_theta_smooth[m]],[5e1,1e7], linewidth=0.5)
fig


#Import Crystals from Materials Project!
# MATPROJ_API_KEY = "hsmr2YLOxLu9a09rP88JKgrFgsUEQcIu"
# from mp_api.client import MPRester
# mpr = MPRester(MATPROJ_API_KEY)
#
# ###Search for materials:
# #graphene
# graphene_lattices = mpr.summary.search(elements=["C"], crystal_system="Hexagonal", num_elements=1, band_gap=(0,0), fields=["material_id", "structure", "efermi", "nsites"])
# len(graphene_lattices)
# graphene_lattices[2]
#
# graphene_lattice = mpr.get_structure_by_material_id("mp-568806")
# graphene_lattice


#calculate miller indicies spacing (d).
peaks = {}
for mat in crystallography:
    mat_peaks = {}
    for i in range(-15,15):
        for j in range(-15,15):
            for k in range(-15,15):
                if not i+j+k == 0: #cannot be all zero
                    label = "" + str(i) + "," + str(j) + "," + str(k)
                    d = crystallography[mat][1](crystallography[mat][0],(i,j,k))
                    mat_peaks[label] = d
    peaks[mat] = mat_peaks


# Try to match peaks to condition n \lambda = 2d sin(\theta).
# Copper X-Ray Source:
lambda_alpha1 = float(xr_params["*WAVE_LENGTH1"][0]) # angstrom^-1
lambda_alpha2 = float(xr_params["*WAVE_LENGTH2"][0]) # angstrom^-1
lambda_alpha1
lambda_alpha2
lambda_beta = 1.3922
lambas = [lambda_alpha1, lambda_alpha2]

xr_theta_smooth[valid_maxima]

matched_peaks = {}
for m in valid_maxima:
    m_th = xr_theta_smooth[m] / 2

    closest = (None, None, None) #Mat, Indexes, Spacing.

    #Get First Itme
    init_Mat = list(crystallography)[0]
    init_Ind = list(peaks[init_Mat])[0]
    init_D = peaks[init_Mat][init_Ind]
    closest = (init_Mat,init_Ind,0,init_D,0) #Material, Miller Index, Bragg N, Peak proximity,

    for mat in crystallography:
        for peak in peaks[mat]:
            d = peaks[mat][peak]
            for n in range(15):
                try:
                    peak2 = 2*math.asin(n*lambda_alpha1/(2*d)) * 180 / math.pi
                except ValueError:
                    peak2 = math.nan

                bragg = abs((n*lambda_alpha1 - 2*d*math.sin(m_th * math.pi / 180))) # n * theta = 2 d sin(theta).
                if  bragg < closest[3]:
                    #New Closest!
                    closest = (mat, peak, n, bragg, peak2)
    matched_peaks[np.around(xr_theta_smooth[m],decimals=2)] = closest
pd.DataFrame(matched_peaks, ["Material","Miller Indx.", "Bragg n", "Δd","θ"])


# %%
