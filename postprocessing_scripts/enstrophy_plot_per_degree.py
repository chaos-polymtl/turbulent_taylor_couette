# SPDX-FileCopyrightText: Copyright (c) 2025 The Lethe Authors
# SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later

import numpy as np
import matplotlib.pyplot as plt
from numpy import loadtxt
import pandas as pd
import math
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

########################################################################################################################
# SET PATHS

#Path to save figures
folder_to_save_figure = '../plots/'

#File folder
folder_data = '../'

#Path to reference data
path_to_reference = '../reference_data/fig12_wang2021.dat'

# READ DATA

#Read reference data
reference_values = np.loadtxt(path_to_reference, skiprows = 1, unpack = True, delimiter = ',')

#Read lethe data
t_ens = []
ens = []

rho = 1.0
omega = 1.0
ri = 0.5
d = 0.5
U = omega*ri

coeff = (2*d*d)/(rho*U*U)

# PLOT

#Settings
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['font.size'] = 38

fig,ax = plt.subplots(1,3)

#For zoom plots
sub_axes=list()
for i in range(3):
    sub_axes.append(inset_axes(ax[i], height=1.8, width="40%", loc='upper left'))

#Colors and labels
colors=['#a6cee3','#1f78b4','#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6']
titles=['$Q_1Q_1$', '$Q_2Q_2$', '$Q_3Q_3$']
lethe_markers = ['^', 's', 'o']
markers_frequence=[70, 140]
markers_frequence_zoom=[140, 280]

#Plot relevant values
for j in range(3):
    if(j == 0):
        sim = ['Q1Q1_l4', 'Q1Q1_l5', 'Q1Q1_l6']
    elif(j == 1):
        sim = ['Q2Q2_l4', 'Q2Q2_l5', 'Q2Q2_l6']
    elif(j == 2):
        sim = ['Q3Q3_l4', 'Q3Q3_l5', 'Q3Q3_l6']

    for i in range(len(sim)):

        #File folder
        folder_data = '../' + str(sim[i]) 
        print(folder_data)

        # EXTRACT RELEVANT DATA
        temp_t_ens,temp_ens=np.loadtxt(folder_data + '/enstrophy.dat', skiprows=10, unpack=True)
        temp_ens=temp_ens*coeff
        ax[j].plot(temp_t_ens, temp_ens, color=colors[i], lw=4)
        sub_axes[j].plot(temp_t_ens, temp_ens, color=colors[i], lw=4)

    ax[j].plot(reference_values[0], reference_values[8], ':', color='black', lw=4)
    sub_axes[j].plot(reference_values[0], reference_values[8], ':', color='black', lw=4)

#Adjust output settings
for i in range(3):
    ax[i].set_xlabel('Time [s]')
    ax[i].set_ylabel('Enstrophy [-]')
    ax[i].set_xlim([0, 33.5])
    ax[i].set_ylim([0, 18])
    ax[i].set_title(titles[i])
    ax[i].set_xticks([0,10,20,30])
    sub_axes[i].set_xticks([])
    sub_axes[i].set_yticks([])
    sub_axes[i].axis([19,30,10,18])

fig.set_size_inches(21,8)
fig.tight_layout()

fig.subplots_adjust(top=0.75)

#Legend
legend_elements = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4),
                    Line2D([0], [0], color='k', lw=4, linestyle=':')]

labels = ['153K cells', '942K cells', '6M cells', 'Wang DNS reference']

fig.legend(legend_elements, labels, loc='upper center', ncol = 4,  facecolor = 'white', framealpha = 0.75,  edgecolor = 'black', fancybox = False, shadow = False, fontsize=38, handlelength=1.5, columnspacing=1)

#Save fig
plt.savefig(folder_to_save_figure + 'enstrophy_comparison_per_degree.pdf')
# plt.show()