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

#With reference enstrophy superposed
enstrophy = True

# READ DATA

rho = 1.0
omega = 1.0
re = 1.0
ri = 0.5
d = 0.5
U = omega*ri
vol =np.pi**2*(re**2 - ri**2)

# PLOT

#Settings
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family']='DejaVu Serif'
plt.rcParams['font.serif']='cm'
plt.rcParams['font.size'] = 38

fig,ax = plt.subplots(1,3)

#Colors and labels
colors=['#a6cee3','#1f78b4','#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6']
titles=['$Q_1Q_1$', '$Q_2Q_2$', '$Q_3Q_3$']
lethe_markers = ['^', 's', 'o']
markers_frequence=[70, 140]
markers_frequence_zoom=[140, 280]

#Read reference data
reference_values = np.loadtxt(path_to_reference, skiprows = 1, unpack = True, delimiter = ',')

max_ylim = 0
min_ylim = 1000

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
        folder_data = '../' + sim[i] 
        print(folder_data)

        # EXTRACT RELEVANT DATA
        name_files = ['/ke_rate.dat', '/viscous_dissipation.dat', '/torque.00.dat']
        t_1, ke_rate = np.loadtxt(folder_data + name_files[0], skiprows=10, unpack=True)
        t_2, visc = np.loadtxt(folder_data + name_files[1], skiprows=10, unpack=True)
        t_3, tx, ty, tz = np.loadtxt(folder_data + name_files[2], skiprows=10, unpack=True)
        size = np.min([len(t_1),len(t_2),len(t_3)])
        
        # Plot enstrophy behind
        t_6, e = np.loadtxt(folder_data + '/enstrophy.dat', skiprows=10, unpack=True)
        e = e * (2*d*d)/(rho*U*U)
        ax_two=ax[j].twinx()
        # ax_two.plot(t_6, e, ':', label="Enstrophy", color=colors[i], lw=3, alpha=0.7, zorder=1)
        ax_two.tick_params(axis='y', labelcolor=colors[1])
        ax_two.set_yticks([])

        #Plot reference enstrophy
        if (enstrophy):
            ax_two.plot(reference_values[0], reference_values[8], ':', color='black', lw=3, zorder=2, alpha=0.3)

        # Calculate balance:
        scale=np.max(visc)
        torque_power = tz[:size]*omega/vol
        balance = -(-ke_rate[:size]+torque_power[:size]+visc[:size])/scale
        ax[j].plot(t_1[:size], balance, label="Balance",color=colors[i], lw=4, zorder=0)

        max_ylim = max(max_ylim, max(balance))
        min_ylim = min(min_ylim, min(balance))

        #Print values for each mesh:
        print("Energy Balance")
        print("Sim ", sim[i])
        print("Max y value: ", max(balance))
        print("Min y value: ", min(balance))

        #Reference vertical peak lines
        # if (enstrophy):
        #     times = [22.058, 24.016, 27.2, 33.5]
        #     for time in times:
        #         ax[j].axvline(x = time, color= colors[3], alpha =0.2, lw=1)

print("Overall maximum and minimum:")
print("Max: ", max_ylim)
print("Min: ", min_ylim)

#Adjust output settings
for i in range(3):
    ax[i].set_xlabel('Time [s]')
    ax[i].set_ylabel(r'${\epsilon}_{\textrm{n}}$' +' [-]')
    ax[i].set_xlim([0, 33.5])
    ax[i].set_xticks([0,10,20,33.5])
    ax[i].set_ylim([min_ylim + (0.1*min_ylim), max_ylim +(0.1*max_ylim)])
    ax[i].set_title(titles[i])
    ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax[i].yaxis.get_offset_text().set_visible(False)
    ax[i].text(0, 0.29, r'$\times 10^{-1}$')

fig.set_size_inches(21,8)
fig.tight_layout()

fig.subplots_adjust(top=0.75)

#Legend
if (enstrophy):
    legend_elements = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4),
                    Line2D([0], [0], color='black', lw=4, linestyle=':', alpha=0.7)]

    labels = ['153K cells', '942K cells', '6M cells', 'Wang DNS reference enstrophy']

    n_columns = 4
else:
    legend_elements = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4)]

    labels = ['153K cells', '942K cells', '6M cells']  

    n_columns = 3

fig.legend(legend_elements, labels, loc='upper center', ncol = n_columns,  facecolor = 'white', framealpha = 0.75,  edgecolor = 'black', fancybox = False, shadow = False, fontsize=38, handlelength=1.5, columnspacing=1)

#Save fig
if (enstrophy):
   plt.savefig(folder_to_save_figure + 'energy_balance_per_degree_with_enstrophy.pdf')
else:
   plt.savefig(folder_to_save_figure + 'energy_balance_per_degree.pdf') 
# plt.show()