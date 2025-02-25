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
path_to_reference = '../reference_data/fig3a_wang2021.dat'

# READ DATA

#Read reference data
reference_values = np.loadtxt(path_to_reference, skiprows = 1, unpack = True, delimiter = ',')


#Read lethe data
t_ke = []
ke = []

rho = 1.0
omega = 1.0
ri = 0.5
d = 0.5
U = omega*ri

coeff = (2)/(U*U)

name_ke_file = '/kinetic_energy.dat'
files = ['Q1Q1_l4', 'Q1Q1_l5', 'Q1Q1_l6', 'Q2Q2_l4', 'Q2Q2_l5', 'Q2Q2_l6', 'Q3Q3_l4', 'Q3Q3_l5', 'Q3Q3_l6']

for index,filename in enumerate(files):
    temp_t_ke,temp_ke=np.loadtxt(folder_data + filename + name_ke_file ,skiprows=1,unpack=True)
    temp_ke=temp_ke*coeff
    t_ke.append(temp_t_ke)
    ke.append(temp_ke)

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
    sub_axes.append(inset_axes(ax[i], height=1.8, width="40%", loc='lower left'))

#Colors and labels
colors=['#a6cee3','#1f78b4','#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6']
titles=['$Q_1Q_1$', '$Q_2Q_2$', '$Q_3Q_3$']
lethe_markers = ['^', 's', 'o']
markers_frequence=[70, 140]
markers_frequence_zoom=[140, 280]

#Plot relevant values
count = 0
for j in range(3):
    for i in range(3):
        index = i + (j*3)
        ax[j].plot(t_ke[index], ke[index], color=colors[i], lw=4)
        count = count + 1
    ax[j].plot(reference_values[0], reference_values[1], ':', color='black', lw=4)


#Plot same values in sub axes
count = 0
for j in range(3):
    for i in range(3):
        index = i + (j*3)
        sub_axes[j].plot(t_ke[index], ke[index], color=colors[i], lw=4)
        count = count + 1
    sub_axes[j].plot(reference_values[0], reference_values[1], ':', color='black', lw=4)

#Adjust output settings
for i in range(3):
    ax[i].set_xlabel('Time [s]')
    ax[i].set_ylabel('Kinetic energy [-]')
    ax[i].set_xlim([0, 33.5])
    ax[i].set_ylim([0.12, 0.22])
    ax[i].set_title(titles[i])
    ax[i].set_xticks([0,10,20,30])

    ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
    ax[i].yaxis.get_offset_text().set_visible(False)
    ax[i].text(0, 0.225, r'$\times 10^{-1}$')

    sub_axes[i].set_xticks([])
    sub_axes[i].set_yticks([])
    sub_axes[i].axis([18,30,0.15,0.2])

fig.set_size_inches(21,8)
fig.tight_layout()

fig.subplots_adjust(top=0.75)

#Legend
legend_elements = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4),
                    Line2D([0], [0], color='k', lw=4, linestyle='dashed')]

labels = ['153K', '942K', '6M', 'Wang DNS reference']

fig.legend(legend_elements, labels, loc='upper center', ncol = 4,  facecolor = 'white', framealpha = 0.75,  edgecolor = 'black', fancybox = False, shadow = False, fontsize=38, handlelength=1.5, columnspacing=1)

#Save fig
plt.savefig(folder_to_save_figure + 'kinetic_energy_comparison_per_degree.pdf')
# plt.show()