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

#Read lethe data
t_ens = []
ens = []

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

#For zoom plots
loc_sub_axes = ['upper left', 'lower left', 'lower left']
sub_axes=list()
for i in range(2):
    sub_axes.append(inset_axes(ax[i], height=1.8, width="40%", loc=loc_sub_axes[i]))

#Colors and labels
colors=['#a6cee3','#1f78b4','#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6']
ylabels=['Enstrophy [-]', "Kinetic energy [-]", r'${\epsilon}_{\textrm{n}}$' +' [-]']
lethe_markers = ['^', 's', 'o']
max_ylim = 0
min_ylim = 1000

#Plot relevant values
for j in range(3):
    if(j == 0):
        field_name = '/enstrophy.dat'
        data_index = 8
        path_to_reference = '../reference_data/fig12_wang2021.dat'
        coeff = (2*d*d)/(rho*U*U)
    elif(j == 1):
        field_name = '/kinetic_energy.dat'
        data_index = 1
        path_to_reference = '../reference_data/fig3a_wang2021.dat'
        coeff = (2)/(U*U)
    
    sim = ['Q2Q2_l5', 'Q2Q2_l5_cfl_2', 'Q2Q2_l5_cfl_4']

    for i in range(len(sim)):

        #File folder
        folder_data = '../' + sim[i] 
        print(folder_data)

        # EXTRACT RELEVANT DATA
        if  (j != 2):
            t,value=np.loadtxt(folder_data + field_name, skiprows=10, unpack=True)
            value=value*coeff
            ax[j].plot(t, value, color=colors[i], lw=4, zorder = 0)
            sub_axes[j].plot(t, value, color=colors[i], lw=4)
        
            #Read reference data
            reference_values = np.loadtxt(path_to_reference, skiprows = 1, unpack = True, delimiter = ',')
            ax[j].plot(reference_values[0], reference_values[data_index], ':', color='black', lw=4, alpha=0.3, zorder = 1)

            sub_axes[j].plot(reference_values[0], reference_values[data_index], ':', color='black', lw=4)
        
        # Balance
        if(j == 2):
            name_files = ['/ke_rate.dat', '/viscous_dissipation.dat', '/torque.00.dat']
            t_1, ke_rate = np.loadtxt(folder_data + name_files[0], skiprows=2, unpack=True)
            t_2, visc = np.loadtxt(folder_data + name_files[1], skiprows=2, unpack=True)
            t_3, tx, ty, tz = np.loadtxt(folder_data + name_files[2], skiprows=2, unpack=True)
            size = np.min([len(t_1),len(t_2),len(t_3)])
            scale=np.max(visc)
            torque_power = tz[:size]*omega/vol
            balance = -(-ke_rate[:size]+torque_power[:size]+visc[:size])/scale
            ax[j].plot(t_1[:size], balance, label="Balance",color=colors[i], lw=3, zorder=2)

            max_ylim = max(max_ylim, max(balance))
            min_ylim = min(min_ylim, min(balance))


#Adjust output settings
for i in range(3):
    ax[i].set_xlabel('Time [s]')
    ax[i].set_ylabel(ylabels[i])
    ax[i].set_xticks([0,10,20,30])
    ax[i].set_xlim([0, 33.5])

    if (i != 2):
        sub_axes[i].set_xticks([])
        sub_axes[i].set_yticks([])
    if(i == 0):
        ax[i].set_ylim([0, 18])
        sub_axes[i].axis([19,30,10,18])
    elif(i == 1):
        ax[i].set_ylim([0.12, 0.22])
        sub_axes[i].axis([18,30,0.15,0.2])
        ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax[i].yaxis.get_offset_text().set_visible(False)
        ax[i].text(0, 0.225, r'$\times 10^{-1}$')
    elif(i == 2):
        ax[i].set_xticks([0,10,20,33.5])
        ax[i].set_ylim([min_ylim + (0.1*min_ylim), max_ylim +(0.1*max_ylim)])
        ax[i].ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        ax[i].yaxis.get_offset_text().set_visible(False)
        ax[i].text(0, 0.105 , r'$\times 10^{-1}$')

fig.set_size_inches(21,8)
fig.tight_layout()

fig.subplots_adjust(top=0.75)

#Legend
legend_elements = [Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4),
                    Line2D([0], [0], color='k', lw=4, linestyle=':', alpha=0.7)]

labels = ['CFL = 1', 'CFL = 2', 'CFL = 4', 'Wang DNS reference']

fig.legend(legend_elements, labels, loc='upper center', ncol = 4,  facecolor = 'white', framealpha = 0.75,  edgecolor = 'black', fancybox = False, shadow = False, fontsize=38, handlelength=1.5, columnspacing=1)

#Save fig
plt.savefig(folder_to_save_figure + 'sim_5_cfl_comparison.pdf')
# plt.show()