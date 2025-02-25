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
from numpy import trapz
# from scipy.integrate import simpson
from numpy.linalg import norm

########################################################################################################################
# SET PATHS

#Path to save figures
folder_to_save_figure = '../plots/'

#File folder
folder_data = '../'

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

fig,ax = plt.subplots(1,2)

#Colors and labels
colors=['#a6cee3','#1f78b4','#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6']
labels=['$Q_1Q_1$', '$Q_2Q_2$', '$Q_3Q_3$']
lethe_markers = ['s', 'd', 'o']


areas = []
l2_norms = []
linf_norms =[]

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

        # EXTRACT RELEVANT DATA
        name_files = ['/ke_rate.dat', '/viscous_dissipation.dat', '/torque.00.dat']
        t_1, ke_rate = np.loadtxt(folder_data + name_files[0], skiprows=10, unpack=True)
        t_2, visc = np.loadtxt(folder_data + name_files[1], skiprows=10, unpack=True)
        t_3, tx, ty, tz = np.loadtxt(folder_data + name_files[2], skiprows=10, unpack=True)
        size = np.min([len(t_1),len(t_2),len(t_3)])
        
        # Calculate balance
        scale=np.max(visc)
        torque_power = tz[:size]*omega/vol
        balance = -(-ke_rate[:size]+torque_power[:size]+visc[:size])/scale

        # Calculate area
        area = trapz(balance, dx=t_1[-1] - t_1[-2])
        # area2 = simpson(balance, dx=t_1[-1] - t_1[-2])

        # Append values to arrays
        # areas.append(area2)
        l2_norms.append(norm(balance, 2))
        linf_norms.append(norm(balance, np.inf))

    

# Number of DoFs and times extracted from simulations
dofs = [765760, 4368000, 27946240, 5642240, 33034240, 215982080, 18274240, 108443520, 716803840]
times = [782.8, 4976.2, 155983.6, 5157, 52057.4, 603477.8, 34982.1, 344198.4, 3625603.2]
times = [x / 3600 for x in times]

y_labels = [r'$L^2$' + "-norm" + r'$({\epsilon}_{\textrm{n}})$']
x_labels = ["No. DoFs", "Total time " + r'$\times$'+ " Nodes [h]"]
for j in range(2):
    if(j == 0):
        x = dofs
        y = l2_norms
    elif(j == 1):
        x = times
        y = l2_norms

    ax[j].plot(x[0:3], y[0:3], color=colors[0], marker=lethe_markers[0], lw=4, markeredgecolor='k', markeredgewidth=2, markersize=15, zorder=3)
    ax[j].plot(x[3:6], y[3:6], color=colors[1], marker=lethe_markers[1], lw=4, markeredgecolor='k', markeredgewidth=2, markersize=15, zorder=2, linestyle = 'dashed')
    ax[j].plot(x[6:9], y[6:9], color=colors[2], marker=lethe_markers[2], lw=4, markeredgecolor='k', markeredgewidth=2, markersize=15, zorder=1)
    ax[j].set_xlabel(x_labels[j])
    ax[j].set_ylabel(y_labels[0])
    ax[j].grid(True)
    ax[j].set_xscale('log')
    if (j == 0):
        ax[j].set_xticks([10**6, 10**7, 10**8, 10**9])
        ax[j].add_patch(plt.Rectangle((7*10**6,0.9),7*10**7,0.35,
                                  fill=False, color='#fb9a99', alpha=1, zorder=0, linewidth=4))      
    elif (j == 1):
        ax[j].set_xticks([10**-1, 10**0, 10**1, 10**2, 10**3])
        ax[j].add_patch(plt.Rectangle((5*10**0,0.9),80,0.35,
                                  fill=False, color='#fb9a99', alpha=1, zorder=0, linewidth=4))
    ax[j].get_xaxis().set_minor_locator(plt.LogLocator(base=10, subs='all', numticks=100))

fig.set_size_inches(15,8)
fig.tight_layout()

fig.subplots_adjust(top=0.8)

legend_elements = [Line2D([0], [0], color=colors[0], lw=4, marker=lethe_markers[0], markeredgecolor='k', markeredgewidth=2, markersize=15),
                    Line2D([0], [0], color=colors[1], lw=4, marker=lethe_markers[1], markeredgecolor='k', markeredgewidth=2, markersize=15, linestyle = 'dashed'),
                    Line2D([0], [0], color=colors[2], lw=4, marker=lethe_markers[2], markeredgecolor='k', markeredgewidth=2, markersize=15)]


n_columns = 3

fig.legend(legend_elements, labels, loc='upper center', ncol = n_columns,  facecolor = 'white', framealpha = 0.75,  edgecolor = 'black', fancybox = False, shadow = False, fontsize=38, handlelength=1.5, columnspacing=1)

plt.savefig(folder_to_save_figure + 'energy_balance_l2_norm_vs_time_and_dofs.pdf')

