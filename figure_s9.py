import matplotlib.pyplot as plt
import os
from utils import read_force_displacement_data
import numpy as np

'''
THREE COUNTERSNAPPING ELEMENTS COUPLED IN SERIES EXPERIMENT (FIGURE S9, PANELS C,D)
'''
colors = ['#025196', '#4C8E73', '#A0AE54', '#fdb338']
t, u, f = read_force_displacement_data(os.path.join('FIGURE S9', 'Specimen_RawData_1.csv'))
branch_start_indices = [3001, 5441, 6702, 8238, 9443, 9904, 11693, 11825, 13022, 13479, 13854]
branch_end_indices = [5416, 6688, 8231, 9432, 9880, 11679, 11816, 13011, 13459, 13826, 14335]
branch_cols = [colors[0], colors[3], colors[2], colors[3], colors[2], colors[1], colors[2], colors[3], colors[2],
               colors[1], colors[0]]


# plot
_, axs = plt.subplots(2, 1)
for i in range(len(branch_start_indices)):
    # 300 data points per force-displacement branch, otherwise figure to heavy
    indices = np.round(np.linspace(branch_start_indices[i], branch_end_indices[i], 300)).astype(int)
    axs[1].plot(u[indices], f[indices], 'o', color=branch_cols[i], markersize=2)
    if i < len(branch_start_indices) - 1:
        axs[1].plot(u[branch_end_indices[i]:branch_start_indices[i + 1]], f[branch_end_indices[i]:branch_start_indices[i + 1]],
                'o', color='#aaaaaa', markersize=1)
axs[1].set_xlabel('displacement (mm)')
axs[1].set_ylabel('force (N)')

indices = np.round(np.linspace(branch_start_indices[0], branch_end_indices[-1], 2000)).astype(int)
axs[0].plot(t[indices] - t[branch_start_indices[0]], u[indices],
            'o', markersize=2, lw=.5, color='#000000', zorder=4)
for i in range(len(branch_start_indices)):
    axs[0].axvspan(t[branch_start_indices[i]] - t[branch_start_indices[0]],
                    t[branch_end_indices[i]] - t[branch_start_indices[0]],
                    facecolor=branch_cols[i], alpha=0.4)
for i in range(len(branch_start_indices) - 1):
    axs[0].axvspan(t[branch_end_indices[i]] - t[branch_start_indices[0]],
                    t[branch_start_indices[i + 1]] - t[branch_start_indices[0]],
                    facecolor='#aaaaaa')
axs[0].set_ylabel('displacement (mm)')
axs[0].set_xlabel('time (s)')

plt.show()