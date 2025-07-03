import matplotlib.pyplot as plt
import os
from utils import read_force_displacement_data
import numpy as np

'''
TWO COUNTERSNAPPING ELEMENTS COUPLED IN PARALLEL EXPERIMENT (FIGURE 5, PANELS C,D,G)
'''
data_dir = os.path.join('FIGURE 5', 'parallel_coupling')
colors = ['#025196', '#71ac62', '#fdb338']
t, u, f = read_force_displacement_data(os.path.join(data_dir, 'Specimen_RawData_1.csv'))

branch_start_indices = [27390, 32463, 35427, 39638, 39830, 43127, 45991, 49309, 49433]
branch_end_indices = [32455, 35412, 39631, 39821, 43110, 45982, 49295, 49424, 51776]
branch_colors = [colors[0], colors[1], colors[0], colors[1], colors[2], colors[1], colors[2], colors[1], colors[0]]

# plot
_, axs = plt.subplots(2, 1)

indices = np.round(np.linspace(branch_start_indices[0], branch_end_indices[-1], 2000)).astype(int)
axs[0].plot(t[indices] - t[branch_start_indices[0]], u[indices], 'o', markersize=2, lw=.75, color='#000000', zorder=4)

for i in range(len(branch_start_indices)):
    axs[0].axvspan(t[branch_start_indices[i]] - t[branch_start_indices[0]],
                    t[branch_end_indices[i]] - t[branch_start_indices[0]],
                    facecolor=branch_colors[i], alpha=0.4)

for i in range(len(branch_start_indices) - 1):
    axs[0].axvspan(t[branch_end_indices[i]] - t[branch_start_indices[0]],
                    t[branch_start_indices[i + 1]] - t[branch_start_indices[0]],
                    facecolor='#aaaaaa')
axs[0].set_xlabel('time (s)')
axs[0].set_ylabel('displacement (mm)')



for i in range(len(branch_start_indices)):
    # only 400 data points are plotted per force-displacement branch otherwise, plot is too heavy
    indices = np.round(np.linspace(branch_start_indices[i], branch_end_indices[i], 400)).astype(int)
    axs[1].plot(u[indices], f[indices], 'o', color=branch_colors[i], markersize=2)
    if i < len(branch_start_indices) - 1:
        axs[1].plot(u[branch_end_indices[i]:branch_start_indices[i + 1]],
                f[branch_end_indices[i]:branch_start_indices[i + 1]],
                'o', color='#aaaaaa', markersize=2)
axs[1].set_ylabel('force (N)')
axs[1].set_xlabel('displacement (mm)')
plt.show()


'''
TWO COUNTERSNAPPING ELEMENTS COUPLED IN SERIES EXPERIMENT (FIGURE 5, PANELS J,K,N)
'''
data_dir = os.path.join('FIGURE 5', 'series_coupling')
colors = ['#025196', '#71ac62', '#fdb338']
t, u, f = read_force_displacement_data(os.path.join(data_dir, 'Specimen_RawData_1.csv'))

branch_start_indices = [2807, 5358, 6748, 8328, 9755, 10344]
branch_end_indices = [5345, 6731, 8319, 9740, 10325, 11023]
branch_colors = [colors[0], colors[2], colors[1], colors[2], colors[1], colors[0]]

# plot
_, axs = plt.subplots(2, 1)

indices = np.round(np.linspace(branch_start_indices[0], branch_end_indices[-1], 2000)).astype(int)
axs[0].plot(t[indices] - t[branch_start_indices[0]], u[indices], 'o', markersize=2, lw=.75, color='#000000', zorder=4)

for i in range(len(branch_start_indices)):
    axs[0].axvspan(t[branch_start_indices[i]] - t[branch_start_indices[0]],
                    t[branch_end_indices[i]] - t[branch_start_indices[0]],
                    facecolor=branch_colors[i], alpha=0.4)

for i in range(len(branch_start_indices) - 1):
    axs[0].axvspan(t[branch_end_indices[i]] - t[branch_start_indices[0]],
                    t[branch_start_indices[i + 1]] - t[branch_start_indices[0]],
                    facecolor='#aaaaaa')
axs[0].set_xlabel('time (s)')
axs[0].set_ylabel('displacement (mm)')



for i in range(len(branch_start_indices)):
    # only 400 data points are plotted per force-displacement branch otherwise, plot is too heavy
    indices = np.round(np.linspace(branch_start_indices[i], branch_end_indices[i], 400)).astype(int)
    axs[1].plot(u[indices], f[indices], 'o', color=branch_colors[i], markersize=2)
    if i < len(branch_start_indices) - 1:
        axs[1].plot(u[branch_end_indices[i]:branch_start_indices[i + 1]],
                f[branch_end_indices[i]:branch_start_indices[i + 1]],
                'o', color='#aaaaaa', markersize=2)
axs[1].set_ylabel('force (N)')
axs[1].set_xlabel('displacement (mm)')
plt.show()

