import matplotlib.pyplot as plt
import os
import numpy as np

'''
FREE VIBRATION EXPERIMENT (FIGURE 4, PANELS C,E)
'''
data_dir = os.path.join('FIGURE 4', 'free_vibrations')
px_per_mm = 6.67  # see info.txt
fps = 50.0  # see info.txt

pull_release_indices = np.loadtxt(os.path.join(data_dir, 'pull_release_interval_frame_indices.csv')).astype(int)
switch_indices = np.loadtxt(os.path.join(data_dir, 'state_switch_frame_indices.csv')).astype(int)
switching_phases_indices = np.loadtxt(os.path.join(data_dir, 'switching_phases_interval_frame_indices.csv')).astype(int)
position = -np.loadtxt(os.path.join(data_dir, 'aruko_position_in_px.csv'), delimiter=',') / px_per_mm

start_index = round(10 * fps)  # skip 10 first seconds of video where nothing significant happens
displacement = position - position[start_index]
t = np.arange(displacement.shape[0]) / fps
t -= t[start_index]

# plot
_, ax = plt.subplots()
ax.plot(t[start_index:], displacement[start_index:], '-', color='#000000', markersize=2, lw=0.75, zorder=4)

for i in range(switch_indices.shape[0]):
    if i % 2 == 0:
        if i >= switch_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = switch_indices[i + 1]
        ax.axvspan(t[switch_indices[i]], t[end_index], facecolor='#ccdcea', label='state 0' if i == 0 else '')
    else:
        if i >= switch_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = switch_indices[i + 1]
        ax.axvspan(t[switch_indices[i]], t[end_index], facecolor='#fee8c3', label='state 1' if i == 1 else '')

for i in range(0, pull_release_indices.shape[0], 2):
    end_index = pull_release_indices[i + 1]
    ax.axvspan(t[pull_release_indices[i]], t[end_index],facecolor='#ffffff', alpha=0.5)

# dummy for legend
ax.axvspan(t[1], t[0], facecolor='#FFFFFF', alpha=0.5, edgecolor='#BBBBBB', label='pull-release')

for i in range(0, switching_phases_indices.shape[0], 2):
    end_index = switching_phases_indices[i + 1]
    ax.axvspan(t[switching_phases_indices[i]], t[end_index], facecolor='#000000', alpha=0.2, label='manual switch' if i == 0 else '')

ax.set_xlabel('time (s)')
ax.set_ylabel('change in extension (mm)')
ax.set_xlim((t[start_index], 32.0))
ax.legend(numpoints=5, loc='lower left', ncol=4, bbox_to_anchor=(0, 0.925))
plt.show()

'''
FORCED VIBRATIONS SOFT TO STIFF EXPERIMENT, FIGURE 4, PANEL G-TOP
'''
data_dir = os.path.join('FIGURE 4', 'forced_vibrations_soft_to_stiff')

px_per_mm = 8.18
fps = 50.0
start_index = round(7.50 * fps)
switch_indices = np.loadtxt(os.path.join(data_dir, 'switch_frame_indices.csv'), delimiter=',').astype(int)
uy = np.loadtxt(os.path.join(data_dir, 'tracking_results', 'uy.csv'), delimiter=',') / px_per_mm
ux = np.loadtxt(os.path.join(data_dir, 'tracking_results', 'ux.csv'), delimiter=',') / px_per_mm
t = np.arange(uy.shape[0]) / fps
t -= t[start_index]

# plot
_, ax = plt.subplots()
ax.plot(t[start_index:], uy[start_index:, 1], '-o', color='#888888', markersize=2, label='output')
ax.plot(t[start_index:], uy[start_index:, 0], '-o', color='#000000', markersize=2, label='input')

for i in range(switch_indices.shape[0]):
    if i % 2 == 0:
        if i >= switch_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = switch_indices[i + 1]
        ax.axvspan(t[switch_indices[i]], t[end_index],
                        facecolor='#ccdcea', zorder=0.0)
    else:
        if i >= switch_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = switch_indices[i + 1]
        ax.axvspan(t[switch_indices[i]], t[end_index],
                        facecolor='#fee8c3', zorder=0.0)
ax.legend(numpoints=5, loc='upper left')
ax.set_ylim((-11, 11))
ax.set_xlim(t[start_index], 13.5)
ax.set_ylabel('I/O (mm)')
ax.set_xlabel('time (s)')
plt.show()

'''
FORCED VIBRATIONS STIFF TO SOFT EXPERIMENT, FIGURE 4, PANEL G-BOTTOM
'''
data_dir = os.path.join('FIGURE 4', 'forced_vibrations_stiff_to_soft')

px_per_mm = 8.18
fps = 50.0
start_index = round(6.50 * fps)
switch_indices = np.loadtxt(os.path.join(data_dir, 'switch_frame_indices.csv'), delimiter=',').astype(int)
uy = np.loadtxt(os.path.join(data_dir, 'tracking_results', 'uy.csv'), delimiter=',') / px_per_mm
ux = np.loadtxt(os.path.join(data_dir, 'tracking_results', 'ux.csv'), delimiter=',') / px_per_mm
t = np.arange(uy.shape[0]) / fps
t -= t[start_index]

# plot
_, ax = plt.subplots()
ax.plot(t[start_index:], uy[start_index:, 1], '-o', color='#888888', markersize=2, label='output')
ax.plot(t[start_index:], uy[start_index:, 0], '-o', color='#000000', markersize=2, label='input')

for i in range(switch_indices.shape[0]):
    if i % 2 == 0:
        if i >= switch_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = switch_indices[i + 1]
        ax.axvspan(t[switch_indices[i]], t[end_index],
                        facecolor='#ccdcea', zorder=0.0)
    else:
        if i >= switch_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = switch_indices[i + 1]
        ax.axvspan(t[switch_indices[i]], t[end_index],
                        facecolor='#fee8c3', zorder=0.0)
ax.legend(numpoints=5, loc='upper left')
ax.set_ylim((-11, 11))
ax.set_xlim(t[start_index], 8.2)
ax.set_ylabel('I/O (mm)')
ax.set_xlabel('time (s)')
plt.show()
