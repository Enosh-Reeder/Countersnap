import matplotlib.pyplot as plt
import os
from utils import read_force_displacement_data, extract_loading_sequence, extract_unloading_sequence
import numpy as np
from scipy.stats import linregress
from scipy.interpolate import interp1d

'''
DISPLACEMENT DRIVEN EXPERIMENT (FIGURE 2, PANELS A,C)
'''
normalized = False  # set to True to plot the normalized data instead
data_dir = os.path.join('FIGURE 2', 'displacement_driven')
t, u, f = read_force_displacement_data(os.path.join(data_dir, 'Specimen_RawData_1.csv'))

start_l, end_l = extract_loading_sequence(u, 1)
start_u, end_u = extract_unloading_sequence(u, 1)

_, ax = plt.subplots()
if normalized:
    uc = 25.2
    fc = 0.877
else:
    uc = fc = 1.0
ax.plot(u[start_l:end_l]/uc, f[start_l:end_l]/fc, 'o', markersize=1, color='#000000')
ax.plot(u[start_u:end_u]/uc, f[start_u:end_u]/fc, 'o', markersize=1, color='#999999')
ax.set_xlabel(f'{"normalized " if normalized else ""}displacement {"(mm)" if not normalized else "(-)"}')
ax.set_ylabel(f'{"normalized " if normalized else ""}force {"(N)" if not normalized else "(-)"}')
plt.show()



'''
FORCE DRIVEN EXPERIMENT (FIGURE 2, PANEL E)
'''
data_dir = os.path.join('FIGURE 2', 'force_driven')
delay = 7.21  # video started before instron 7.21 s before measurement started
px_per_mm = 6.68  # see info.txt
initial_u = 8.2934  # see info.txt
initial_f = 0.459285  # see info.txt
fps = 50  # see info.txt
initial_phase_frame_index = 1000  # see info.txt
initial_phase_instron_index = 1000  # see info.txt
linreg_start_instron_index = 7000  # see info.txt
linreg_end_instron_index = 7800  # see info.txt
correction_start_instron_index = 7300  # see info.txt
correction_end_instron_index = 7400  # see info.txt
dynamic_start_instron_index = 7331  # see info.txt
dynamic_end_instron_index = 7373  # see info.txt

# FORCE DATA FROM INSTRON
t_instron, _, f_instron = read_force_displacement_data(os.path.join(data_dir, 'Specimen_RawData_1.csv'))

# initial force measured by instron != initial weight applied --> re-calibration
f_calibrated = f_instron - np.median(f_instron[:initial_phase_instron_index]) + initial_f

# correcting the f_calibrated vector for the data point measured during the dynamic transition
# (which measured the forced felt by the instron sensor != applied force during the dynamic snapping)
# to estimate the applied force (weight of the amount of water) --> f == absolute applied force
slope, intercept, rvalue, _, _ = linregress(t_instron[linreg_start_instron_index:linreg_end_instron_index],
                                            f_calibrated[linreg_start_instron_index:linreg_end_instron_index])
f = f_calibrated.copy()
f[correction_start_instron_index:correction_end_instron_index] = slope * t_instron[correction_start_instron_index:correction_end_instron_index] + intercept

# ELONGATION DATA FROM VIDEO
position_video = np.loadtxt(os.path.join(data_dir, 'end_node_position_in_px', 'y.csv'), delimiter=',') / px_per_mm
initial_position = np.mean(position_video[:initial_phase_frame_index])

u_video = (position_video - initial_position) + initial_u
t_video = np.arange(u_video.shape[0]) / fps
u = interp1d(t_video - delay, u_video)(t_instron)

uc = 23.92  # critical displacement
fc = 0.809  # critical force

_, ax  = plt.subplots()
ax.plot(u[:dynamic_start_instron_index] / uc, f[:dynamic_start_instron_index] / fc,
        color='#000000', ls='', marker='o', markersize=2)
ax.plot(u[dynamic_start_instron_index:dynamic_end_instron_index] / uc,
        f[dynamic_start_instron_index:dynamic_end_instron_index] / fc,
        color='#000000', alpha=0.3, ls='',marker='o', markersize=2)
ax.plot(u[dynamic_end_instron_index:] / uc,
        f[dynamic_end_instron_index:] / fc,
        color='#000000', ls='',marker='o', markersize=2)
ax.set_xlabel('normalized displacement (-)')
ax.set_ylabel('normalized force (-)')
plt.show()

'''
MIXED DRIVEN EXPERIMENT (FIGURE 2, PANEL H)
'''
data_dir = os.path.join('FIGURE 2', 'mixed_driven')

px_per_mm = 155 / 24
fps = 50
delay = 0.608

weight_position_raw = np.loadtxt(os.path.join(data_dir, 'tracking_results_weight', 'y.csv'), delimiter=',') / px_per_mm
t_instron, u_instron, _ = read_force_displacement_data(os.path.join(data_dir, 'Specimen_RawData_1.csv'))

# calibrating for the second loading cycle
before_second_cs_data_index = 100831  # just before snapping (used for plotting with different colors)
after_second_cs_data_index = 100915  # just after snapping (used for plotting with different colors)
start, end = extract_loading_sequence(u_instron, 1)
t_video = np.arange(weight_position_raw.shape[0]) / fps

# normalized weight elevation and instron elevation
uc = 23.47  # critical_instron_elevation
weight_position = interp1d(t_video - delay, -weight_position_raw)(t_instron)
weight_elevation = weight_position - np.median(weight_position[before_second_cs_data_index - 4000:before_second_cs_data_index])
weight_elevation /= uc
u_instron /= uc

# plot
_, ax = plt.subplots()
ax.plot(u_instron[start:before_second_cs_data_index],
        weight_elevation[start:before_second_cs_data_index],
        'o', color='#000000', markersize=1, markeredgewidth=0.5)

ax.plot(u_instron[before_second_cs_data_index:after_second_cs_data_index],
        weight_elevation[before_second_cs_data_index:after_second_cs_data_index],
        'o', color='#999999', markersize=1, markeredgewidth=0.5)

ax.plot(u_instron[after_second_cs_data_index:end],
        weight_elevation[after_second_cs_data_index:end],
        'o', color='#000000', markersize=1, markeredgewidth=0.5)

ax.set_xlim((22.50 / uc, 25 / uc))
ax.set_ylim((-.01, 0.175))
ax.set_xlabel('normalized weight elevation (-)')
ax.set_ylabel('normalized instron elevation (-)')
plt.show()


