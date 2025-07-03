from springable.readwrite.fileio import read_behavior
from springable.mechanics.mechanical_behavior import BezierBehavior
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import os
import matplotlib.cm as mcm
import matplotlib.colors as colors


'''
ESTIMATED FORCE-DISPLACEMENT CURVES OF THE INDIVIDUAL BUILDING BLOCKS (FIGURE S5, PANELS A,B,C,D,E)
'''
block_library_dir = os.path.join('FIGURE S2', 'block_library')
_, axs = plt.subplots(2, 3)

# softening f-d curves as a function of a (panel A)
block_folders = ['soft_0', 'soft_1', 'soft_2']
behaviors_paths = [os.path.join(block_library_dir, block_folder, 'bezier_behavior.csv')
                   for block_folder in block_folders]
a = [4, 5, 6]
_soft_behaviors = [read_behavior(behavior_path) for behavior_path in behaviors_paths]
control_point_sets = [soft_behavior.get_control_points() for soft_behavior in _soft_behaviors]
u1 = interp1d(a, [cp[0][1] for cp in control_point_sets], fill_value='extrapolate')
u2 = interp1d(a, [cp[0][2] for cp in control_point_sets], fill_value='extrapolate')
u3 = interp1d(a, [cp[0][3] for cp in control_point_sets], fill_value='extrapolate')
f1 = interp1d(a, [cp[1][1] for cp in control_point_sets], fill_value='extrapolate')
f2 = interp1d(a, [cp[1][2] for cp in control_point_sets], fill_value='extrapolate')
f3 = interp1d(a, [cp[1][3] for cp in control_point_sets], fill_value='extrapolate')
uu = np.linspace(0, 25, 1000)
a_vector = np.arange(3, 7.5, 0.5)
cmap = 'GnBu'
cn = colors.Normalize(vmin=1, vmax=8)  # 50% / 25%
scm = mcm.ScalarMappable(norm=cn, cmap=cmap)
for a_i in a_vector:
    extrapolation = a_i < 4 or a_i > 6
    soft_behavior = BezierBehavior(1.0, u_i=[u1(a_i), u2(a_i), u3(a_i)], f_i=[f1(a_i), f2(a_i), f3(a_i)])
    axs[0, 0].plot(uu, soft_behavior.gradient_energy(uu+1.0)[0],
                   color=scm.to_rgba(a_i), lw=1, ls='--' if extrapolation else '-')
for _soft_behavior in _soft_behaviors:
    axs[0, 0].plot(uu, _soft_behavior.gradient_energy(uu+1.0)[0], color='k', lw=0.5)
axs[0, 0].set_ylabel('force (N)')
axs[0, 0].set_xlabel('displacement (mm)')
axs[0, 0].set_xlim((-1.0, 20))
axs[0, 0].set_ylim((-0.1, 1.5))
cbar = plt.colorbar(scm, cax=None, ax=axs[0, 0])
cbar.ax.set_title('a (mm)')
cbar.ax.set_ylim((3, 7))

# softening f-d curves as function of b (panel B)
block_folders = ['soft_0', 'soft_3', 'soft_6']
behaviors_paths = [os.path.join(block_library_dir, block_folder, 'bezier_behavior.csv')
                   for block_folder in block_folders]
b = [1, 1.5, 2]
u1 = interp1d(b, [cp[0][1] for cp in control_point_sets], fill_value='extrapolate')
u2 = interp1d(b, [cp[0][2] for cp in control_point_sets], fill_value='extrapolate')
u3 = interp1d(b, [cp[0][3] for cp in control_point_sets], fill_value='extrapolate')
f1 = interp1d(b, [cp[1][1] for cp in control_point_sets], fill_value='extrapolate')
f2 = interp1d(b, [cp[1][2] for cp in control_point_sets], fill_value='extrapolate')
f3 = interp1d(b, [cp[1][3] for cp in control_point_sets], fill_value='extrapolate')
uu = np.linspace(0, 25, 1000)
b_vector = np.arange(0.5, 2.75, 0.25)
cmap = 'GnBu'
cn = colors.Normalize(vmin=-0.5, vmax=3)  # 50% / 25%
scm = mcm.ScalarMappable(norm=cn, cmap=cmap)
for b_i in b_vector:
    extrapolation = b_i > 2 or b_i < 1
    soft_behavior = BezierBehavior(1.0, u_i=[u1(b_i), u2(b_i), u3(b_i)], f_i=[f1(b_i), f2(b_i), f3(b_i)])
    axs[0, 1].plot(uu, soft_behavior.gradient_energy(uu + 1.0)[0], color=scm.to_rgba(b_i), lw=1,
             ls='--' if extrapolation else '-')
for _soft_behavior in _soft_behaviors:
    axs[0, 1].plot(uu, _soft_behavior.gradient_energy(uu + 1.0)[0], color='k', lw=0.4)
axs[0, 1].set_ylabel('force (N)')
axs[0, 1].set_xlabel('displacement (mm)')
axs[0, 1].set_xlim((-1.0, 20))
axs[0, 1].set_ylim((-0.1, 1.5))
cbar = plt.colorbar(scm, cax=None, ax=axs[0, 1])
cbar.ax.set_title('b (mm)')
cbar.ax.set_ylim((0.5, 2.5))

# stiffening f-d curves as function of c (panel C)
block_folders = ['stiff_0', 'stiff_1', 'stiff_2', 'stiff_3']
behaviors_paths = [os.path.join(block_library_dir, block_folder, 'bezier_behavior.csv')
                   for block_folder in block_folders]
c = [12, 13, 14, 15]
_stiff_behaviors = [read_behavior(behavior_path) for behavior_path in behaviors_paths]
control_point_sets = [stiff_behavior.get_control_points() for stiff_behavior in _stiff_behaviors]
u1 = interp1d(c, [cp[0][1] for cp in control_point_sets], fill_value='extrapolate')
u2 = interp1d(c, [cp[0][2] for cp in control_point_sets], fill_value='extrapolate')
u3 = interp1d(c, [cp[0][3] for cp in control_point_sets], fill_value='extrapolate')
f1 = interp1d(c, [cp[1][1] for cp in control_point_sets], fill_value='extrapolate')
f2 = interp1d(c, [cp[1][2] for cp in control_point_sets], fill_value='extrapolate')
f3 = interp1d(c, [cp[1][3] for cp in control_point_sets], fill_value='extrapolate')
uu = np.linspace(0, 25, 1000)
c_vector = np.arange(11, 16.5, .5)
cmap = 'GnBu'
cn = colors.Normalize(vmin=8.5, vmax=17.25)
scm = mcm.ScalarMappable(norm=cn, cmap=cmap)
for c_i in c_vector:
    extrapolation = c_i < 12 or c_i > 15
    stiff_behavior = BezierBehavior(1.0, u_i=[u1(c_i), u2(c_i), u3(c_i)], f_i=[f1(c_i), f2(c_i), f3(c_i)])
    axs[0, 2].plot(uu, stiff_behavior.gradient_energy(uu + 1.0)[0], color=scm.to_rgba(c_i), lw=1,
             ls='--' if extrapolation else '-')
for _stiff_behavior in _stiff_behaviors:
    axs[0, 2].plot(uu, _stiff_behavior.gradient_energy(uu + 1.0)[0], color='k', lw=0.4)
axs[0, 2].set_ylabel('force (N)')
axs[0, 2].set_xlabel('displacement (mm)')
axs[0, 2].set_xlim((-1.0, 25))
axs[0, 2].set_ylim((-0.1, 1.5))
cbar = plt.colorbar(scm, cax=None, ax=axs[0, 2])
cbar.ax.set_title('c (mm)')
cbar.ax.set_ylim(11, 16)

# nonmonotonic f-d curves as function of d (panel D)
block_folders = ['snap_3', 'snap_4', 'snap_5']
behaviors_paths = [os.path.join(block_library_dir, block_folder, 'bezier_behavior.csv')
                   for block_folder in block_folders]
d = [5, 7, 9]
_snap_behaviors = [read_behavior(behavior_path) for behavior_path in behaviors_paths]
control_point_sets = [soft_behavior.get_control_points() for soft_behavior in _snap_behaviors]
u1 = interp1d(d, [cp[0][1] for cp in control_point_sets], fill_value='extrapolate')
u2 = interp1d(d, [cp[0][2] for cp in control_point_sets], fill_value='extrapolate')
u3 = interp1d(d, [cp[0][3] for cp in control_point_sets], fill_value='extrapolate')
f1 = interp1d(d, [cp[1][1] for cp in control_point_sets], fill_value='extrapolate')
f2 = interp1d(d, [cp[1][2] for cp in control_point_sets], fill_value='extrapolate')
f3 = interp1d(d, [cp[1][3] for cp in control_point_sets], fill_value='extrapolate')
uu = np.linspace(0, 25, 1000)
d_vector = np.arange(4, 10.5, 0.5)
cmap = 'GnBu'
cn = colors.Normalize(vmin=1, vmax=11.5)  # 50% / 25%
scm = mcm.ScalarMappable(norm=cn, cmap=cmap)

for d_i in d_vector:
    extrapolation = d_i > 9 or d_i < 5
    snap_behavior = BezierBehavior(1.0, u_i=[u1(d_i), u2(d_i), u3(d_i)], f_i=[f1(d_i), f2(d_i), f3(d_i)])
    axs[1, 0].plot(uu, snap_behavior.gradient_energy(uu + 1.0)[0], color=scm.to_rgba(d_i), lw=1,
             ls='--' if extrapolation else '-')
for _snap_behavior in _snap_behaviors:
    axs[1, 0].plot(uu, _snap_behavior.gradient_energy(uu + 1.0)[0], color='k', lw=0.4)

axs[1, 0].set_ylabel('force (N)')
axs[1, 0].set_xlabel('displacement (mm)')
axs[1, 0].set_xlim((-1.0, 20))
axs[1, 0].set_ylim((-0.1, 1.5))
cbar = plt.colorbar(scm, cax=None, ax=axs[1, 0])
cbar.ax.set_title('d (mm)')
cbar.ax.set_ylim((4, 10))

# nonmonotonic f-d curves as a function of theta
block_folders = ['snap_2', 'snap_5', 'snap_8']
behaviors_paths = [os.path.join(block_library_dir, block_folder, 'bezier_behavior.csv')
                   for block_folder in block_folders]
theta = [70, 60, 50]
_snap_behaviors = [read_behavior(behavior_path) for behavior_path in behaviors_paths]
control_point_sets = [soft_behavior.get_control_points() for soft_behavior in _snap_behaviors]
u1 = interp1d(theta, [cp[0][1] for cp in control_point_sets], fill_value='extrapolate')
u2 = interp1d(theta, [cp[0][2] for cp in control_point_sets], fill_value='extrapolate')
u3 = interp1d(theta, [cp[0][3] for cp in control_point_sets], fill_value='extrapolate')
f1 = interp1d(theta, [cp[1][1] for cp in control_point_sets], fill_value='extrapolate')
f2 = interp1d(theta, [cp[1][2] for cp in control_point_sets], fill_value='extrapolate')
f3 = interp1d(theta, [cp[1][3] for cp in control_point_sets], fill_value='extrapolate')
uu = np.linspace(0, 25, 1000)
theta_vector = np.arange(80, 35, -5)
cmap = 'GnBu'
cn = colors.Normalize(vmin=20, vmax=90)  # 50% / 25%
scm = mcm.ScalarMappable(norm=cn, cmap=cmap)

for theta_i in theta_vector:
    extrapolation = theta_i < 50 or theta_i > 70
    snap_behavior = BezierBehavior(1.0, u_i=[u1(theta_i), u2(theta_i), u3(theta_i)],
                                   f_i=[f1(theta_i), f2(theta_i), f3(theta_i)])
    axs[1, 2].plot(uu, snap_behavior.gradient_energy(uu + 1.0)[0], color=scm.to_rgba(theta_i), lw=1,
             ls='--' if extrapolation else '-')
for _snap_behavior in _snap_behaviors:
    axs[1, 2].plot(uu, _snap_behavior.gradient_energy(uu + 1.0)[0], color='k', lw=0.4)

axs[1, 2].set_ylabel('force (N)')
axs[1, 2].set_xlabel('displacement (mm)')
axs[1, 2].set_xlim((-1.0, 20))
axs[1, 2].set_ylim((-0.1, 1.5))
cbar = plt.colorbar(scm, cax=None, ax=axs[1, 2])
cbar.ax.set_title('theta (deg)')
cbar.ax.set_ylim((40, 80))
plt.show()

'''
COUNTERSNAPPING MAGNITUDE AS A FUNCTION OF THE DESIGN PARAMETERS (FIGURE S5, PANELS F,G,H,I,J)
'''
data_dir = os.path.join('FIGURE S5', 'sensitivity_data')
_, axs = plt.subplots(2, 3)

# wrt to a (panel F)
a_vector = np.loadtxt(os.path.join(data_dir, 'a.csv'), delimiter=',')
snap_magnitudes = np.loadtxt(os.path.join(data_dir, 'snap_magnitude_data_vs_a.csv'), delimiter=',')
low_extrapolation_indices = a_vector < 4
high_extrapolation_indices = a_vector > 6
interpolation_indices = np.logical_and(np.logical_not(low_extrapolation_indices),
                                       np.logical_not(high_extrapolation_indices))
axs[0, 0].plot(a_vector[interpolation_indices], snap_magnitudes[interpolation_indices], 'k-', lw=.75)
axs[0, 0].plot(a_vector[low_extrapolation_indices], snap_magnitudes[low_extrapolation_indices], 'k--', lw=.75)
axs[0, 0].plot(a_vector[high_extrapolation_indices], snap_magnitudes[high_extrapolation_indices], 'k--', lw=.75)
current_index = np.argmin(np.abs(a_vector - 4))
axs[0, 0].plot(a_vector[current_index], snap_magnitudes[current_index], 'r^', markersize=2)
axs[0, 0].set_ylabel('force jump (N)')
axs[0, 0].set_xlabel('a (mm)')
axs[0, 0].set_xlim((np.min(a_vector), np.max(a_vector)))
axs[0, 0].set_ylim((-0.5, 0.3))

# wrt to b (panel G)
b_vector = np.loadtxt(os.path.join(data_dir, 'b.csv'), delimiter=',')
snap_magnitudes = np.loadtxt(os.path.join(data_dir, 'snap_magnitude_data_vs_b.csv'), delimiter=',')
low_extrapolation_indices = b_vector < 1
high_extrapolation_indices = b_vector > 2
interpolation_indices = np.logical_and(np.logical_not(low_extrapolation_indices),
                                       np.logical_not(high_extrapolation_indices))
axs[0, 1].plot(b_vector[interpolation_indices], snap_magnitudes[interpolation_indices], 'k-', lw=.75)
axs[0, 1].plot(b_vector[low_extrapolation_indices], snap_magnitudes[low_extrapolation_indices], 'k--', lw=.75)
axs[0, 1].plot(b_vector[high_extrapolation_indices], snap_magnitudes[high_extrapolation_indices], 'k--', lw=.75)
current_index = np.argmin(np.abs(b_vector - 1))
axs[0, 1].plot(b_vector[current_index], snap_magnitudes[current_index], 'r^', markersize=2)
axs[0, 1].set_ylabel('force jump (N)')
axs[0, 1].set_xlabel('b (mm)')
axs[0, 1].set_xlim((np.min(b_vector), np.max(b_vector)))
axs[0, 1].set_ylim((-0.5, 0.3))

# wrt to c (panel H)
c_vector = np.loadtxt(os.path.join(data_dir, 'c.csv'), delimiter=',')
snap_magnitudes = np.loadtxt(os.path.join(data_dir, 'snap_magnitude_data_vs_c.csv'), delimiter=',')
low_extrapolation_indices = c_vector < 12
high_extrapolation_indices = c_vector > 15
interpolation_indices = np.logical_and(np.logical_not(low_extrapolation_indices),
                                       np.logical_not(high_extrapolation_indices))
axs[0, 2].plot(c_vector[interpolation_indices], snap_magnitudes[interpolation_indices], 'k-', lw=.75)
axs[0, 2].plot(c_vector[low_extrapolation_indices], snap_magnitudes[low_extrapolation_indices], 'k--', lw=.75)
axs[0, 2].plot(c_vector[high_extrapolation_indices], snap_magnitudes[high_extrapolation_indices], 'k--', lw=.75)
current_index = np.argmin(np.abs(c_vector - 14))
axs[0, 2].plot(c_vector[current_index], snap_magnitudes[current_index], 'r^', markersize=2)
axs[0, 2].set_ylabel('force jump (N)')
axs[0, 2].set_xlabel('c (mm)')
axs[0, 2].set_xlim((11, 16))
axs[0, 2].set_ylim((-0.5, 0.3))

# wrt to d (panel I)
d_vector = np.loadtxt(os.path.join(data_dir, 'd.csv'), delimiter=',')
snap_magnitudes = np.loadtxt(os.path.join(data_dir, 'snap_magnitude_data_vs_d.csv'), delimiter=',')
low_extrapolation_indices = d_vector < 5
high_extrapolation_indices = d_vector > 9
interpolation_indices = np.logical_and(np.logical_not(low_extrapolation_indices),
                                       np.logical_not(high_extrapolation_indices))
axs[1, 0].plot(d_vector[interpolation_indices], snap_magnitudes[interpolation_indices], 'k-', lw=.75)
axs[1, 0].plot(d_vector[low_extrapolation_indices], snap_magnitudes[low_extrapolation_indices], 'k--', lw=.75)
axs[1, 0].plot(d_vector[high_extrapolation_indices], snap_magnitudes[high_extrapolation_indices], 'k--', lw=.75)
current_index = np.argmin(np.abs(d_vector - 9))
axs[1, 0].plot(d_vector[current_index], snap_magnitudes[current_index], 'r^', markersize=2)
axs[1, 0].set_ylabel('force jump (N)')
axs[1, 0].set_xlabel('d (mm)')
axs[1, 0].set_xlim((np.min(d_vector), np.max(d_vector)))
axs[1, 0].set_ylim((-0.5, 0.3))

# wrt to theta (panel J)
theta_vector = np.loadtxt(os.path.join(data_dir, 'theta.csv'), delimiter=',')
snap_magnitudes = np.loadtxt(os.path.join(data_dir, 'snap_magnitude_data_vs_theta.csv'), delimiter=',')
low_extrapolation_indices = theta_vector < 50
high_extrapolation_indices = theta_vector > 70
interpolation_indices = np.logical_and(np.logical_not(low_extrapolation_indices),
                                       np.logical_not(high_extrapolation_indices))
axs[1, 2].plot(theta_vector[interpolation_indices], snap_magnitudes[interpolation_indices], 'k-', lw=.75)
axs[1, 2].plot(theta_vector[low_extrapolation_indices], snap_magnitudes[low_extrapolation_indices], 'k--', lw=.75)
axs[1, 2].plot(theta_vector[high_extrapolation_indices], snap_magnitudes[high_extrapolation_indices], 'k--', lw=.75)
current_index = np.argmin(np.abs(theta_vector - 60))
axs[1, 2].plot(theta_vector[current_index], snap_magnitudes[current_index], 'r^', markersize=2)
axs[1, 2].set_ylabel('force jump (N)')
axs[1, 2].set_xlabel('theta (deg)')
axs[1, 2].set_xlim((np.min(theta_vector), np.max(theta_vector)))
axs[1, 2].set_ylim((-0.5, 0.3))
plt.show()

