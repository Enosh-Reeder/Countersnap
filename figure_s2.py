from springable.readwrite.fileio import read_behavior, read_design_parameters, read_results
from springable.mechanics.static_solver import Result
from springable.graphics.plot import extract_loading_path, force_displacement_curve_in_ax
from springable.graphics.default_graphics_settings import PlotOptions
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import os
import numpy as np
from utils import read_force_displacement_data, extract_loading_sequence

'''
FITTING BEZIER CURVES VERSUS EXPERIMENTAL DATA COMPARISON (FIGURE S2, PANELS A,B,C)
'''
data_dir = os.path.join('FIGURE S2', 'block_library')
_, axs = plt.subplots(1, 3)

for i, folder in enumerate(os.listdir(data_dir)):
    match folder[:-2]:
        case 'soft':
            ax_ix = 0
        case 'stiff':
            ax_ix = 1
        case 'snap':
            ax_ix = 2
        case _:
            continue
    exp_u = np.loadtxt(os.path.join(data_dir, folder, 'experimental_displacements.csv'), delimiter=',')
    exp_f = np.loadtxt(os.path.join(data_dir, folder, 'experimental_forces.csv'), delimiter=',')
    final_u = np.loadtxt(os.path.join(data_dir, folder, 'final_u.csv'), delimiter=',')
    u = np.linspace(0, final_u, 500)
    fitting_behavior = read_behavior(behavior_path=os.path.join(data_dir, folder, 'bezier_behavior.csv'))
    f = fitting_behavior.gradient_energy(fitting_behavior.get_natural_measure() + u)[0]
    axs[ax_ix].plot(exp_u, exp_f, 'o', markersize=2, color='#999999', label='experiment')
    axs[ax_ix].plot(u, f, '-', color='#000000', label='fitting Bezier curve')
axs[0].set_xlabel('displacement (mm)')
axs[1].set_xlabel('displacement (mm)')
axs[2].set_xlabel('displacement (mm)')
axs[0].set_ylabel('force (N)')
plt.show()


'''
RESULTS OF THE BRUTE-FORCE SEARCH FOR COUNTERSNAPPING VIA NUMERICAL SIMULATION (FIGURE S2, PANEL E)
'''
data_dir = os.path.join('FIGURE S2', 'scanning_results')

def compute_snapping_magnitude(result: Result, drive_mode):
    (loading_path_indices,
     loading_critical_indices,
     loading_restabilization_indices) = extract_loading_path(result, drive_mode=drive_mode)
    if loading_critical_indices:
        u, f = result.get_equilibrium_path()
        return f[loading_restabilization_indices[0]] - f[loading_critical_indices[0]]
    else:
        return np.nan

norm = mcolors.Normalize(vmin=-0.22, vmax=0.22, clip=True)
sm = cm.ScalarMappable(norm=norm, cmap='PRGn_r')
colors = [np.empty(shape=(9, 9, 4)) for _ in range(4)]
processed = 0
for folder in os.listdir(data_dir):
    if not folder.startswith('sim'):
        continue
    sim_dir = os.path.join(data_dir, folder)
    param = read_design_parameters(sim_dir)
    i, j, k = int(param['snap'][-1]), int(param['soft'][-1]), int(param['stiff'][-1])
    res = read_results(sim_dir)
    fd_magnitude = compute_snapping_magnitude(res, 'force')
    if np.isnan(fd_magnitude):
        c = (252/255, 184/255, 86/255, 1.0)
    else:
        ud_magnitude = compute_snapping_magnitude(res, 'displacement')
        if np.isnan(ud_magnitude):
            c = (249/255, 213/255, 118/255, 1.0)
        else:
            c = sm.to_rgba(ud_magnitude)
    colors[k][j,i, :] = c
    processed += 1
    print(f'# simulations post-processed: {processed}')


# plot
fig, axs = plt.subplots(1, 4)

for stiff_ix in range(4):
    axs[stiff_ix].imshow(colors[stiff_ix], origin='lower', extent=[-0.5, 8.5, -0.5, 8.5])
    axs[stiff_ix].set_xticks(np.arange(-0.5, 9), minor=True)
    axs[stiff_ix].set_yticks(np.arange(-0.5, 9), minor=True)
    axs[stiff_ix].set_xticks(np.arange(0, 9))
    axs[stiff_ix].set_yticks(np.arange(0, 9))
    axs[stiff_ix].set_xticklabels(np.arange(0, 9))
    axs[stiff_ix].set_yticklabels(np.arange(0, 9))
    axs[stiff_ix].grid(True, which='minor', color='black')
    axs[stiff_ix].tick_params(which='minor', size=0)
cbar_ax = fig.add_axes([0.92, 0.4, 0.02, 0.2])
cbar = fig.colorbar(sm, cax=cbar_ax, orientation='vertical')
cbar.set_label("force jump (N)")
plt.show()


'''
COMPARISON SIMULATION WITH EXPERIMENT (FIGURE S2, PANEL G)
'''
sim_data_dir = os.path.join('FIGURE S2', 'scanning_results')
for folder in os.listdir(sim_data_dir):
    if not folder.startswith('sim'):
        continue
    params = read_design_parameters(os.path.join(sim_data_dir, folder))
    if params['soft'] == 'soft_0' and params['stiff'] == 'stiff_2' and params['snap'] == 'snap_5':
        res = read_results(os.path.join(sim_data_dir, folder))
        break
else:
    print('cannot find the selected combination')
    exit()


exp_data_dir = os.path.join('FIGURE 2', 'displacement_driven')
_, exp_u, exp_f = read_force_displacement_data(os.path.join(exp_data_dir, 'Specimen_rawData_1.csv'))
start, end = extract_loading_sequence(exp_u, 1)

# plot
_, ax = plt.subplots()
ax.plot(exp_u, exp_f, 'o', markersize=2, color='#ce4499', zorder=0, label='experiment')
plot_options = PlotOptions()
extra_options = {'color_for_stable_points': '#000000',
                 'color_for_stabilizable_points': '#000000',
                 'color_for_unstable_points': '#999999'}
plot_options.update(**extra_options)
force_displacement_curve_in_ax(res, ax, plot_options=plot_options, label='simulation')
ax.set_xlabel('displacement (mm)')
ax.set_ylabel('force (N)')
ax.legend()
plt.show()











