import matplotlib.pyplot as plt
from matplotlib import cm
import os
from matplotlib.colors import to_rgba
from utils import read_force_displacement_data, extract_loading_sequence, rgba_to_hex


'''
FORCE-DISPLACEMENT CURVES OF THE BUILDING BLOCKS (FIGURE 1, PANELS J,K,L)
'''
data_dir = os.path.join('FIGURE 1', 'block_tensile_tests')
viridis = cm.get_cmap('viridis')

block_parameters = {
    'soft_0': {'a': 4, 'b': 1.0},
    'soft_1': {'a': 5, 'b': 1.0},
    'soft_2': {'a': 6, 'b': 1.0},
    'soft_3': {'a': 4, 'b': 1.5},
    'soft_4': {'a': 5, 'b': 1.5},
    'soft_5': {'a': 6, 'b': 1.5},
    'soft_6': {'a': 4, 'b': 2.0},
    'soft_7': {'a': 5, 'b': 2.0},
    'soft_8': {'a': 6, 'b': 2.0},
    'stiff_0': {'c': 12},
    'stiff_1': {'c': 13},
    'stiff_2': {'c': 14},
    'stiff_3': {'c': 15},
    'snap_0': {'d': 5, 'theta': 70},
    'snap_1': {'d': 7, 'theta': 70},
    'snap_2': {'d': 9, 'theta': 70},
    'snap_3': {'d': 5, 'theta': 60},
    'snap_4': {'d': 7, 'theta': 60},
    'snap_5': {'d': 9, 'theta': 60},
    'snap_6': {'d': 5, 'theta': 50},
    'snap_7': {'d': 7, 'theta': 50},
    'snap_8': {'d': 9, 'theta': 50},
}

_, axs = plt.subplots(3, 1)
for folder in os.listdir(data_dir):
    t, u, f = read_force_displacement_data(os.path.join(data_dir, folder, 'Specimen_RawData_1.csv'))
    start, end = extract_loading_sequence(u, 1)  # extract second loading sequence
    color = None
    match folder[:-2]:
        case 'soft':
            alpha = 1 - (block_parameters[folder]['a'] - 4) / 4
            rgb  = viridis( (block_parameters[folder]['b'] -1) / 1.25)
            color = rgba_to_hex(to_rgba(rgb, alpha))
            ax_ix = 0
        case 'stiff':
            rgb = viridis((block_parameters[folder]['c'] - 12) / 4)
            alpha = 1
            color = rgba_to_hex(to_rgba(rgb, alpha))
            ax_ix = 1
        case 'snap':
            alpha = 1 - (block_parameters[folder]['d'] - 5) / 8
            rgb = viridis((block_parameters[folder]['theta'] - 50) / 25)
            color = rgba_to_hex(to_rgba(rgb, alpha))
            ax_ix = 2
        case _:
            continue
    axs[ax_ix].plot(u[start:end], f[start:end], 'o', color=color, markersize=2)
    axs[ax_ix].set_ylabel('force (N)')
axs[2].set_xlabel('displacement (mm)')
plt.show()

