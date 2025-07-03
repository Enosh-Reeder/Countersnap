import matplotlib.pyplot as plt
import os
from utils import read_force_displacement_data, extract_loading_sequence, extract_unloading_sequence


'''
COMPARISON COUNTERSNAPPING WITH PURELY AND PURELY PARALLEL COUNTERPARTS (FIGURE S7, PANEL A)
'''
data_dir = 'FIGURE S7'
_, up, fp = read_force_displacement_data(os.path.join(data_dir, 'purely_parallel', 'Specimen_rawData_1.csv'))
_, us, fs = read_force_displacement_data(os.path.join(data_dir, 'purely_series', 'Specimen_RawData_1.csv'))
_, u, f = read_force_displacement_data(os.path.join('FIGURE 2', 'displacement_driven', 'Specimen_RawData_1.csv'))

start, _ = extract_loading_sequence(up, 1)
_, end = extract_unloading_sequence(up, 1)
up = up[start:end].copy()
fp = fp[start:end].copy()

start, _ = extract_loading_sequence(us, 1)
_, end = extract_unloading_sequence(us, 1)
us = us[start:end].copy()
fs = fs[start:end].copy()

start, _ = extract_loading_sequence(u, 1)
_, end = extract_unloading_sequence(u, 1)
u = u[start:end].copy()
f = f[start:end].copy()


# plot
_, ax = plt.subplots()
ax.plot(u, f, 'o', markersize=2, color='#000000', label='countersnapping structure')
ax.plot(up, fp, 'o', markersize=2, color='#fdb338', label='purely parallel')
ax.plot(us, fs, 'o', markersize=2, color='#025196', label='purely series')
ax.legend()
plt.show()
