from springable.readwrite.fileio import read_results, read_model
from springable.graphics.plot import force_displacement_curve_in_ax
from springable.graphics.default_graphics_settings import PlotOptions
import os
import matplotlib.pyplot as plt
import numpy as np


'''
INDIVIDUAL NONLINEAR SPRING BEHAVIORS LEADING TO COUNTERSNAPPING UPON UNLOADING (FIGURE S4, PANEL B)
'''
data_dir = os.path.join('FIGURE S4')

mdl = read_model(os.path.join(data_dir, 'counterunsnapping_spring_network.csv'))
soft_spring_behavior = mdl.get_assembly().get_elements()[0].get_behavior()
stiff_spring_behavior = mdl.get_assembly().get_elements()[1].get_behavior()
snap_spring_behavior = mdl.get_assembly().get_elements()[2].get_behavior()
u = np.linspace(0, 3.5, 500)

# plot
_, ax = plt.subplots()
soft_spring_natural_length = soft_spring_behavior.get_natural_measure()
ax.plot(u, soft_spring_behavior.gradient_energy(soft_spring_natural_length + u)[0], 'b-')
stiff_spring_natural_length = stiff_spring_behavior.get_natural_measure()
ax.plot(u, stiff_spring_behavior.gradient_energy(stiff_spring_natural_length + u)[0], 'g-')
snap_spring_natural_length = snap_spring_behavior.get_natural_measure()
ax.plot(u, snap_spring_behavior.gradient_energy(snap_spring_natural_length + u)[0], 'r-')
ax.set_xlabel('displacement (mm)')
ax.set_ylabel('force (N)')
plt.show()

'''
COUNTER-UNSNAPPING SPRING NETWORK FORCE-DISPLACEMENT CURVE (FIGURE S4, PANEL C)
'''
data_dir = os.path.join('FIGURE S4')
res = read_results(os.path.join(data_dir, 'sim_results'))

# plot
_, ax = plt.subplots()
plot_options = PlotOptions()
force_displacement_curve_in_ax(res, ax, plot_options=plot_options)
ax.set_xlabel('displacement (mm)')
ax.set_ylabel('force (N)')
ax.legend()
plt.show()






