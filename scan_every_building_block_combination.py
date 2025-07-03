"""
By running this script (with Python library springable version 0.1.1 installed),
the same scanning results shown in the article will be replicated.
Results have already been computed and stored in the folder FIGURE S2/scanning_results.
"""
from springable.simulation import scan_parameter_space
scan_parameter_space('FIGURE S2/spring_network.csv', 'FIGURE S2/scanning_results',
                     graphics_settings_path='FIGURE S2/graphics_settings.toml',
                     solver_settings_path='FIGURE S2/solver_settings.toml',
                     scan_parameters_one_by_one=False)