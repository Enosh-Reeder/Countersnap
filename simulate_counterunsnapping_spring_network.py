"""
By running this script (with Python library springable version 0.1.1 installed),
the equilibrium path of the counter-unsnapping spring network is computed.
Results have already been computed and stored in the folder FIGURE S4/sim_results.
"""
from springable.simulation import simulate_model
simulate_model('FIGURE S4/counterunsnapping_spring_network.csv', 'FIGURE S4/sim_results')


