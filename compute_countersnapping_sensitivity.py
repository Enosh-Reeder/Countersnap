"""
By running this script (with Python library springable version 0.1.1 installed), the data used to analyze the
sensitivity of the countersnapping force jump with respect to the design parameters is computed via simulation.
Results have already been computed and stored in the folder FIGURE S5/sensitivity_data.
"""

from springable.simulation import solve_model
from springable.readwrite.fileio import read_model, read_parameters_from_model_file, read_behavior
from springable.graphics.plot import extract_loading_path
from springable.mechanics.static_solver import Result
import numpy as np
from scipy.interpolate import interp1d
import os

compute_sensitivity_wrt_a = True
compute_sensitivity_wrt_b = False
compute_sensitivity_wrt_c = True
compute_sensitivity_wrt_d = True
compute_sensitivity_wrt_theta = True

save_dir = os.path.join('FIGURE S5', 'sensitivity_data')
model_path = os.path.join('FIGURE S5', 'parametric_spring_network.csv')
block_library_dir = os.path.join('FIGURE S2', 'block_library')
n = 200  # sampling


def compute_snapping_magnitude(result: Result):
    (loading_path_indices,
     loading_critical_indices,
     loading_restabilization_indices) = extract_loading_path(result, drive_mode='displacement')
    if loading_critical_indices:
        u, f = result.get_equilibrium_path()
        return f[loading_restabilization_indices[0]] - f[loading_critical_indices[0]]
    else:
        return np.nan

if compute_sensitivity_wrt_a:
    # SENSITIVITY W.R.T PARAMETER A
    filename = 'snap_magnitude_data_vs_a.csv'
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
    a_vector = np.linspace(3, 7, n)
    initial_parameters, _ = read_parameters_from_model_file(model_path=model_path)
    snapping_magnitude = np.empty_like(a_vector)
    i = 0
    for a_i in a_vector:
        print(f'progress {i + 1}/{n}')
        parameters = initial_parameters.copy()
        parameters['u1so'] = u1(a_i)
        parameters['u2so'] = u2(a_i)
        parameters['u3so'] = u3(a_i)
        parameters['f1so'] = f1(a_i)
        parameters['f2so'] = f2(a_i)
        parameters['f3so'] = f3(a_i)
        mdl = read_model(model_path, parameters)
        res = solve_model(mdl, solver_settings=os.path.join('FIGURE S5', 'solver_settings.toml'))
        snapping_magnitude[i] = compute_snapping_magnitude(res)
        i += 1
    np.savetxt(os.path.join(save_dir, filename), snapping_magnitude, delimiter=',')
    np.savetxt(os.path.join(save_dir, 'a.csv'), a_vector, delimiter=',')

if compute_sensitivity_wrt_b:
    filename = 'snap_magnitude_data_vs_b.csv'
    block_folders = ['soft_0', 'soft_3', 'soft_6']
    behaviors_paths = [os.path.join(block_library_dir, block_folder, 'bezier_behavior.csv')
                       for block_folder in block_folders]
    b = [1, 1.5, 2]

    _soft_behaviors = [read_behavior(behavior_path) for behavior_path in behaviors_paths]

    control_point_sets = [soft_behavior.get_control_points() for soft_behavior in _soft_behaviors]

    u1 = interp1d(b, [cp[0][1] for cp in control_point_sets], fill_value='extrapolate')
    u2 = interp1d(b, [cp[0][2] for cp in control_point_sets], fill_value='extrapolate')
    u3 = interp1d(b, [cp[0][3] for cp in control_point_sets], fill_value='extrapolate')
    f1 = interp1d(b, [cp[1][1] for cp in control_point_sets], fill_value='extrapolate')
    f2 = interp1d(b, [cp[1][2] for cp in control_point_sets], fill_value='extrapolate')
    f3 = interp1d(b, [cp[1][3] for cp in control_point_sets], fill_value='extrapolate')
    b_vector = np.linspace(0.5, 2.5, n)
    initial_parameters, _ = read_parameters_from_model_file(model_path=model_path)
    snapping_magnitude = np.empty_like(b_vector)
    i = 0
    for b_i in b_vector:
        print(f'progress {i + 1}/{n}')
        parameters = initial_parameters.copy()
        parameters['u1so'] = u1(b_i)
        parameters['u2so'] = u2(b_i)
        parameters['u3so'] = u3(b_i)
        parameters['f1so'] = f1(b_i)
        parameters['f2so'] = f2(b_i)
        parameters['f3so'] = f3(b_i)
        mdl = read_model(model_path, parameters)
        res = solve_model(mdl, solver_settings=os.path.join('FIGURE S5', 'solver_settings.toml'))
        snapping_magnitude[i] = compute_snapping_magnitude(res)
        i += 1
    np.savetxt(os.path.join(save_dir, filename), snapping_magnitude, delimiter=',')
    np.savetxt(os.path.join(save_dir, 'b.csv'), b_vector, delimiter=',')

if compute_sensitivity_wrt_c:
    # SENSITIVITY W.R.T. PARAMETER C
    filename = 'snap_magnitude_data_vs_c.csv'
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
    c_vector = np.linspace(11.5, 18.5, n)
    initial_parameters, _ = read_parameters_from_model_file(model_path=model_path)
    snapping_magnitude = np.empty_like(c_vector)
    i = 0
    for a_i in c_vector:
        print(f'progress {i + 1}/{n}')
        parameters = initial_parameters.copy()
        parameters['u1st'] = u1(a_i)
        parameters['u2st'] = u2(a_i)
        parameters['u3st'] = u3(a_i)
        parameters['f1st'] = f1(a_i)
        parameters['f2st'] = f2(a_i)
        parameters['f3st'] = f3(a_i)
        mdl = read_model(model_path, parameters)
        res = solve_model(mdl, solver_settings=os.path.join('FIGURE S5', 'solver_settings.toml'))
        snapping_magnitude[i] = compute_snapping_magnitude(res)
        i += 1
    np.savetxt(os.path.join(save_dir, filename), snapping_magnitude, delimiter=',')
    np.savetxt(os.path.join(save_dir, 'c.csv'), c_vector, delimiter=',')

if compute_sensitivity_wrt_d:
    # SENSITIVITY W.R.T. PARAMETER D
    filename = 'snap_magnitude_data_vs_d.csv'
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
    d_vector = np.linspace(4, 10, n)
    initial_parameters, _ = read_parameters_from_model_file(model_path=model_path)
    snapping_magnitude = np.empty_like(d_vector)
    i = 0
    for d_i in d_vector:
        print(f'progress {i + 1}/{n}')
        parameters = initial_parameters.copy()
        parameters['u1sn'] = u1(d_i)
        parameters['u2sn'] = u2(d_i)
        parameters['u3sn'] = u3(d_i)
        parameters['f1sn'] = f1(d_i)
        parameters['f2sn'] = f2(d_i)
        parameters['f3sn'] = f3(d_i)

        mdl = read_model(model_path, parameters)
        res = solve_model(mdl, solver_settings=os.path.join('FIGURE S5', 'solver_settings.toml'))
        snapping_magnitude[i] = compute_snapping_magnitude(res)
        i += 1
    np.savetxt(os.path.join(save_dir, filename), snapping_magnitude, delimiter=',')
    np.savetxt(os.path.join(save_dir, 'd.csv'), d_vector, delimiter=',')

if compute_sensitivity_wrt_theta:
    # SENSITIVITY W.R.T. PARAMETER THETA
    filename = 'snap_magnitude_data_vs_theta.csv'
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
    theta_vector = np.linspace(40, 80, n)
    initial_parameters, _ = read_parameters_from_model_file(model_path=model_path)
    snapping_magnitude = np.empty_like(theta_vector)
    i = 0
    for theta_i in theta_vector:
        print(f'progress {i + 1}/{n}')
        parameters = initial_parameters.copy()
        parameters['u1sn'] = u1(theta_i)
        parameters['u2sn'] = u2(theta_i)
        parameters['u3sn'] = u3(theta_i)
        parameters['f1sn'] = f1(theta_i)
        parameters['f2sn'] = f2(theta_i)
        parameters['f3sn'] = f3(theta_i)

        mdl = read_model(model_path, parameters)
        res = solve_model(mdl, solver_settings=os.path.join('FIGURE S5', 'solver_settings.toml'))
        snapping_magnitude[i] = compute_snapping_magnitude(res)
        i += 1
    np.savetxt(os.path.join(save_dir, filename), snapping_magnitude, delimiter=',')
    np.savetxt(os.path.join(save_dir, 'theta.csv'), theta_vector, delimiter=',')

