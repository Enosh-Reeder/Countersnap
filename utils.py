import csv
import numpy as np
from matplotlib.colors import to_rgba, to_hex


def rgba_to_hex(rgba, background=(1, 1, 1)):
    r_fg, g_fg, b_fg, alpha = rgba
    r_bg, g_bg, b_bg = background
    r = alpha * r_fg + (1 - alpha) * r_bg
    g = alpha * g_fg + (1 - alpha) * g_bg
    b = alpha * b_fg + (1 - alpha) * b_bg
    return to_hex((r, g, b))


def read_force_displacement_data(data_dir):
    exp_time = []
    exp_displacements = []
    exp_forces = []
    with open(data_dir, newline='') as file_object:
        reader = csv.reader(file_object, delimiter=';')
        reading_data = False
        for row in reader:
            if not reading_data:
                if row:
                    try:
                        x = float(row[0].replace(',', '.'))
                    except ValueError:
                        continue
                    else:
                        reading_data = True
            if reading_data:
                exp_time.append(float(row[0].replace(',', '.')))
                exp_displacements.append(float(row[1].replace(',', '.')))
                exp_forces.append(float(row[2].replace(',', '.')))
    return np.array(exp_time), np.array(exp_displacements), np.array(exp_forces)

def extract_loading_sequence(u, loading_sequence_index):
    loading_start_indices, loading_end_indices, _, _ = _get_sequence_indices_no_pause(u)
    start_index = loading_start_indices[loading_sequence_index]
    end_index = loading_end_indices[loading_sequence_index]
    return start_index, end_index


def extract_unloading_sequence(u, unloading_sequence_index):
    _, _, unloading_start_indices, unloading_end_indices = _get_sequence_indices_no_pause(u)
    start_index = unloading_start_indices[unloading_sequence_index]
    end_index = unloading_end_indices[unloading_sequence_index]
    return start_index, end_index


def _get_sequence_indices_no_pause(u):
    """ assumes no pause after unloading sequence """
    start_index = np.argmax(u > 1e-3)
    loading_start_indices = []
    loading_end_indices = []
    unloading_start_indices = []
    unloading_end_indices = []
    loading = True
    unloading = False
    latest_index = start_index
    loading_start_indices.append(start_index)
    while True:
        if loading:
            relative_index = np.argmax(np.diff(u[latest_index:]) < -1e-3)
            if relative_index == 0:
                loading_end_indices.append(u.shape[0] - 1)
                break
            loading_end_indices.append(latest_index + relative_index)
            unloading_start_indices.append(latest_index + relative_index)
            loading = False
            unloading = True
            latest_index += relative_index
            continue
        if unloading:
            relative_index = np.argmax(np.diff(u[latest_index:]) > -1e-9)
            if relative_index == 0:
                unloading_end_indices.append(u.shape[0] - 1)
                break
            unloading_end_indices.append(latest_index + relative_index)
            loading_start_indices.append(latest_index + relative_index)
            loading = True
            unloading = False
            latest_index += relative_index
            continue
    return (loading_start_indices, loading_end_indices,
            unloading_start_indices, unloading_end_indices)