import matplotlib.pyplot as plt
import os
import numpy as np

'''
MONOTONIC BLOCK STICK SLIP EXPERIMENT (data not shown in article)
'''
data_dir = os.path.join('FIGURE 3', 'monotonic')

frame_delay = 105  # side cam started before
top_px_per_mm = 260 / 50  # top camera px_per_mm
side_slider_px_per_mm = 14.5
fps = 50.0
radius_connector_ball = 5.0  # mm
distance_main_hole_ball_center = 2.55  # mm

loading_start_end_indices = np.loadtxt(os.path.join(data_dir, 'side_cam', 'loading_direction_change_indices.csv')).astype(int)
slider_small_holes_x = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_small_holes', 'x.csv'), delimiter=',') / side_slider_px_per_mm
slider_small_holes_y = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_small_holes', 'y.csv'), delimiter=',') / side_slider_px_per_mm
slider_main_hole_x = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_main_hole', 'x.csv'), delimiter=',') / side_slider_px_per_mm
slider_main_hole_y = np.loadtxt(os.path.join(data_dir, 'side_cam','tracking_results_slider_main_hole', 'y.csv'), delimiter=',') / side_slider_px_per_mm
paper_top_x = np.loadtxt(os.path.join(data_dir, 'top_cam','tracking_results_top_paper_markers', 'x.csv'), delimiter=',') / top_px_per_mm
clamped_node_x = np.loadtxt(os.path.join(data_dir, 'top_cam', 'tracking_results_top_clamp_node', 'x.csv'), delimiter=',') / top_px_per_mm

# trigonometry to calculate the contact point position based on the connector holes positions
# and correct for rolling

# angle of the connector with ground
theta = np.arcsin((slider_main_hole_x - slider_small_holes_x[:, 0]) / np.sqrt(
    (slider_main_hole_x - slider_small_holes_x[:, 0]) ** 2 + (
            slider_main_hole_y - slider_small_holes_y[:, 0]) ** 2))
contact_point_x = slider_main_hole_x - distance_main_hole_ball_center * np.sin(theta) - radius_connector_ball * theta
contact_point_displacement_x = contact_point_x - contact_point_x[0]

side_view_nb_frames = slider_main_hole_x.shape[0]
top_view_nb_frames = paper_top_x.shape[0]
paper_top_position_x = np.mean(paper_top_x, axis=1) - clamped_node_x
paper_top_displacement_x = paper_top_position_x - paper_top_position_x[0]
nb_useful_frames = min(side_view_nb_frames - frame_delay, top_view_nb_frames)
t_combined = np.arange(nb_useful_frames) / fps
drift_combined = contact_point_displacement_x[
                 frame_delay:frame_delay + nb_useful_frames] + paper_top_displacement_x[:nb_useful_frames]

# plot
_, ax = plt.subplots()
ax.scatter(drift_combined, t_combined, c='#000000', s=2, zorder=4)
ax.set_ylabel('time (s)')
ax.set_xlabel('drift (mm)')

loading_start_end_indices -= frame_delay
for i in range(loading_start_end_indices.shape[0]):
    if i % 2 == 0:
        if i >= loading_start_end_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = loading_start_end_indices[i + 1]
        ax.axhspan(t_combined[loading_start_end_indices[i]], t_combined[end_index],
                        facecolor='#dddddd', label='loading' if i == 0 else '')
    else:
        if i >= loading_start_end_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = loading_start_end_indices[i + 1]
        ax.axhspan(t_combined[loading_start_end_indices[i]], t_combined[end_index],
                        facecolor='#eeeeee', label='unloading' if i == 1 else '')
ax.set_xlim((-25, 25))
ax.invert_yaxis()
ax.invert_xaxis()
ax.xaxis.set_label_position('top')
ax.legend(numpoints=3, loc='upper left')
plt.show()

'''
REGULAR SNAPPING STICK SLIP EXPERIMENT (FIGURE 3, PANEL C)
'''
data_dir = os.path.join('FIGURE 3', 'regular_snapping')

frame_delay = 102  # side cam started before
top_px_per_mm = 260 / 50
side_slider_px_per_mm = 14.5
fps = 50.0
radius_connector_ball = 5.0  # mm
distance_main_hole_ball_center = 2.55  # mm

loading_start_end_indices = np.loadtxt(os.path.join(data_dir, 'side_cam', 'loading_direction_change_indices.csv')).astype(int)
slider_small_holes_x = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_small_holes', 'x.csv'), delimiter=',') / side_slider_px_per_mm
slider_small_holes_y = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_small_holes', 'y.csv'), delimiter=',') / side_slider_px_per_mm
slider_main_hole_x = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_main_hole', 'x.csv'), delimiter=',') / side_slider_px_per_mm
slider_main_hole_y = np.loadtxt(os.path.join(data_dir, 'side_cam','tracking_results_slider_main_hole', 'y.csv'), delimiter=',') / side_slider_px_per_mm
paper_top_x = np.loadtxt(os.path.join(data_dir, 'top_cam','tracking_results_top_paper_markers', 'x.csv'), delimiter=',') / top_px_per_mm
clamped_node_x = np.loadtxt(os.path.join(data_dir, 'top_cam', 'tracking_results_top_clamp_node', 'x.csv'), delimiter=',') / top_px_per_mm
paper_top_x = np.loadtxt(os.path.join(data_dir, 'top_cam','tracking_results_top_paper_markers', 'x.csv'), delimiter=',') / top_px_per_mm
clamped_node_x = np.loadtxt(os.path.join(data_dir, 'top_cam', 'tracking_results_top_clamp_node', 'x.csv'), delimiter=',') / top_px_per_mm

# trigonometry to calculate the contact point position based on the connector holes positions
# and correct for rolling

# angle of the connector with ground
theta = np.arcsin((slider_main_hole_x - slider_small_holes_x[:, 0]) / np.sqrt(
    (slider_main_hole_x - slider_small_holes_x[:, 0]) ** 2 + (
            slider_main_hole_y - slider_small_holes_y[:, 0]) ** 2))
contact_point_x = slider_main_hole_x - distance_main_hole_ball_center * np.sin(theta) - radius_connector_ball * theta
contact_point_displacement_x = contact_point_x - contact_point_x[0]

side_view_nb_frames = slider_main_hole_x.shape[0]
top_view_nb_frames = paper_top_x.shape[0]
paper_top_position_x = np.mean(paper_top_x, axis=1) - clamped_node_x
paper_top_displacement_x = paper_top_position_x - paper_top_position_x[0]
nb_useful_frames = min(side_view_nb_frames - frame_delay, top_view_nb_frames)
t_combined = np.arange(nb_useful_frames) / fps
drift_combined = contact_point_displacement_x[
                 frame_delay:frame_delay + nb_useful_frames] + paper_top_displacement_x[:nb_useful_frames]

# plot
_, ax = plt.subplots()
ax.scatter(drift_combined, t_combined, c='#000000', s=2, zorder=4)
ax.set_ylabel('time (s)')
ax.set_xlabel('drift (mm)')

loading_start_end_indices -= frame_delay
for i in range(loading_start_end_indices.shape[0]):
    if i % 2 == 0:
        if i >= loading_start_end_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = loading_start_end_indices[i + 1]
        ax.axhspan(t_combined[loading_start_end_indices[i]], t_combined[end_index],
                        facecolor='#dddddd', label='loading' if i == 0 else '')
    else:
        if i >= loading_start_end_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = loading_start_end_indices[i + 1]
        ax.axhspan(t_combined[loading_start_end_indices[i]], t_combined[end_index],
                        facecolor='#eeeeee', label='unloading' if i == 1 else '')
ax.set_xlim((-25, 25))
ax.invert_yaxis()
ax.invert_xaxis()
ax.xaxis.set_label_position('top')
ax.legend(numpoints=3, loc='upper left')
plt.show()

'''
COUNTERSNAPPING STICK SLIP EXPERIMENT (FIGURE 3, PANEL D)
'''
data_dir = os.path.join('FIGURE 3', 'countersnapping')

frame_delay = 162  # side cam started before
top_px_per_mm = 260 / 50
side_slider_px_per_mm = 226 / 15
side_paper_px_per_mm = 789 / 50
fps = 50.0

radius_connector_ball = 5.0  # mm
distance_main_hole_ball_center = 2.55  # mm

loading_start_end_indices = np.loadtxt(os.path.join(data_dir, 'side_cam', 'loading_direction_change_indices.csv')).astype(int)
slider_small_holes_x = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_small_holes', 'x.csv'), delimiter=',') / side_slider_px_per_mm
slider_small_holes_y = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_small_holes', 'y.csv'), delimiter=',') / side_slider_px_per_mm
slider_main_hole_x = np.loadtxt(os.path.join(data_dir, 'side_cam', 'tracking_results_slider_main_hole', 'x.csv'), delimiter=',') / side_slider_px_per_mm
slider_main_hole_y = np.loadtxt(os.path.join(data_dir, 'side_cam','tracking_results_slider_main_hole', 'y.csv'), delimiter=',') / side_slider_px_per_mm
paper_top_x = np.loadtxt(os.path.join(data_dir, 'top_cam','tracking_results_top_paper_markers', 'x.csv'), delimiter=',') / top_px_per_mm
clamped_node_x = np.loadtxt(os.path.join(data_dir, 'top_cam', 'tracking_results_top_clamp_node', 'x.csv'), delimiter=',') / top_px_per_mm
paper_top_x = np.loadtxt(os.path.join(data_dir, 'top_cam','tracking_results_top_paper_markers', 'x.csv'), delimiter=',') / top_px_per_mm
clamped_node_x = np.loadtxt(os.path.join(data_dir, 'top_cam', 'tracking_results_top_clamp_node', 'x.csv'), delimiter=',') / top_px_per_mm

# trigonometry to calculate the contact point position based on the connector holes positions
# and correct for rolling

# angle of the connector with ground
theta = np.arcsin((slider_main_hole_x - slider_small_holes_x[:, 0]) / np.sqrt(
    (slider_main_hole_x - slider_small_holes_x[:, 0]) ** 2 + (
            slider_main_hole_y - slider_small_holes_y[:, 0]) ** 2))
contact_point_x = slider_main_hole_x - distance_main_hole_ball_center * np.sin(theta) - radius_connector_ball * theta
contact_point_displacement_x = contact_point_x - contact_point_x[0]

side_view_nb_frames = slider_main_hole_x.shape[0]
top_view_nb_frames = paper_top_x.shape[0]
paper_top_position_x = np.mean(paper_top_x, axis=1) - clamped_node_x
paper_top_displacement_x = paper_top_position_x - paper_top_position_x[0]
nb_useful_frames = min(side_view_nb_frames - frame_delay, top_view_nb_frames)
t_combined = np.arange(nb_useful_frames) / fps
drift_combined = contact_point_displacement_x[
                 frame_delay:frame_delay + nb_useful_frames] + paper_top_displacement_x[:nb_useful_frames]

# plot
_, ax = plt.subplots()
ax.scatter(drift_combined, t_combined, c='#000000', s=2, zorder=4)
ax.set_ylabel('time (s)')
ax.set_xlabel('drift (mm)')

loading_start_end_indices -= frame_delay
for i in range(loading_start_end_indices.shape[0]):
    if i % 2 == 0:
        if i >= loading_start_end_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = loading_start_end_indices[i + 1]
        ax.axhspan(t_combined[loading_start_end_indices[i]], t_combined[end_index],
                        facecolor='#dddddd', label='loading' if i == 0 else '')
    else:
        if i >= loading_start_end_indices.shape[0] - 2:
            end_index = -1
        else:
            end_index = loading_start_end_indices[i + 1]
        ax.axhspan(t_combined[loading_start_end_indices[i]], t_combined[end_index],
                        facecolor='#eeeeee', label='unloading' if i == 1 else '')
ax.set_xlim((-5, 45))
ax.invert_yaxis()
ax.invert_xaxis()
ax.xaxis.set_label_position('top')
ax.legend(numpoints=3, loc='upper left')
plt.show()