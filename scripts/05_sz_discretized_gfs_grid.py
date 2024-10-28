import pickle as pkl
import numpy as np
import cutde.halfspace as HS
import os
from time import time

"""Calculates greens functions along map at grid points
Outputs a dictionary of greens functions solutions"""

##### USER INPUTS #####
sz_mesh_version = "_multi50_gentlerdip"     # must match mesh_gf_outfiles subdirectory name
out_files_directory = "mesh_gf_outfiles_r1_testing"
cell_size = 5000            # in meters; spacing between grid points. ~4k makes ok map for Wellington, but it's slow.
x, y = 1760934, 5431096     # central location of grid; Seaview
buffer_size = 12.e4         # in meters (area around Wellington to calculate displacements)

#######################
gf_type = "grid"

# load files: open discretized dict of subduction interface to calculate grenns functions over
out_files_path = f"../{out_files_directory}/sz{sz_mesh_version}"
with open(f"{out_files_path}/sz_discretized_dict.pkl", "rb") as f:
    discretized_dict = pkl.load(f)

# Grid of x and y to calculate sea surface displacements at
x_data = np.arange(round(x-buffer_size, -2), round(x+buffer_size, -2), cell_size)
y_data = np.arange(round(y-buffer_size, -2), round(y+buffer_size, -2), cell_size)

xmesh, ymesh = np.meshgrid(x_data, y_data)
points_x = xmesh.flatten()
points_y = ymesh.flatten()
points_xyz = np.vstack((points_x, points_y, points_x * 0.)).T


# this just numbers the grid points to be consistent with the other named site files. At the moment it's not used.
# Later on it becomes important to keep track of how the grid is reshaped into a list to make sure the point number
# matches the grid cell.
site_name_list = list(range(len(x_data) * len(y_data)))

gf_dict = {}
for fault_id in discretized_dict.keys():
    triangles = discretized_dict[fault_id]["triangles"]
    rake = discretized_dict[fault_id]["rake"]

    # Get DS and SS components for each triangle element, depending on the element rake
    ss_comp = np.cos(np.radians(rake))
    ds_comp = np.sin(np.radians(rake))
    total_slip_array = np.zeros([triangles.shape[0], 3])

    begin = time()
    # Calculate the slip components for each triangle element
    for tri in range(triangles.shape[0]):
        ss, ds = np.linalg.lstsq(np.array([ss_comp[tri], ds_comp[tri]]).reshape([1, 2]), np.array([1]).reshape([1, 1]),
                                 rcond=None)[0]
        total_slip_array[tri, :2] = np.array([ss[0], ds[0]])

    disps = HS.disp_free(obs_pts=points_xyz, tris=triangles, slips=total_slip_array, nu=0.25)

    # Set rake to 90 so that in future functions total displacement is just equal to DS
    # therefore, the saved SS/DS/rake items in the dictionary are not "true" to the mesh, but they are
    # equivalent to prevent script changes for the displacement step
    disp_dict = {"ss": disps[:, -1] * 0, "ds": disps[:, -1], "rake": 90, "site_coords": points_xyz,
                 "site_name_list": site_name_list, "x_data": x_data, "y_data": y_data}

    gf_dict[fault_id] = disp_dict
    if fault_id % 10 == 0:
        print(
            f'discretized dict {fault_id} of {len(discretized_dict.keys())} processing ({triangles.shape[0]} triangles per patch)')

    # if fault_id % 1 == 0:
    #     print(f'discretized dict {fault_id} of {len(discretised_dict.keys())} done in {time() - begin:.2f} seconds ({triangles.shape[0]} triangles per patch)')


# save greens function dictionary as pickle file to outfiles path
with open(f"{out_files_path}/sz_gf_dict_{gf_type}.pkl", "wb") as f:
    pkl.dump(gf_dict, f)

