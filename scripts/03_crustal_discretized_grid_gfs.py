import pickle as pkl
import numpy as np
import cutde.halfspace as HS

### USER INPUTS ###
crustal_mesh_version = "CFM"
out_files_directory = "mesh_gf_outfiles"
steeper_dip, gentler_dip = False, False
# smaller cells make nicer maps but take forever. ~4km grid ok for figures
cell_size = 10000  # in meters,
x, y = 1760934, 5431096  # central location of grid; Seaview
buffer = 12.e4  # in meters (area around Wellington to calculate displacements)

#########################
gf_type = "grid"
out_files_path = f"../{out_files_directory}/crustal_{crustal_mesh_version}"

# load files
with open(f"{out_files_path}/crustal_discretized_dict.pkl", "rb") as f:
    discretized_dict = pkl.load(f)

# Grid of x and y to calculate sea surface displacements at
x_data = np.arange(round(x-buffer, -3), round(x+buffer, -3), cell_size)
y_data = np.arange(round(y-buffer, -3), round(y+buffer, -3), cell_size)

xmesh, ymesh = np.meshgrid(x_data, y_data)
# xmesh is organized by grid rows (one item is one row). Each item is a duplicate of the x_data array.
# ymesh is organized by rows too (one item is one row). Each value in the item is the same value (the y_data value).
points_x = xmesh.flatten()
points_y = ymesh.flatten()
points_xyz = np.vstack((points_x, points_y, points_x * 0.)).T

# this just numbers the grid points to be consistent with the other named site files. At the moment it's not used.
# Later on it becomes important to keep track of how the grid is reshaped into a list to make sure the point number
# matches the grid cell.
site_name_list = list(range(len(x_data) * len(y_data)))

#make empty greens funcitons dictionary
gf_dict = {}

# calculates greens for 1 m slip on each fault section.
counter = 0
for fault_id in discretized_dict.keys():
    counter += 1
    triangles = discretized_dict[fault_id]["triangles"]
    rake = discretized_dict[fault_id]["rake"]

    zero_slip_array = np.zeros((triangles.shape[0],))
    ones_slip_array = np.ones((triangles.shape[0],))

    dip_slip_array = np.vstack([zero_slip_array, ones_slip_array, zero_slip_array]).T
    strike_slip_array = np.vstack([ones_slip_array, zero_slip_array, zero_slip_array]).T

    disps_ss = HS.disp_free(obs_pts=points_xyz, tris=triangles, slips=strike_slip_array, nu=0.25)
    disps_ds = HS.disp_free(obs_pts=points_xyz, tris=triangles, slips=dip_slip_array, nu=0.25)


    disp_dict = {"ss": disps_ss[:, -1], "ds": disps_ds[:, -1], "rake": rake, "site_coords": points_xyz,
                 "site_name_list": site_name_list, "x_data": x_data, "y_data": y_data}

    gf_dict[fault_id] = disp_dict
    if counter % 10 == 0:
        print(f"discretized dict fault {counter} of {len(discretized_dict.keys())}")


with open(f"{out_files_path}/crustal_gf_dict_{gf_type}.pkl", "wb") as f:
    pkl.dump(gf_dict, f)
