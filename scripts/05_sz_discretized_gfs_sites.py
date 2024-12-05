import pickle as pkl
import numpy as np
import cutde.halfspace as HS
from time import time
from shapely.geometry import MultiPoint
import os

"""Calculates greens functions along map at specified site coordinates
Outputs a dictionary of greens functions solutions"""

###### USER INPUTS ########
sz_mesh_version = "multi50"
out_files_directory = "mesh_gf_outfiles_EXAMPLE"
steeper_dip, gentler_dip = False, False

# in list form for one coord or list of lists for multiple (in NZTM)
site_1_coord = np.array([1749376, 5427530, 0])   # downtown Wellington
site_2_coord = np.array([1736455, 5427195, 0])   # South Coast, sea/shore electricity cable location
site_3_coord = np.array([1751064, 5423128, 0])   # Wellington airport
site_4_coord = np.array([1754357, 5445716, 0])   # Porirua CBD north of Ohariu fault
site_5_coord = np.array([1754557, 5445119, 0])   # Porirua CBD south of Ohariu fault
site_6_coord = np.array([1757199, 5434207, 0])   # Petone (Hutt Valley); people and office buildings
site_7_coord = np.array([1759240, 5432111, 0])   # Seaview (Hutt Valley); oil tankers loading/unloading
site_8_coord = np.array([1766726, 5471342, 0])   # Paraparaumu; west coast
site_9_coord = np.array([1758789, 5427418, 0])   # Eastbourne (eastern wellington harbour)
site_10_coord = np.array([1760183, 5410911, 0])   # Turakirae Head
site_11_coord = np.array([1779348, 5415831, 0])   # Lake Ferry (small settlement, flood infrustructure)
site_12_coord = np.array([1789451, 5391086, 0])      # Cape Palliser (marine terraces to compare)
site_13_coord = np.array([1848038, 5429751, 0])     # Flat Point (round out point spacing)

site_coords = np.vstack((site_1_coord, site_2_coord, site_3_coord, site_4_coord, site_5_coord, site_6_coord,
                   site_7_coord, site_8_coord, site_9_coord, site_10_coord, site_11_coord, site_12_coord, site_13_coord))
site_name_list = ["Wellington CBD", "South Coast", "Wellington Airport", "Porirua CBD north", "Porirua CBD south",
                  "Petone", "Seaview", "Paraparaumu", "Eastbourne", "Turakirae Head", "Lake Ferry", "Cape Palliser",
                  "Flat Point"]

#############################################
out_files_path = f"../{out_files_directory}/sz_{sz_mesh_version}"

x_data = site_coords[:, 0]
y_data = site_coords[:, 1]
gf_type = "sites"

if steeper_dip == True and gentler_dip == False:
    sz_mesh_version += "_steeperdip"
elif gentler_dip == True and steeper_dip == False:
    sz_mesh_version += "_gentlerdip"
elif gentler_dip == False and steeper_dip == False:
    sz_mesh_version += ""
else:
    print("Dip modifications are wrong. Only one statement can be True at once. Try again.")
    exit()

# Load files
with open(f"{out_files_path}/sz_discretized_dict.pkl",
          "rb") as f:
    discretized_dict = pkl.load(f)

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

    disps = HS.disp_free(obs_pts=site_coords, tris=triangles, slips=total_slip_array, nu=0.25)

    # Set rake to 90 so that in future functions total displacement is just equal to DS
    # therefore, the saved SS/DS/rake items in the dictionary are not "true" to the mesh, but they are
    # equivalent to prevent script changes for the displacament step
    disp_dict = {"ss": disps[:, -1] * 0, "ds": disps[:, -1], "rake": 90, "site_coords": site_coords,
                 "site_name_list": site_name_list, "x_data": x_data, "y_data": y_data}

    gf_dict[fault_id] = disp_dict
    if fault_id % 10 == 0:
        print(
            f'discretized dict {fault_id} of {len(discretized_dict.keys())} processing ({triangles.shape[0]} triangles per patch)')

    # if fault_id % 1 == 0:
    #     print(f'discretized dict {fault_id} of {len(discretised_dict.keys())} done in {time() - begin:.2f} seconds ({triangles.shape[0]} triangles per patch)')


with open(f"{out_files_path}/sz_gf_dict_{gf_type}.pkl", "wb") as f:
    pkl.dump(gf_dict, f)

