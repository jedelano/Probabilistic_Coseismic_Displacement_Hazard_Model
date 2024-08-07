import pickle as pkl
import numpy as np
import geopandas as gpd
import cutde.halfspace as HS
from shapely.geometry import MultiPoint


# Calculates greens functions along map at specified interval
# Read in the geojson file from the NSHM inversion solution
sz_mesh_version = "_multi50"
out_files_directory = "mesh_gf_outfiles"
steeper_dip, gentler_dip = False, False
point_dist = 1000   # in meters

########
gf_type = "coast"

if steeper_dip == True and gentler_dip == False:
    sz_mesh_version += "_steeperdip"
elif gentler_dip == True and steeper_dip == False:
    sz_mesh_version += "_gentlerdip"
elif gentler_dip == False and steeper_dip == False:
    sz_mesh_version += ""
else:
    print("Dip modifications are wrong. Only one statement can be True at once. Try again.")
    exit()

out_files_path = f"../{out_files_directory}/sz{sz_mesh_version}"

# Load files
target_coastline = gpd.read_file("../data/map/wellington_coast_merged.geojson")
# distance between points along map

with open(f"{out_files_path}/sz_discretized_dict.pkl",
          "rb") as f:
    discretized_dict = pkl.load(f)

# find points at point_dist intervals along each line
coast_points = []
for line in target_coastline.geometry:  # shapely multiline#
    # step through each line at increasing distances
    # starts at 0 km and goes up to largest whole km
    for i in range(0, int(line.length), point_dist):
        coast_point_i = line.interpolate(i)
        coast_points.append(coast_point_i)

#### For filtering scenarios based on disp
# get x and y data for points along map
x_data, y_data = [], []
for point in coast_points:
    x_data.append(point.x)
    y_data.append(point.y)

points_x, points_y = np.array(x_data), np.array(y_data)
#np.save(f"{out_files_path}/xpoints{sz_mesh_version}.npy", xpoints)
#np.save(f"{out_files_path}/ypoints{sz_mesh_version}.npy", ypoints)
points_xyz = np.vstack((points_x, points_y, points_x * 0.)).T

# this just numbers the map points to be consistent with the other named site files.
site_name_list = range(len(points_x))

gf_dict_coast = {}
for fault_id in discretized_dict.keys():
    triangles = discretized_dict[fault_id]["triangles"]
    rake = discretized_dict[fault_id]["rake"]
    if len(triangles) > 0:
        vertices = triangles.reshape(triangles.shape[0] * triangles.shape[1], 3)
        vertex_multipoint = MultiPoint(vertices)

        zero_slip_array = np.zeros((triangles.shape[0],))
        ones_slip_array = np.ones((triangles.shape[0],))

        dip_slip_array = np.vstack([zero_slip_array, ones_slip_array, zero_slip_array]).T
        strike_slip_array = np.vstack([ones_slip_array, zero_slip_array, zero_slip_array]).T

        # calculate displacements
        disps_ss = HS.disp_free(obs_pts=points_xyz, tris=triangles, slips=strike_slip_array, nu=0.25)
        disps_ds = HS.disp_free(obs_pts=points_xyz, tris=triangles, slips=dip_slip_array, nu=0.25)

        # make displacement dictionary for outputs. only use vertical disps (last column)
        disp_dict = {"ss": disps_ss[:, -1], "ds": disps_ds[:, -1], "rake": rake, "site_name_list": site_name_list,
                     "site_coords": points_xyz}

        gf_dict_coast[fault_id] = disp_dict


with open(f"{out_files_path}/sz_gf_dict_{gf_type}.pkl", "wb") as f:
    pkl.dump(gf_dict_coast, f)

