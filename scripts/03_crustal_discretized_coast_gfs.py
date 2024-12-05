import pickle as pkl
import numpy as np
import cutde.halfspace as HS
from shapely.geometry import MultiPoint
import geopandas as gpd

"""calculates green's functions at regularly-spaced points along a line vector files (geojson)"""
###### USER INPUTS ######
crustal_mesh_version = "CFM"
out_files_directory = "mesh_gf_outfiles_EXAMPLE"
steeper_dip, gentler_dip = False, False
# distance between points along map
point_dist = 500   # in meters

#########################
gf_type = "coast"
out_files_path = f"../{out_files_directory}/crustal_{crustal_mesh_version}"
# load files
with open(f"{out_files_path}/crustal_discretized_dict.pkl", "rb") as f:
    discretised_dict = pkl.load(f)

target_coastline = gpd.read_file("../data/map/wellington_coast_merged.geojson")

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
points_xyz = np.vstack((points_x, points_y, points_x * 0.)).T

# this just numbers the map points to be consistent with the other named site files.
site_name_list = range(len(points_x))

gf_dict_coast = {}
for i in discretised_dict.keys():
    triangles = discretised_dict[i]["triangles"]
    rake = discretised_dict[i]["rake"]
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

        # make displacement dictionary for outputs. only use the vertical disps. (last column)
        disp_dict = {"ss": disps_ss[:, -1], "ds": disps_ds[:, -1], "rake": rake, "site_name_list": site_name_list,
                     "site_coords": points_xyz}

        gf_dict_coast[i] = disp_dict

with open(f"{out_files_path}/crustal_gf_dict_{gf_type}.pkl", "wb") as f:
    pkl.dump(gf_dict_coast, f)


