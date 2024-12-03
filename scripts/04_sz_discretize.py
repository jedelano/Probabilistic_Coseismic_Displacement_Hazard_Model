import geopandas as gpd
import pandas as pd
import numpy as np
import meshio
from shapely.geometry import Polygon, Point
import pickle as pkl
import os
from pcdhm.shared_helper import check_triangle_normal

#### USER INPUT #####
# specify in the mesh version string if the dip angle has been modified e.g. "_multi50_steeperdip" or "_gentlerdip"
sz_mesh_version = "multi50_steeperdip"
out_files_directory = "mesh_gf_outfiles"

# Sensitivity testing for subduction interface depth
steeper_dip = True
gentler_dip = False

# this can be any working branch, should be the same for all.
NSHM_directory = "NZSHM22_AveragedInversionSolution-QXV0b21hdGlvblRhc2s6MTA3MzMx"

# De-blobify outputs
deblobify = True

#######################
if steeper_dip and gentler_dip:
    print("Dip modifications are wrong. Only one statement can be True at once. Try again.")
    exit()

# make outfiles directory if it doesn't already exist
if not os.path.exists(f"../{out_files_directory}"):
    os.mkdir(f"../{out_files_directory}")

out_files_path = f"../{out_files_directory}/sz_{sz_mesh_version}"
if not os.path.exists(out_files_path):
    os.mkdir(out_files_path)

# Read in the geojson file from the NSHM inversion solution
traces = gpd.read_file(f"../data/sz_solutions/{NSHM_directory}/ruptures/fault_sections.geojson").to_crs(epsg=2193)

# get rectangle centroid and polygon vertices based on NSHM fault data
all_rectangle_centroids = []
all_rectangle_polygons = []
# make dataframe forrectangle polygon attribtes (for export later)
df1 = pd.DataFrame()

df_rectangle_centroid = pd.DataFrame()

# Turn section traces into rectangular patches using the metadata in the GeoJSON file
for i, trace in traces.iterrows():
    # Convert the trace to a numpy array
    trace_array = np.array(trace.geometry.coords)
    # Convert the depths to metres
    trace_array[:, -1] *= -1000.

    initial_patch_height = (trace.LowDepth - trace.UpDepth) * 1000.
    initial_tandip = np.tan(np.radians(trace.DipDeg))
    patch_horizontal_dist = initial_patch_height / initial_tandip

    # uses geometry from NSHM
    low_depth, up_depth = trace.LowDepth, trace.UpDepth
    # Calculate the centroid of the trace
    trace_centroid = np.array([*trace.geometry.centroid.coords[0], trace.UpDepth * -1000])
    # Calculate the height of the patch
    patch_height = (trace.LowDepth - trace.UpDepth) * 1000.
    dip_angle = trace.DipDeg

    # write patch attributes to dictionary and add to bottom of data frame
    df2 = pd.DataFrame({'fault_id': [int(trace.FaultID)], 'dip_deg': [dip_angle], 'patch_height_m': [patch_height],
                        'up_depth_km': [up_depth], 'low_depth_km': [low_depth]}, index=[int(trace.FaultID)])
    df1 = pd.concat([df1, df2])

    df2_rectangle_centroid = pd.DataFrame({'fault_id': [trace.FaultID], 'depth': [np.mean([trace.UpDepth, trace.LowDepth])]}, index=[int(trace.FaultID)])
    df_rectangle_centroid = pd.concat([df_rectangle_centroid, df2_rectangle_centroid])
    #######################

    # Calculate the strike of the trace
    trace_strike = trace_array[1] - trace_array[0]
    # Normalise the strike vector
    trace_strike = trace_strike / np.linalg.norm(trace_strike)
    # Calculate the across-strike vector, by rotating the strike vector 90 degrees
    # aka dip direction vector (x, y), Normalized.
    across_strike = np.matmul(np.array([[0, 1], [-1, 0]]), trace_strike[:-1])

    # Calculate the down-dip vector by incorporating the dip angle
    cosdip = np.cos(np.radians(dip_angle))
    sindip = np.sin(np.radians(dip_angle))
    down_dip = np.array([cosdip * across_strike[0], cosdip * across_strike[1], - 1 * sindip])

    # Calculate the centroid of the patch
    rectangle_centroid = trace_centroid + (patch_height / sindip) / 2 * down_dip

    # Calculate the corners of the patch and make a shapely polygon
    dd1 = trace_array[1] + (patch_height / sindip) * down_dip
    dd2 = trace_array[0] + (patch_height / sindip) * down_dip
    rectangle_polygon = Polygon(np.vstack([trace_array, dd1, dd2]))

    # Append the patch centroid and polygon to lists
    all_rectangle_centroids.append(rectangle_centroid)
    all_rectangle_polygons.append(rectangle_polygon)

# write rectangle centroid and rectangle polygons to geojson
all_rectangle_centroids = np.array(all_rectangle_centroids)
all_rectangle_centroids_gs = gpd.GeoSeries([Point(centroid) for centroid in all_rectangle_centroids], crs=2193)
all_rectangle_centroids_gdf = gpd.GeoDataFrame(df_rectangle_centroid, geometry=all_rectangle_centroids_gs.geometry, crs=2193)
all_rectangle_centroids_gdf.to_file(
    f"{out_files_path}/all_rectangle_centroids.geojson", driver="GeoJSON")

all_rectangle_outline_gs = gpd.GeoSeries(all_rectangle_polygons, crs=2193)
all_rectangle_outline_gdf = gpd.GeoDataFrame(df1, geometry=all_rectangle_outline_gs.geometry, crs=2193)
all_rectangle_outline_gdf.to_file(f"{out_files_path}/all_rectangle_outlines.geojson", driver="GeoJSON", engine="fiona")

# read in triangle mesh and add the patch centroids as points
mesh = meshio.read(f"../data/hik_kerk3k_with_rake.vtk")
mesh_rake = mesh.cell_data["rake"][0]
mesh_triangles = mesh.cells_dict["triangle"]    # indices of vertices that make up triangles
mesh_vertices = mesh.points              # xyz of vertices
mesh_centroids = np.mean(mesh_vertices[mesh_triangles], axis=1)

mesh_vertices = mesh.points # xyz of vertices in mesh. indexed.
# array of 3 xyz arrays. (three sets of vertices to make a triangle)
triangle_vertex_arrays = mesh_vertices[mesh_triangles]

###### ensuring correct strike convention for accurate displacment calculations later
# calculate the normal vector to each triangle. If negative, reverse ordering of triangle of vertices.
ordered_triangle_vertex_arrays = []
for triangle in triangle_vertex_arrays:
    ordered_triangle_array = check_triangle_normal(triangle)
    ordered_triangle_vertex_arrays.append(ordered_triangle_array)
ordered_triangle_vertex_arrays = np.array(ordered_triangle_vertex_arrays)

triangle_centroids = np.mean(ordered_triangle_vertex_arrays, axis=1)   # three part arrray. mean of x, mean of y, mean of z.

# Find the rake for each rectangle patch by finding closest triangle centroid to the rectangle centroid
rectangle_rake = []
for rectangle_centroid in all_rectangle_centroids:
    distances = np.linalg.norm(rectangle_centroid - mesh_centroids, axis=1)
    if distances.min() < 2.e4:
        closest_triangle = np.argmin(distances)
        rectangle_rake.append(mesh_rake[closest_triangle])
    else:
        rectangle_rake.append(np.nan)
rectangle_rake = np.array(rectangle_rake)

# find the closest rectangle to each triangle centroid
closest_rectangles = []
for ix, triangle_centroid in enumerate(triangle_centroids):
    distances = np.linalg.norm(all_rectangle_centroids - triangle_centroid, axis=1)
    if distances.min() < 2.2e4:
        nearest = np.where(distances < 2.2e4)[0]
        vsep = all_rectangle_centroids[:, 2] - triangle_centroid[2]
        abs_sep = abs(vsep)
        if np.sum(distances < 2.2e4) == 1 or not deblobify:  # If only one option, use that option. Alternatively, use geographically nearest if using original, blobify method
            closest_rectangle = np.argmin(distances)
        else:
            nearest2 = nearest[np.argsort(abs_sep[nearest])[:2]]
            if np.diff(abs_sep[nearest2])[0] < 1.5e3:  # If nearest 2 are within 1.5 km vertical seperation, use the geographically nearest
                closest_rectangle = nearest2[np.argmin(distances[nearest2])]
            else:
                closest_rectangle = nearest2[np.argmin(abs_sep[nearest2])]  # If nearest 2 are > 1 km vertical seperation, use the structurally nearest
        closest_rectangles.append(closest_rectangle)
    else:
        closest_rectangles.append(-1)

# Manually correct some triangles
if os.path.exists('../data/mesh_corrections.csv'):
    print('Manually correcting some triangles')
    with open('../data/mesh_corrections.csv', 'r') as f:
        corrections = [[int(val) for val in line.strip().split(',')] for line in f.readlines()]

    for tri, closest_rectangle in corrections:
        closest_rectangles[tri] = closest_rectangle

# Prevent isolated triangles
if os.path.exists('../data/hik_kerk3k_neighbours.txt'):
    print('Removing isolated triangles')
    with open('../data/hik_kerk3k_neighbours.txt', 'r') as f:
        neighbours = [[int(tri) for tri in line.strip().split()] for line in f.readlines()]

    for tri in range(len(closest_rectangles)):
        rect = closest_rectangles[tri]
        if rect != -1:
            if len(neighbours[tri]) == 3:
                neigh = [closest_rectangles[neigh] for neigh in neighbours[tri] if closest_rectangles[neigh] != -1]
                if sum([ix != rect for ix in neigh]) >= 2 and len(neigh) == 3:
                    rects, count = np.unique(neigh, return_counts=True)
                    closest_rectangles[tri] = rects[np.argmax(count)]
elif deblobify:
    print('No neighbour file found - final output may include isolated triangles')

closest_rectangles = np.array(closest_rectangles)

# Create polygons from triangles
discretized_polygons = []
discretized_dict = {}
for index in traces.index:
    triangles_locs = np.where(closest_rectangles == index)[0]
    triangles = ordered_triangle_vertex_arrays[triangles_locs]
    # modify the depths to test dip changes
    if steeper_dip:
        triangles[:, :, 2] *= 1.15
    if gentler_dip:
        triangles[:, :, 2] *= 0.85

    # make dictionary of triangles that go with each polygon
    triangle_polygons = [Polygon(triangle) for triangle in triangles]
    dissolved_triangle_polygons = gpd.GeoSeries(triangle_polygons).union_all()
    discretized_polygons.append(dissolved_triangle_polygons)
    discretized_dict[index] = {"triangles": triangles, "rake": mesh_rake[triangles_locs]}

# Create a geodataframe and geospon file from the polygons
gdf = gpd.GeoDataFrame({"rake": rectangle_rake, "geometry": discretized_polygons, "fault_id": traces.index})
gdf.crs = traces.crs
gdf.to_file(f"{out_files_path}/sz_discretized_polygons.geojson", driver="GeoJSON", engine="fiona")

pkl.dump(discretized_dict, open(f"{out_files_path}/sz_discretized_dict.pkl", "wb"))
