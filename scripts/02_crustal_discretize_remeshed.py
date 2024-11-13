import geopandas as gpd
import pandas as pd
import numpy as np
import meshio
from shapely.geometry import Polygon, Point
import pickle as pkl
from pcdhm.shared_helper import read_rake_csv, check_triangle_normal

# This script discretizes the fault meshes into patches based on the "sheet in the wind" NSHM faults.
# only run once per fault geometry (can reuse for any branch) because all the faults are the same

###### Inputs
# Doesn't really matter which inversion solution because all the NSHM fault files are the same.
NSHM_directory = "NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDEz"
# provide model extension to match the mesh directory and name output directory
crustal_mesh_version = "_CFM"    #_CFM, _Model1, or _Model2
crustal_model_directory = "_CFM_tapered" # can be the same as the mesh version, add "taperd" if tapered slip.
out_files_directory = "mesh_gf_outfiles_r1"

######## script
# this must be the same length as the number of meshes and have some value that matches all the target fault sections
# will need to come up with a better way to do this in the future when more faults/meshes are used
target_NSHM_fault_names = ["Aotea|Evans Bay", "Dry River|Huangarua", "Fisherman", "Honeycomb",
                           "Moonshine|Otaki", "Ohariu", "Okupe", "Opouawe", "Otaraia",
                           "Pahaua", "Palliser|Kaiwhata", "Pukerua", "Riversdale","Mana|Otaheke",
                           "Wairarapa", "Wellington Hutt", "Whareama", "Wharekauhau", "Whitemans"]

mesh_directory = f"../data/wellington_alt_geom/meshes{crustal_mesh_version}/STL_remeshed"

out_files_path = f"../{out_files_directory}/crustal{crustal_model_directory}"
target_traces = gpd.read_file(f"{out_files_path}/name_filtered_fault_sections.geojson" "").to_crs(
    epsg=2193)
all_traces = gpd.read_file(f"../data/crustal_solutions/{NSHM_directory}/ruptures/fault_sections.geojson"
                       "").to_crs(epsg=2193)
# read in rake data
rake_dict = read_rake_csv(f"../data/wellington_alt_geom/alt_geom_rakes.csv")
if crustal_mesh_version == "_Model2":
    rake_col = "model2_rake"
elif crustal_mesh_version == "_Model1":
    rake_col = "model1_rake"
elif crustal_mesh_version == "_CFM":
    rake_col = "cfm_rake"

# read in fault meshes and assign unique names to each one
# there is probably a better way to do this, but the names matter so this is it for now.
AoteaEvansBay_mesh = meshio.read(f"{mesh_directory}/Aotea-EvansBay_remeshed.stl")
DryRiver_Huangarua_mesh = meshio.read(f"{mesh_directory}/DryRiver-Huangaruacombined_remeshed.stl")
Fisherman_mesh = meshio.read(f"{mesh_directory}/Fishermancombined_remeshed.stl")
Honeycomb_mesh = meshio.read(f"{mesh_directory}/Honeycomb_remeshed.stl")
Moonshine_OtakiForks_mesh = meshio.read(f"{mesh_directory}/Moonshine-OtakiForkscombined_remeshed.stl")
Ohariu_mesh = meshio.read(f"{mesh_directory}/Ohariucombined_remeshed.stl")
Okupe_mesh = meshio.read(f"{mesh_directory}/Okupe1_remeshed.stl")
OpouaweUruti_mesh = meshio.read(f"{mesh_directory}/Opouawe-Uruti_remeshed.stl")
Otaraia_mesh = meshio.read(f"{mesh_directory}/Otaraia_remeshed.stl")
Pahaua_mesh = meshio.read(f"{mesh_directory}/Pahaua_remeshed.stl")
PalliserKaiwhata_mesh = meshio.read(f"{mesh_directory}/Palliser-Kaiwhata_remeshed.stl")
PukeruaShepherdsGully_mesh = meshio.read(f"{mesh_directory}/Pukerua-ShepherdsGullycombined_remeshed.stl")
Riversdale_mesh = meshio.read(f"{mesh_directory}/Riversdale_remeshed.stl")
ShepherdsGullyMana_Otaheke_mesh = meshio.read(f"{mesh_directory}/ShepherdsGully-Mana-Otahekecombined_remeshed.stl")
Wairarapa_mesh = meshio.read(f"{mesh_directory}/Wairarapacombined_remeshed.stl")
WellingtonHutt_mesh = meshio.read(f"{mesh_directory}/WellingtonHuttValleycombined_remeshed.stl")
WhareamaBank_mesh = meshio.read(f"{mesh_directory}/WhareamaBank_remeshed.stl")
Wharekauhau_mesh = meshio.read(f"{mesh_directory}/Wharekauhau_remeshed.stl")
WhitemansValley_mesh = meshio.read(f"{mesh_directory}/WhitemansValley_remeshed.stl")

#make mesh list with same order as fault name list
mesh_list = [AoteaEvansBay_mesh, DryRiver_Huangarua_mesh, Fisherman_mesh, Honeycomb_mesh, Moonshine_OtakiForks_mesh,
                Ohariu_mesh, Okupe_mesh, OpouaweUruti_mesh, Otaraia_mesh, Pahaua_mesh, PalliserKaiwhata_mesh,
                PukeruaShepherdsGully_mesh, Riversdale_mesh, ShepherdsGullyMana_Otaheke_mesh, Wairarapa_mesh,
                WellingtonHutt_mesh, WhareamaBank_mesh, Wharekauhau_mesh, WhitemansValley_mesh]
# one unique name per fault mesh. Matches the names and order of the mesh list.
mesh_name_list = ["AoteaEvansBay_mesh", "DryRiver_Huangarua_mesh", "Fisherman_mesh", "Honeycomb_mesh",
                  "Moonshine_OtakiForks_mesh", "Ohariu_mesh", "Okupe_mesh", "OpouaweUruti_mesh", "Otaraia_mesh",
                  "Pahaua_mesh", "PalliserKaiwhata_mesh", "PukeruaShepherdsGully_mesh", "Riversdale_mesh",
                  "ShepherdsGullyMana_Otaheke_mesh", "Wairarapa_mesh", "WellingtonHutt_mesh", "WhareamaBank_mesh",
                  "Wharekauhau_mesh", "WhitemansValley_mesh"]
#########

# Sensitivity testing for patch edge depth/dip angle. Skip for crustal for now.
steeper_dip = False
gentler_dip = False

####################

named_rectangle_centroids = []
named_rectangle_polygons = []
all_rectangle_polygons = []
# make dataframe for fault section rectangle attributes (for export later)
df_named_rectangle = pd.DataFrame(columns=['fault_id', 'fault_name', 'dip_deg', 'patch_height_m', 'up_depth_km',
                                  'low_depth_km', 'dip_dir_deg'])
df_all_rectangle = pd.DataFrame(columns=['fault_id', 'fault_name', 'dip_deg', 'patch_height_m', 'up_depth_km',
                                  'low_depth_km', 'dip_dir_deg'])

df_named_rectangle_centroid = pd.DataFrame(columns=['fault_id', 'fault_name', 'dip_deg', 'patch_height_m',
                                                    'up_depth_km','low_depth_km', 'dip_dir_deg'])
df_all_rectangle_centroid = pd.DataFrame(columns=['fault_id', 'fault_name'])

# turn NSHM traces into rectangular patches using the metadata in the GeoJSON file
for i, trace in all_traces.iterrows():
    # Convert the trace to a numpy array
    trace_xy_array = np.array(trace.geometry.coords)
    # NSHM crustal doesn't have z coords on the traces, need to add in as zero.
    trace_array = np.zeros(shape =(2, 3))
    # pull first and last vertex from the trace
    trace_array[0] = np.append(trace_xy_array[0], trace.UpDepth)
    trace_array[1] = np.append(trace_xy_array[-1], trace.UpDepth)
    # Convert the depths to metres (not needed here since it's zero, but good practice)
    trace_array[:, -1] *= -1000.

    low_depth, up_depth = trace.LowDepth, trace.UpDepth
    # Calculate the centroid of the trace
    trace_centroid = np.array([*trace.geometry.centroid.coords[0], trace.UpDepth * -1000])
    # Calculate the height of the patch
    rectangle_height = (trace.LowDepth - trace.UpDepth) * 1000.
    dip_angle = trace.DipDeg
    extension2 = ""

    # write patch attributes to dictionary and add to bottom of data frame
    df2_rectangle = pd.DataFrame({'fault_id': [int(trace.FaultID)], 'fault_name': [trace.FaultName],
                                  'dip_deg': [dip_angle], 'rectangle_patch_height_m': [rectangle_height],
                                  'up_depth_km': [up_depth], 'low_depth_km': [low_depth],
                                  'dip_dir_deg': [trace.DipDir]}, index=[i])
    df_all_rectangle = pd.concat([df_all_rectangle, df2_rectangle])

    # # Calculate the across-strike vector, (orthogonal to strike) aka dip direction xy vector
    across_strike = np.array([np.sin(np.radians(trace.DipDir)), np.cos(np.radians(trace.DipDir))])

    # Calculate the down-dip vector by incorporating the dip angle
    cosdip = np.cos(np.radians(dip_angle))
    sindip = np.sin(np.radians(dip_angle))
    down_dip = np.array([cosdip * across_strike[0], cosdip * across_strike[1], - 1 * sindip])


    # Calculate the corners of the patch and make a shapely polygon
    dd1 = trace_array[1] + (rectangle_height / sindip) * down_dip
    dd2 = trace_array[0] + (rectangle_height / sindip) * down_dip
    rectangle_polygon = Polygon(np.vstack([trace_array, dd1, dd2]))

    # Append the patch and polygon to lists
    all_rectangle_polygons.append(rectangle_polygon)
####
all_rectangle_outline_gs = gpd.GeoSeries(all_rectangle_polygons, crs=2193)
all_rectangle_outline_gdf = gpd.GeoDataFrame(df_all_rectangle, geometry=all_rectangle_outline_gs.geometry, crs=2193)
all_rectangle_outline_gdf.to_file(f"{out_files_path}/all_rectangle_outlines.geojson",
                              driver="GeoJSON")

# Turn name-filtered section traces into rectangular patches using the metadata in the GeoJSON file
for i, trace in target_traces.iterrows():
    # Convert the trace to a numpy array
    trace_xy_array = np.array(trace.geometry.coords)
    # NSHM crustal doesn't have z coords on the traces, need to add in as zero.
    trace_array = np.zeros(shape =(2, 3))
    # pull first and last vertex from the trace
    trace_array[0] = np.append(trace_xy_array[0], trace.UpDepth)
    trace_array[1] = np.append(trace_xy_array[-1], trace.UpDepth)
    # Convert the depths to metres (not needed here since it's zero, but good practice)
    trace_array[:, -1] *= -1000.

    # dip/depth sensitivity testing
    initial_rectangle_height = (trace.LowDepth - trace.UpDepth) * 1000.
    initial_tandip = np.tan(np.radians(trace.DipDeg))
    rectangle_horizontal_dist = initial_rectangle_height / initial_tandip

    # set conditions for rectangular patch geometry parameters, which change with dip angle
    # for steeper dip ( fault section depths * 1.85)
    # probably is a little bit broken if steeper or gentler are True
    if steeper_dip == True and gentler_dip == False:
        low_depth, up_depth = trace.LowDepth * 1.15, trace.UpDepth * 1.15
        rectangle_height = (low_depth - up_depth) * 1000.
        trace_centroid = np.array([*trace.geometry.centroid.coords[0], trace.UpDepth * -1000])
        dip_angle = np.degrees(np.arctan(rectangle_height / rectangle_horizontal_dist))
        crustal_mesh_version += "_steeperdip"
    # for gentler dip ( fault section depths * 0.85)
    elif gentler_dip == True and steeper_dip == False:
        low_depth, up_depth = trace.LowDepth * 0.85, trace.UpDepth * 0.85
        rectangle_height = (low_depth - up_depth) * 1000.
        trace_centroid = np.array([*trace.geometry.centroid.coords[0], trace.UpDepth * -1000])
        dip_angle = np.degrees(np.arctan(rectangle_height / rectangle_horizontal_dist))
        crustal_mesh_version += "_gentlerdip"
    # uses geometry from NSHM with no modifications
    elif gentler_dip is False and steeper_dip is False:
        low_depth, up_depth = trace.LowDepth, trace.UpDepth
        # Calculate the centroid of the trace
        trace_centroid = np.array([*trace.geometry.centroid.coords[0], trace.UpDepth * -1000])
        # Calculate the height of the patch
        rectangle_height = (trace.LowDepth - trace.UpDepth) * 1000.
        dip_angle = trace.DipDeg
        crustal_mesh_version += ""
    else:
        print("Dip modifications are wrong. Only one statement can be True at once. Try again.")
        exit()

    # write patch attributes to dictionary and add to bottom of data frame
    df2_rectangle = pd.DataFrame({'fault_id': [int(trace.FaultID)], 'fault_name': [trace.FaultName],
                                  'dip_deg': [dip_angle], 'rectangle_patch_height_m': [rectangle_height],
                                  'up_depth_km': [up_depth], 'low_depth_km': [low_depth],
                                  'dip_dir_deg': [trace.DipDir]}, index=[i])
    df_named_rectangle = pd.concat([df_named_rectangle, df2_rectangle])

    df2_rectangle_centroid = pd.DataFrame({'fault_id': [trace.FaultID], 'fault_name': [trace.FaultName]}, index=[i])
    df_named_rectangle_centroid = pd.concat([df_named_rectangle_centroid, df2_rectangle_centroid])
    #######################
    # # Calculate the across-strike vector, (orthogonal to strike) aka dip direction xy vector
    across_strike = np.array([np.sin(np.radians(trace.DipDir)), np.cos(np.radians(trace.DipDir))])

    # Calculate the down-dip vector by incorporating the dip angle
    cosdip = np.cos(np.radians(dip_angle))
    sindip = np.sin(np.radians(dip_angle))
    down_dip = np.array([cosdip * across_strike[0], cosdip * across_strike[1], - 1 * sindip])

    # Calculate the centroid of the patch
    rectangle_centroid = trace_centroid + (rectangle_height / sindip) / 2 * down_dip

    # Calculate the corners of the patch and make a shapely polygon
    dd1 = trace_array[1] + (rectangle_height / sindip) * down_dip
    dd2 = trace_array[0] + (rectangle_height / sindip) * down_dip
    rectangle_polygon = Polygon(np.vstack([trace_array, dd1, dd2]))

    # Append the patch centroid and polygon to lists
    named_rectangle_centroids.append(rectangle_centroid)
    named_rectangle_polygons.append(rectangle_polygon)

# write patch centroid and rectangle polygons to geojson
# probably an opportunity to reduce code length by combining this with above
named_rectangle_centroids = np.array(named_rectangle_centroids)
named_rectangle_centroids_gs = gpd.GeoSeries([Point(centroid) for centroid in named_rectangle_centroids], crs=2193)
named_rectangle_centroids_gdf = gpd.GeoDataFrame(
    df_named_rectangle_centroid, geometry=named_rectangle_centroids_gs.geometry, crs=2193)
named_rectangle_centroids_gdf.to_file(
    f"{out_files_path}/named_rectangle_centroids.geojson", driver="GeoJSON")

named_rectangle_outline_gs = gpd.GeoSeries(named_rectangle_polygons, crs=2193)
named_rectangle_outline_gdf = gpd.GeoDataFrame(df_named_rectangle, geometry=named_rectangle_outline_gs.geometry, crs=2193)
named_rectangle_outline_gdf.to_file(
    f"{out_files_path}/named_rectangle_polygons.geojson", driver="GeoJSON")


#### mesh stuff
#make output geodata frame to add all discretized patch/fault info to
out_gdf = gpd.GeoDataFrame({"rake":[], "geometry":[], "fault_name":[], "fault_id":[]})

# set up dictionaries. Key will be the fault id
discretised_dict = {}       #values are xyz or triangle coords for each fault ID (i.e., rectangle)

# loop over individual fault meshes, find closest rectangle patch centroid for same fault name
#######important to subset by name so you don't get the closest patch centroid from a different fault
for i in range(len(mesh_list)):
    # fault mesh name, deal with some mesh at a time
    mesh = mesh_list[i]
    mesh_name = mesh_name_list[i]
    print(f"making discretized mesh for {mesh_name}")
    mesh_triangles_indices = mesh.cells_dict["triangle"]    # indices of vertices that make up triangles

    # NOT CURRENTLY USING FOR CRUSTAL. Multiply depths by steeper/gentler constant
    if steeper_dip == True:
        mesh_vertices = mesh.points * [1, 1, 1.15]
    elif gentler_dip == True:
        mesh_vertices = mesh.points * [1, 1, 0.85]
    else:
        mesh_vertices = mesh.points # xyz of vertices in mesh. indexed.

    # array of 3 xyz arrays. (three sets of vertices to make a triangle)
    triangle_vertex_arrays = mesh_vertices[mesh_triangles_indices]

    ###### ensuring correct strike convention for accurate displacment calculations later
    # calculate the normal vector to each triangle. If negative, reverse ordering of triangle of vertices.
    ordered_triangle_vertex_arrays = []
    for triangle in triangle_vertex_arrays:
        ordered_triangle_array = check_triangle_normal(triangle)
        ordered_triangle_vertex_arrays.append(ordered_triangle_array)
    ordered_triangle_vertex_arrays = np.array(ordered_triangle_vertex_arrays)

    triangle_centroids = np.mean(ordered_triangle_vertex_arrays, axis=1)   # three part arrray. means of x, y, and z.

    # find the closest fault rectangle patch to each triangle centroid
    closest_rectangles = []
    closest_rectangle_fault_ids = []    # index (fault ID) for closest rectangle centroid to each triangle centroid
    for triangle_centroid in triangle_centroids:
        #search_terms = target_NSHM_fault_names[i]
        # subset centroids by fault name that matches mesh name (prevents grabbing the wrong fault)
        rectangle_centroid_gdf_i = named_rectangle_centroids_gdf[
            named_rectangle_centroids_gdf['fault_name'].str.contains(target_NSHM_fault_names[i], case=False)]
        # get coordinates of patch centroids
        # needs to be in array format to work with triangle centroids
        patch_centroids_i = np.vstack([np.array(value.coords) for value in
                                       rectangle_centroid_gdf_i.geometry.values])

        # find closest patch centroid to each triangle
        distances = np.linalg.norm(patch_centroids_i - triangle_centroid, axis=1)
        if distances.min() < 10.e4:
            closest_patch = np.argmin(distances)    # gets list index from subset of patches
            # finds fault id of closest patch within centroids that are subset by name
            closest_fault_id = rectangle_centroid_gdf_i.fault_id.values[closest_patch]
            #closest_fault_id = named_rectangle_centroids_gdf.fault_id.values[closest_patch]
            closest_rectangle_fault_ids.append(closest_fault_id)
        else:
            closest_rectangle_fault_ids.append(-1) # ??

    # array of: index (fault ID) for closest rectangle centroid to each triangle centroid
    closest_rectangle_fault_ids = np.array(closest_rectangle_fault_ids)

    # Create fault polygons from triangles
    mesh_polygons = []

    # for adding to discretized polygon geojson attributes
    polygon_rakes = []
    polygon_fault_ids = []

    # make discretized patches (patch of triangles)
    for j, trace in target_traces.iterrows():
        # get the triangles that are closest to the trace of interest (i.e., index of closest rectangle matches fault
        # index)
        triangles_locs = np.where(closest_rectangle_fault_ids == trace.FaultID)
        triangles = ordered_triangle_vertex_arrays[triangles_locs]

        # skip the trace if there's no matching mesh in the rake dictionary
        trace_id = trace.FaultID
        if trace_id in rake_dict.keys():
            # extract rake value from rake dictionary based on FaultID
            rake = rake_dict[trace.FaultID][rake_col]

            # skip discretizing a rectangular patch if we don't have a mesh that matches it
            if len(triangles) > 0:
                # make dictionary of triangles that go with each polygon
                triangle_polygons = [Polygon(triangle) for triangle in triangles]
                dissolved_polygons = gpd.GeoSeries(triangle_polygons).union_all()
                mesh_polygons.append(dissolved_polygons)

                # Find the rake for each patch
                polygon_rakes.append(rake)

                # add patch index
                polygon_fault_ids.append(trace.FaultID)

                # set the triangles and rake for each fault ID
                discretised_dict[trace.FaultID] = {"triangles": triangles, "rake": rake}

    # get patch rakes and fault ids to add to output geojson file
    polygon_rakes = np.array(polygon_rakes)
    polygon_fault_ids = np.array(polygon_fault_ids)
    mesh_name_array = np.full(polygon_fault_ids.shape, mesh_name)

    # Create a geodataframe from the discretized polygons
    gdf = gpd.GeoDataFrame({"rake": polygon_rakes, "geometry": mesh_polygons, "fault_name": mesh_name_array, "fault_id":
        polygon_fault_ids})
    out_gdf = pd.concat([out_gdf, gdf])

# make pickle with triangle vertices
pkl.dump(discretised_dict, open(f"{out_files_path}/crustal_discretized_dict.pkl", "wb"))

# write discretized polygons to geojson
out_gdf.crs = target_traces.crs
out_gdf.to_file(f"{out_files_path}/crustal_discretized_polygons.geojson",
                driver="GeoJSON")







