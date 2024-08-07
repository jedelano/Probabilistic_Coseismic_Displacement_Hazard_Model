import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import os
from functools import reduce
import numpy as np
import pickle as pkl
from scipy.interpolate import griddata
import matplotlib.pyplot as plt


def get_plot_order_list(dictionary, skipped_sites=None):
    """takes a dictionary and returns a list of keys in the dictionary
    can remove sites names that you wish to skip, so they aren't plotted later
    orders them from west to east along map distance"""

    site_names = list(dictionary.keys())

    if skipped_sites is not None:
        for site in skipped_sites:
            site_names.remove(site)

    # some thing to do with ordering the sites based on dictionary information

    # this next part is a stop gap solution
    plot_order = ["Paraparaumu", "Porirua CBD north", "Porirua CBD south", "South Coast",
                      "Wellington Airport", "Wellington CBD", "Petone", "Seaview",  "Eastbourne",
                       "Turakirae Head", "Lake Ferry", "Cape Palliser", "Flat Point"]

    if skipped_sites is not None:
        for site in skipped_sites:
            plot_order.remove(site)

    return plot_order

def get_probability_color(exceed_type):
    """ exceed type can be "total_abs", "up", or "down
    """
    if exceed_type == "total_abs":
        color = "k"
    elif exceed_type == "up":
        color = (189/255, 0, 0)
    elif exceed_type == "down":
        color = (15/255, 72/255, 186/255)

    return color


def make_qualitative_colormap(name, length):
    from collections import namedtuple
    # set up custom color scheme. Uses tab20 as a base, moves the greens to the start (to break up the
    # green/brown/red requence), makes the brown darker, and the red redder.
    if name == "tab20_subbed":
        #up to 20 colors
        all_colors = plt.get_cmap('tab20b')(np.linspace(0, 1, 20))
        color_indices = list(range(0, 20, 2))
        colors = all_colors[color_indices[0:length]]
        blue1, blue2 = colors[0], colors[1]
        green1, green2 = colors[2], colors[3]
        new_brown1 = [140 / 255, 84 / 255, 0 / 255, 1]
        new_brown2 = [112 / 255, 74 / 255, 1 / 255, 1]
        new_dark_red = [133 / 255, 5 / 255, 0 / 255, 1]
        #colors[0], colors[1] = green1, green2
        #colors[2], colors[3] = blue1, blue2
        colors[4], colors[6] = new_brown2, new_dark_red

    elif name == "tab20_reordered":
        # up to 20 colors
        all_colors = plt.get_cmap('tab20b')(np.linspace(0, 1, 20))
        color_indices = list(range(0, 20, 2))
        colors = all_colors[color_indices[0:length]]
        blue1, blue2 = colors[0], colors[1]
        green1, green2 = colors[2], colors[3]
        #new_brown1 = [140 / 255, 84 / 255, 0 / 255, 1]
        new_brown2 = [112 / 255, 74 / 255, 1 / 255, 1]
        new_dark_red = [133 / 255, 5 / 255, 0 / 255, 1]
        colors[0], colors[1] = green1, green2
        colors[2], colors[3] = blue1, blue2
        colors[4], colors[6] = new_brown2, new_dark_red

    elif name == "custom":
        tab20b_all_colors = plt.get_cmap('tab20b')(np.linspace(0, 1, 20))
        tab20c_all_colors = plt.get_cmap('tab20c')(np.linspace(0, 1, 20))
        tab20_all_colors = plt.get_cmap('tab20')(np.linspace(0, 1, 20))
        colors = [tab20b_all_colors[0], tab20b_all_colors[2],
                  tab20c_all_colors[0], tab20c_all_colors[2],
                  tab20b_all_colors[4], tab20b_all_colors[6],
                  tab20_all_colors[10], tab20b_all_colors[10],
                  tab20c_all_colors[4], tab20c_all_colors[6],
                  [133 / 255, 5 / 255, 0 / 255, 1], tab20b_all_colors[14],
                  tab20b_all_colors[16], tab20b_all_colors[18]]
        colors = colors[0:length]

    elif name == "tol_muted_ordered":
        # up to about 8 colors, plus light grey and black
        cset = namedtuple('Mcset',
                          'indigo cyan teal green olive sand rose wine purple pale_grey black')
        colors = cset('#332288','#88CCEE', '#44AA99','#117733', '#999933', '#DDCC77','#CC6677',
                    '#882255', '#AA4499', '#DDDDDD', '#000000')

    return colors


def read_rupture_csv(csv_file: str):
    rupture_dict = {}
    with open(csv_file, "r") as fid:
        index_data = fid.readlines()
    for line in index_data[1:]:
        numbers = [int(num) for num in line.strip().split(",")]
        rupture_dict[numbers[0]] = np.array(numbers[2:])
    return rupture_dict

def read_rake_csv(csv_file: str):
    df = pd.read_csv(csv_file)
    rake_dict = {}
    for i, row in df.iterrows():
        fault_id = int(row["fault_id"])
        rake_dict[fault_id] = {"cfm_rake": row["cfm_rake"], "model1_rake": row["model1_rake"], "model2_rake": row["model2_rake"]}
    return rake_dict

def read_average_slip(csv_file: str):
    df = pd.read_csv(csv_file)
    slip_dic = {}
    for i, row in df.iterrows():
        slip_dic[i] = row["Average Slip (m)"]
    return slip_dic

def cross_3d(a, b):
    """
    Calculates cross product of two 3-dimensional vectors.
    """
    x = ((a[1] * b[2]) - (a[2] * b[1]))
    y = ((a[2] * b[0]) - (a[0] * b[2]))
    z = ((a[0] * b[1]) - (a[1] * b[0]))
    return np.array([x, y, z])
def check_triangle_normal(triangle_vertices):
    """ needed for correct dip direction. check if the normal vector is positive or negative. If negative,
    flip the order of the vertices. Triangle_vertices is a 3x3 array of the vertices of the triangle (3 vertices,
    each with xyz)"""

    vector_a = triangle_vertices[1] - triangle_vertices[0]
    vector_b = triangle_vertices[1] - triangle_vertices[2]
    cross_a_b_vector = cross_3d(vector_a, vector_b)
    # Ensure that normal always points down (for proper strike convention), if not, flip the order of the vertices
    # this results in greens functions where dip slip is reverse motion, which is what we want (positive rake).
    if cross_a_b_vector[-1] > 0:
        ordered_triangle_vertices = triangle_vertices[-1::-1]
    else:
        ordered_triangle_vertices = triangle_vertices
    return ordered_triangle_vertices

def make_total_slip_dictionary(gf_dict_pkl):
    """ calculates total greens function displacement using strike slip gf, dip slip gf, and rake value

    need to run the and crustal_discretized_gfs script first"""

    with open(gf_dict_pkl, "rb") as fid:
        gf_dict = pkl.load(fid)

    # Makes a new total gf displacement dictionary using rake
    gf_adjusted_dict = {}
    for i in gf_dict.keys():
        # greens functions are just for the vertical component
        ss_gf = gf_dict[i]["ss"]
        ds_gf = gf_dict[i]["ds"]
        rake = gf_dict[i]["rake"]

        site_name_list = gf_dict[i]["site_name_list"]
        site_coords = gf_dict[i]["site_coords"]

        # calculate combined vertical from strike slip and dip slip using rake
        # Note: for subduction, the ss disp is "0" and the rake is "90" artificially (see 05_sz_discretized_gfs scritps)
        combined_gf = np.sin(np.radians(rake)) * ds_gf + np.cos(np.radians(rake)) * ss_gf
        gf_adjusted_dict[i] = {"combined_gf": combined_gf, "site_name_list": site_name_list, "site_coords": site_coords}

    return gf_adjusted_dict


def tol_cset(colorset=None):
    """
    Discrete color sets for qualitative data.

    Define a namedtuple instance with the colors.
    Examples for: cset = tol_cset(<scheme>)
      - cset.red and cset[1] give the same color (in default 'bright' colorset)
      - cset._fields gives a tuple with all color names
      - list(cset) gives a list with all colors
    """
    from collections import namedtuple

    namelist = ('bright', 'high-contrast', 'vibrant', 'muted', 'medium-contrast', 'light')
    if colorset is None:
        return namelist
    if colorset not in namelist:
        colorset = 'bright'
        print('*** Warning: requested colorset not defined,',
              'known colorsets are {}.'.format(namelist),
              'Using {}.'.format(colorset))

    if colorset == 'bright':
        cset = namedtuple('Bcset',
                          'blue red green yellow cyan purple grey black')
        return cset('#4477AA', '#EE6677', '#228833', '#CCBB44', '#66CCEE',
                    '#AA3377', '#BBBBBB', '#000000')

    if colorset == 'high-contrast':
        cset = namedtuple('Hcset',
                          'blue yellow red black')
        return cset('#004488', '#DDAA33', '#BB5566', '#000000')

    if colorset == 'vibrant':
        cset = namedtuple('Vcset',
                          'orange blue cyan magenta red teal grey black')
        return cset('#EE7733', '#0077BB', '#33BBEE', '#EE3377', '#CC3311',
                    '#009988', '#BBBBBB', '#000000')

    if colorset == 'muted':
        cset = namedtuple('Mcset',
                          'rose indigo sand green cyan wine teal olive purple pale_grey black')
        return cset('#CC6677', '#332288', '#DDCC77', '#117733', '#88CCEE',
                    '#882255', '#44AA99', '#999933', '#AA4499', '#DDDDDD',
                    '#000000')

    if colorset == 'medium-contrast':
        cset = namedtuple('Mcset',
                          'light_blue dark_blue light_yellow dark_red dark_yellow light_red black')
        return cset('#6699CC', '#004488', '#EECC66', '#994455', '#997700',
                    '#EE99AA', '#000000')

    if colorset == 'light':
        cset = namedtuple('Lcset',
                          'light_blue orange light_yellow pink light_cyan mint pear olive pale_grey black')
        return cset('#77AADD', '#EE8866', '#EEDD88', '#FFAABB', '#99DDFF',
                    '#44BB99', '#BBCC33', '#AAAA00', '#DDDDDD', '#000000')

def filter_ruptures_by_rate(directory):
    """returns rupture indices that have annual rate >0"""

    # load individual files
    rates_df = pd.read_csv(f"../data/{directory}/solution/rates.csv")
    # only keep rupture scenarios with Annual Rates > 0
    trimmed_rates_df = rates_df[rates_df["Annual Rate"] > 0]

    filtered_ruptures = trimmed_rates_df.index.values.tolist()
    print(f"initial scenarios: {len(rates_df)}")
    print(f"rate filtered scenarios: {len(filtered_ruptures)}")
    return filtered_ruptures


# This runs a bit slowly
def filter_ruptures_by_location(NSHM_directory, target_rupture_ids, fault_type,
                                crustal_outfiles_path, sz_outfiles_path,
                                search_radius=2.5e5):
    """ filters the initial rupture scenarios by which patches are involved
        set a distance from interest area and cut out scenarios that don't intersect

        for now, rupture_df can be any list with the targeted rupture indices. For example, a list of ruptures
        that have been filtered by annual rate in the filter_ruptures_by_rate script above

        fault_type = "crustal" or "sz"
        """

    # input: interest location
    Wellington = Point(1749150, 5428092)
    if fault_type == "crustal":
        fault_rectangle_centroids_gdf = gpd.read_file(f"{crustal_outfiles_path}/named_rectangle_centroids.geojson")
    if fault_type == "sz":
        fault_rectangle_centroids_gdf = gpd.read_file(
            f"{sz_outfiles_path}/all_rectangle_centroids.geojson")

    all_ruptures_patch_indices = read_rupture_csv(f"../data/{NSHM_directory}/ruptures/indices.csv")

    # find rupture scenarios that match input target ruptures (e.g., from filter by rate)
    trimmed_rupture_patch_indices = {i: all_ruptures_patch_indices[i] for i in all_ruptures_patch_indices.keys() if i in
                                     target_rupture_ids}

    # find faults patches that are within search radius
    filtered_fault_ids = []
    for i in range(len(fault_rectangle_centroids_gdf.centroid)):
        centroid = fault_rectangle_centroids_gdf.centroid[i]
        if centroid.distance(Wellington) < search_radius:
            #filtered_fault_ids.append(patch_centroids_gdf.index[i])
            filtered_fault_ids.append(int(fault_rectangle_centroids_gdf.fault_id[i]))


    # this can probably be simplified
    # include scenarios that have those patches
    filtered_scenarios = []
    for rupture_index in target_rupture_ids:
        # uses scenarios that include any patch within that search radius
        if any(fault_id in filtered_fault_ids for fault_id in trimmed_rupture_patch_indices[rupture_index]):
            filtered_scenarios.append(rupture_index)

    print(f"location filtered scenarios: {len(filtered_scenarios)}")
    return filtered_scenarios


def calculate_vertical_disps(ruptured_discretized_polygons_gdf, ruptured_rectangle_outlines_gdf, rupture_id,
                             ruptured_fault_ids, slip_taper, rupture_slip_dict, gf_total_slip_dict):
    """ calcualtes displacements for given rupture scenario at a single site
    not yet sure if I should set it up to allow more than one site at a time

    CAVETS/choices:
    - tapered slip assigns one slip value to each discretized polygon (e.g., one fault id). the slip
    values are tapered according to total rupture length of all rectangles.
    - this version changes all very small dispalcements to zero, or if no meshes are used, returns zero displacement
    """

    # find which patches have a mesh and which don't, to use greens functions later just with meshed patches
    ruptured_fault_ids_with_mesh = np.intersect1d(ruptured_fault_ids, list(gf_total_slip_dict.keys()))

    # calculate slip on each discretized polygon
    if slip_taper is False:
        # calculate displacements by multiplying scenario slip by scenario greens function
        # scenario gf sums displacements from all ruptured
        gfs_i = np.sum([gf_total_slip_dict[j]["combined_gf"] for j in ruptured_fault_ids_with_mesh], axis=0)
        disps_scenario = rupture_slip_dict[rupture_id] * gfs_i
        polygon_slips = rupture_slip_dict[rupture_id] * np.ones(len(ruptured_fault_ids_with_mesh))

        # storing zeros is more efficient than nearly zeros. Makes v small displacements = 0
        if len(ruptured_fault_ids_with_mesh) != 0:
            disps_scenario[np.abs(disps_scenario) < 5.e-3] = 0.
        elif len(ruptured_fault_ids_with_mesh) == 0:
            disps_scenario = None

    elif slip_taper is True:
        # get centroid coords of faults discretized polygons with a mesh
        ruptured_polygon_centroid_points = ruptured_discretized_polygons_gdf.centroid
        ruptured_polygon_centroids_x = [point.x for point in ruptured_polygon_centroid_points]
        ruptured_polygon_centroids_y = [point.y for point in ruptured_polygon_centroid_points]
        ruptured_polygon_centroid_coords = np.array([ruptured_polygon_centroids_x, ruptured_polygon_centroids_y]).T

        # get bounds of fault patches, makes np array with 4 coords (minx, miny, maxx, maxy)
        rupture_bounds = ruptured_rectangle_outlines_gdf.total_bounds

        # makes 1000 points along a line between endpoints (bounds of fault rectangles).
        along_rupture_line_x = np.linspace(rupture_bounds[0], rupture_bounds[2], 1000)
        along_rupture_line_y = np.linspace(rupture_bounds[1], rupture_bounds[3], 1000)
        # stack into one column of xy pairs
        along_rupture_line_xy = np.column_stack((along_rupture_line_x, along_rupture_line_y))

        # calculate distance along line for each xy point
        start_point = Point(along_rupture_line_xy[0])
        line_distances = []
        for coord in along_rupture_line_xy:
            next_point = Point(coord)
            distance = start_point.distance(next_point)
            line_distances.append(distance)
        line_length = np.max(line_distances)

        # calculate slip at each interpolated point based on distance
        # this constant is based on the integral of the sin function from 0 to 1 (see NSHM taper)
        max_slip = rupture_slip_dict[rupture_id] / 0.76276
        # apply slip taper function to max slip. slip = sqrt(sin(pi * distance/line_length))
        # making a multiplier list is verbose but helps me keep track of things
        tapered_slip_multipliers = []
        tapered_slip_values = []
        for distance in line_distances:
            if np.sin(np.pi * distance / line_length) < 5.e-5:      # this is to fix error below of sqrt(0)
                tapered_slip_multiplier = 0.
            else:
                tapered_slip_multiplier = np.sqrt(np.sin(np.pi * distance / line_length))
            tapered_slip_multipliers.append(tapered_slip_multiplier)
            tapered_slip_values.append(max_slip * tapered_slip_multiplier)

        # interpolate slip at each discretized polygon (i.e., patch) centroid and corresponding displacement
        polygon_slips = griddata(along_rupture_line_xy, tapered_slip_values, ruptured_polygon_centroid_coords,
                               method="nearest")

        # calculate displacements by multiplying the polygon green's function by slip on each fault
        # this will be a list of lists
        disps_i_list = []
        for i, fault_id in enumerate(ruptured_discretized_polygons_gdf.fault_id):
            disp_i = gf_total_slip_dict[fault_id]["combined_gf"] * polygon_slips[i]
            disps_i_list.append(disp_i)
        #sum displacements from each patch
        disps_scenario = np.sum(disps_i_list, axis=0)
        if len(ruptured_fault_ids_with_mesh) != 0:
            disps_scenario[np.abs(disps_scenario) < 5.e-3] = 0.
        elif len(ruptured_fault_ids_with_mesh) == 0:
            disps_scenario = None


    return disps_scenario, polygon_slips

def get_rupture_disp_dict(NSHM_directory, fault_type, naming_extension, slip_taper, gf_name,
                          results_version_directory, crustal_outfiles_path, sz_outfiles_path):
    """
    inputs: uses extension naming scheme to load NSHM rate/slip data and fault geometry, state slip taper

    extracts site and rupture data, pares down the ruptures based on location (within a buffer) and annual rate >0,
    and passes that to the get_displacements function

    outputs: a dictionary where each key is the rupture id in a pickle file. contains displacements (length = same as
    greens function), annual rate (single value for the rupture), site name list (should be same length as green's
    function), and site coordinates (same length as site name list)

    CAVEATS:
    - current version omits scenarios from the output list if all locations have zero displacement
    - a little clunky because most of the dictionary columns are repeated across all keys.
    """

    # load saved data
    print(f"\nloading data for {naming_extension}")

    rupture_slip_dict = read_average_slip(f"../data/{NSHM_directory}/ruptures/average_slips.csv")
    rates_df = pd.read_csv(f"../data/{NSHM_directory}/solution/rates.csv")

    if fault_type == "crustal":
        discretized_polygons_gdf = gpd.read_file(f"{crustal_outfiles_path}/crustal_discretized_polygons.geojson")
        gf_dict_pkl = f"{crustal_outfiles_path}/crustal_gf_dict_{gf_name}.pkl"
        rectangle_outlines_gdf = gpd.read_file(f"{crustal_outfiles_path}/all_rectangle_outlines.geojson")
    elif fault_type == "sz":
        discretized_polygons_gdf = gpd.read_file(f"{sz_outfiles_path}/sz_discretized_polygons.geojson")
        gf_dict_pkl = f"{sz_outfiles_path}/sz_gf_dict_{gf_name}.pkl"
        rectangle_outlines_gdf = gpd.read_file(f"{sz_outfiles_path}/all_rectangle_outlines.geojson")

    # this line takes ages, only do it once
    all_ruptures = read_rupture_csv(f"../data/{NSHM_directory}/ruptures/indices.csv")

    # for some reason it defaults values to string. Convert to integer.
    discretized_polygons_gdf['fault_id'] = discretized_polygons_gdf['fault_id'].astype('int64')
    rectangle_outlines_gdf['fault_id'] = rectangle_outlines_gdf['fault_id'].astype('int64')

    # filter ruptures by annual rate and location
    filtered_ruptures_annual_rate = filter_ruptures_by_rate(NSHM_directory)
    filtered_ruptures_location = filter_ruptures_by_location(NSHM_directory=NSHM_directory,
                                                             target_rupture_ids=filtered_ruptures_annual_rate,
                                                             fault_type=fault_type,
                                                             crustal_outfiles_path=crustal_outfiles_path,
                                                             sz_outfiles_path=sz_outfiles_path)

    # Makes a new total gf displacement dictionary using rake. If points don't have a name (e.g., for whole map
    # calculations), the site name list is just a list of numbers
    # Note: rake is artificially set to "90" for subduction. See 05_sz_gfs scripts
    gf_total_slip_dict = make_total_slip_dictionary(gf_dict_pkl)
    first_key = list(gf_total_slip_dict.keys())[0]
    site_name_list = gf_total_slip_dict[first_key]["site_name_list"]
    site_coords = gf_total_slip_dict[first_key]["site_coords"]

    # calculate displacements at all the sites by rupture. Output dictionary keys are by rupture ID.
    disp_dictionary = {}
    print(f"calculating displacements for {naming_extension}")
    for rupture_id in filtered_ruptures_location:
        ruptured_fault_ids = all_ruptures[rupture_id]
        ruptured_discretized_polygons_gdf = discretized_polygons_gdf[
            discretized_polygons_gdf.fault_id.isin(ruptured_fault_ids)]
        ruptured_discretized_polygons_gdf = gpd.GeoDataFrame(ruptured_discretized_polygons_gdf, geometry='geometry')
        ruptured_rectangle_outlines_gdf = rectangle_outlines_gdf[
            rectangle_outlines_gdf.fault_id.isin(ruptured_fault_ids)]

        #calculate displacements, output is a list of displacements for each site
        disps_scenario, patch_slips = \
            calculate_vertical_disps(ruptured_discretized_polygons_gdf=ruptured_discretized_polygons_gdf,
                                     ruptured_rectangle_outlines_gdf=ruptured_rectangle_outlines_gdf,
                                     rupture_id=rupture_id, ruptured_fault_ids=ruptured_fault_ids,
                                     slip_taper=slip_taper, rupture_slip_dict=rupture_slip_dict,
                                     gf_total_slip_dict=gf_total_slip_dict)

        # extract annual rate and save data to dictionary. Key is the rupture ID. Ignores scenarios with zero
        # displacement at all sites. Sites can be a grid cell or a specific (named) site.
        if disps_scenario is not None:
            annual_rate = rates_df[rates_df.index == rupture_id]["Annual Rate"].values[0]
            # displacement dictionary for a single rupture scenario at all sites. Key is rupture id.
            rupture_disp_dict = {"rupture_id": rupture_id, "v_disps_m": disps_scenario, "annual_rate": annual_rate,
                                 "site_name_list": site_name_list, "site_coords": site_coords,
                                 "x_data": site_coords[:, 0], "y_data": site_coords[:, 1],
                                 "polygon_slips_m": patch_slips}
            disp_dictionary[rupture_id] = rupture_disp_dict

    # print statement about how many scenarios have displacement > 0 at each site
    print(f"scenarios with displacement > 0: {len(disp_dictionary)}")

    if slip_taper is True:
        extension3 = "_tapered"
    else:
        extension3 = "_uniform"

    # save displacements
    if not os.path.exists(f"../{results_version_directory}/{naming_extension}"):
        os.mkdir(f"../{results_version_directory}/{naming_extension}")

    with open(f"../{results_version_directory}/{naming_extension}/all_rupture_disps_{naming_extension}{extension3}.pkl",
              "wb") as f:
        pkl.dump(disp_dictionary, f)

    return disp_dictionary

def get_figure_bounds(polygon_gdf="", extent=""):
    """sets figure bounds based on key words
    polygon_gdf: either discretized polygon gdf (for displacement figure) or ruptured rectangles gdf (slip figure)
    extent: can specify the extent of figure for interest area """

    if extent == "North Island":    # bounds of whole north island
        buffer = 100000
        plot_xmin, plot_ymin, plot_xmax, plot_ymax = 1525000, 5270000, 2300000, 6176000
        xmin_tick, xmax_tick  = 1600000, plot_xmax
        ymin_tick, ymax_tick = 5400000, 6176000
        tick_separation = 300000
    elif extent == "Wellington":    # bounds around the wellington region, with a 10 km buffer
        x, y, buffer = 1771150, 5428092, 10.e4
        plot_xmin, plot_ymin, plot_xmax, plot_ymax = x - buffer, y - buffer, x + buffer, y + buffer
        tick_separation = 100000.
        xmin_tick, xmax_tick = 1700000, 1800000 + tick_separation / 4
        ymin_tick, ymax_tick = 5380000, 5480000 + tick_separation / 4
    elif extent == "Wellington close":  # bounds around the wellington region, with a 10 km buffer
        x, y, buffer = 1771150, 5428092, 5.e4
        plot_xmin, plot_ymin, plot_xmax, plot_ymax = x - buffer*0.9, y - buffer, x + buffer*1.7, y + buffer
        tick_separation = 50000.
        xmin_tick, xmax_tick = 1750000, 1850000 + tick_separation / 4
        ymin_tick, ymax_tick = 5400000, 5480000 + tick_separation / 4
    elif extent == "ruptured_rectangles":   # intended for ruptured rectangles gdf (slip plot)
        buffer = 20000
        plot_xmin = polygon_gdf.total_bounds[0] - buffer
        plot_ymin = polygon_gdf.total_bounds[1] - buffer
        plot_xmax = polygon_gdf.total_bounds[2] + buffer
        plot_ymax = polygon_gdf.total_bounds[3] + buffer
        xmin_tick, xmax_tick = round(plot_xmin + buffer, -4), plot_xmax
        ymin_tick, ymax_tick = round(plot_ymin + buffer, -4), plot_ymax
        tick_separation = round((plot_ymax - plot_ymin) / 3, -4)
    else:   # bounds of all polyons with a 100 km buffer (intented for discretized polygons, displacement plot)
        plot_xmin, plot_ymin, plot_xmax, plot_ymax = polygon_gdf.total_bounds
        xmin_tick, xmax_tick = round(plot_xmin, -5) - 100000, round(plot_xmax, -5) + 100000
        ymin_tick, ymax_tick = round(plot_ymin, -5) - 100000, round(plot_ymax, -5) + 100000
        tick_separation = 400000.
    return plot_xmin, plot_ymin, plot_xmax, plot_ymax, xmin_tick, xmax_tick, ymin_tick, ymax_tick, tick_separation

def save_target_rates(NSHM_directory, target_rupture_ids, extension1, results_version_directory):
    """get the annual rates from NSHM solution for target ruptures, output a csv file

    NSHM_directory = name of NSHM folder
    target_rupture_ids = list of rupture ids/indices
    out_directory = directory for all the other outputfiles (figures, etc.) """

    # load annual rate file
    rates_df = pd.read_csv(f"../data/{NSHM_directory}/solution/rates.csv")
    # only keep ruptures and rates of interest
    trimmed_rates_df = rates_df[rates_df.index.isin(target_rupture_ids)]

    if not os.path.exists(f"../{results_version_directory}/{extension1}"):
        os.mkdir(f"../{results_version_directory}/{extension1}")

    trimmed_rates_df.to_csv(f"../{results_version_directory}/{extension1}/ruptures_rates.csv", sep=',')
    print(f"{extension1} annual rates written to .csv")

def make_branch_weight_dict(branch_weight_file_path, sheet_name):
    """
    This function reads in the excel file with the branch weights and returns a dictionary with the branch weights
    and other information (scaling values, solution file names, etc.).
    The dictionary keys are the unique ID strings based on the branch parameters

    :param branch_weight_file_path: string; path to the excel file with the branch weights
    :param sheet_name: string; name of the sheet in the excel file with the branch weights
    """

    # read in the Excel file with the branch weights and other metadata
    branch_weights = pd.read_excel(branch_weight_file_path, sheet_name=sheet_name, header=0)

    # make a dictionary with the branch weights and other metadata
    branch_weight_dict = {}
    for row in range(len(branch_weights)):

        N_val = branch_weights["N"][row]
        N_string = str(N_val).replace('.', '')
        b_val = branch_weights["b"][row]
        b_string = str(b_val).replace('.', '')
        C_val = branch_weights["C"][row]
        C_string = str(C_val).replace('.', '')
        S_val = branch_weights["S"][row]
        S_string = str(S_val).replace('.', '')
        def_model  = branch_weights["def_model"][row]
        time_dependence = branch_weights["time_dependence"][row]
        file_suffix = branch_weights["PCDHM_file_suffix"][row]
        total_weight_RN = branch_weights["total_weight_RN"][row]

        # make a unique ID for each branch.
        # The NSHM solution files in 'data' directory do not include the rate scaling factor (S)
        # (i.e., they are all S=1)
        # this simplifies how many solution results are stored and processed
        # These lines use the same solution results (and file suffix) for 3 different S values
        unique_id = f"N{N_string}_b{b_string}_C{C_string}_S{S_string}_{time_dependence}_{def_model}{file_suffix}"

        branch_weight_dict[unique_id] = {"N": N_val, "b": b_val, "C": C_val, "S": S_val,
                                         "def_model": def_model,
                                         "time_dependence": time_dependence,
                                         "file_suffix": file_suffix,
                                         "total_weight_RN": total_weight_RN}

    return branch_weight_dict