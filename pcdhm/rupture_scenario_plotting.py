
import geopandas as gpd
import os
import shutil
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
import matplotlib as mpl
from pcdhm.shared_helper import read_rupture_csv, make_total_slip_dictionary, get_figure_bounds


# not using right now for crustal faults
def plot_vert_difference(NSHM_directory, rupture_slip_dict, target_rupture_ids, patch_polygons, gf_dict_pkl_coast1,
                         gf_dict_pkl_coast2,
                         extension1, extension2):
    """calculates the vertical difference between two greens functions for specific ruptures.
    use for comparing the same rupture scanrio with different geometry or other parameter"""

    # Read in rupture data and map elements
    all_ruptures = read_rupture_csv(f"../data/{NSHM_directory}/ruptures/indices.csv")
    #rupture_slip_dict = read_average_slip(f"{NSHM_directory}/ruptures/average_slips.csv")
    coastline = gpd.read_file("../data/map/nz_coastline.geojson")
    patch_polygons_gdf = gpd.read_file(patch_polygons)

    # Make outfile directory, or delete existing one with same name
    if os.path.exists(f"out_files/{extension1}{extension2}/figures/vert_difference"):
        shutil.rmtree(f"out_files/{extension1}{extension2}/figures/vert_difference")
    os.mkdir(f"out_files/{extension1}{extension2}/figures/vert_difference")

    # Make slip dictionaries for the two coast greens functions (e.g., steeperdip and gentler dip)
    gf_total_slip_dict1 = make_total_slip_dictionary(gf_dict_pkl_coast1)
    gf_total_slip_dict2 = make_total_slip_dictionary(gf_dict_pkl_coast2)

    # Set plot bounds
    x, y, buffer = 1749150, 5428092, 7.e4
    plot_xmin, plot_ymin, plot_xmax, plot_ymax = x - buffer, y - buffer, x + buffer, y + buffer

    # load coastaline xy data
    x_data = np.load("xpoints_coast.npy")
    y_data = np.load("ypoints_coast.npy")

    # Make plot
    for rupture_id in target_rupture_ids:
        plt.close("all")
        ruptured_patches = all_ruptures[rupture_id]
        # sum greens function for all patches to get scenario gf
        gfs1_i = np.sum([gf_total_slip_dict1[j] for j in ruptured_patches], axis=0)
        gfs2_i = np.sum([gf_total_slip_dict2[j] for j in ruptured_patches], axis=0)
        # calculate disps by multiplying scenario avg slip by scenario greens function [X disp, Y disp, VS]
        disps_scenario1 = rupture_slip_dict[rupture_id] * gfs1_i
        disps_scenario2 = rupture_slip_dict[rupture_id] * gfs2_i
        # storing zeros is more efficient than nearly zeros. Makes v small disps = 0
        disps_scenario1[np.abs(disps_scenario1) < 5.e-3] = 0.
        disps_scenario2[np.abs(disps_scenario2) < 5.e-3] = 0.
        # calculate difference between two variations
        disps_scenario_diff = disps_scenario2 - disps_scenario1

        max_vdisp = np.max(np.abs(disps_scenario_diff[:, -1]))
        if max_vdisp < 0.1:
            max_vdisp_cbar = 0.1
        else: max_vdisp_cbar = max_vdisp

        fig, axs = plt.subplots(1, 2, figsize=(6.5, 5))

        # plot ruptured patches in grey
        patch_polygons_gdf[patch_polygons_gdf.index.isin(ruptured_patches)].plot(ax=axs[0], color="0.5")

        # plot vertical diff (last column in disps array)
        disps = axs[1].scatter(x_data, y_data, s=4, c=disps_scenario_diff[:, -1],
                              cmap="seismic", vmin=-max_vdisp_cbar, vmax=max_vdisp_cbar)

        for ax in axs:
            coastline.plot(ax=ax, color="k", linewidth=0.5)
            ax.set_xticks(np.arange(1600000., 1900000., 50000.))
            ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mE'))
            ax.set_yticks(np.arange(5300000., 5600000., 50000.))
            ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mN'))
            plt.setp(ax.get_yticklabels(), rotation=90, ha="center", rotation_mode="anchor")
            ax.tick_params(axis="both", which='major', labelsize=6)
            ax.set_xlim(plot_xmin, plot_xmax)
            ax.set_ylim(plot_ymin, plot_ymax)
            ax.set_aspect("equal")

        divider = make_axes_locatable(axs[0])
        cax1 = divider.append_axes('left', size='5%', pad=0.05)
        cax1.set_visible(False)

        divider2 = make_axes_locatable(axs[1])
        cax2 = divider2.append_axes('right', size='5%', pad=0.05)
        cbar2 = fig.colorbar(disps, cax=cax2, orientation='vertical')
        cbar2.set_label("Vertical deformation difference (m)", fontsize=8)
        cbar2.ax.tick_params(labelsize=6)

        fig.suptitle(f"Rupture {rupture_id}{extension2}")
        fig.tight_layout()
        fig.savefig(f"out_files/{extension1}{extension2}/figures/vert_difference/coast_rupture_{rupture_id}.png",
                    dpi=300)
        print(f"Coastal def difference {rupture_id}")


def vertical_disp_figure(all_ruptures_disp_dict, NSHM_directory, target_rupture_ids, naming_extension, slip_taper,
                         grid, fault_type, results_version_directory, mesh_version,
                         crustal_outfiles_path, sz_outfiles_path,
                         extent="Wellington", file_type_list=["png", "pdf"]):

    """ makes two-part figure with ruptured patches on left and vertical deformation on right
                """
    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    # load data for plotting
    NSHM_traces_gdf = gpd.read_file(f"../data/{NSHM_directory}/ruptures/fault_sections.geojson"
                                    "").to_crs(epsg=2193)
    if fault_type == "crustal":
        discretized_polygons_gdf = gpd.read_file(f"{crustal_outfiles_path}/{fault_type}_discretized_polygons.geojson")
        rectangle_outlines_gdf = gpd.read_file(f"{crustal_outfiles_path}/all_rectangle_outlines.geojson")
        part_a_figure_extent = "ruptured_rectangles"
    elif fault_type == "sz":
        discretized_polygons_gdf = gpd.read_file(f"{sz_outfiles_path}/{fault_type}_discretized_polygons.geojson")
        rectangle_outlines_gdf = gpd.read_file(f"{sz_outfiles_path}/all_rectangle_outlines.geojson")
        part_a_figure_extent = "North Island"
    # for some reason it defaults values to string. Convert to integer.
    discretized_polygons_gdf['fault_id'] = discretized_polygons_gdf['fault_id'].astype('int64')
    rectangle_outlines_gdf['fault_id'] = rectangle_outlines_gdf['fault_id'].astype('int64')
    NSHM_traces_gdf['FaultID'] = NSHM_traces_gdf['FaultID'].astype('int64')

    coastline = gpd.read_file("../data/map/nz_coastline.geojson")
    plate_boundary = gpd.read_file("../data/map/plate_boundary.geojson")

    # load rupture and corresponding fault index data for extracting which faults to plot
    all_ruptures = read_rupture_csv(f"../data/{NSHM_directory}/ruptures/indices.csv")

    # get displacement figure bounds
    plot_b_xmin, plot_b_ymin, plot_b_xmax, plot_b_ymax, \
        xmin_b_tick, xmax_b_tick, ymin_b_tick, ymax_b_tick, tick_b_separation \
        = get_figure_bounds(extent=extent, polygon_gdf=discretized_polygons_gdf)

    ### make figure
    for rupture_id in target_rupture_ids:
        plt.close("all")
        ruptured_fault_ids = all_ruptures[rupture_id]
        # extract fault outline data using ruptured fault ids
        ruptured_discretized_polygons_gdf = discretized_polygons_gdf[
            discretized_polygons_gdf.fault_id.isin(ruptured_fault_ids)]
        ruptured_discretized_polygons_gdf = gpd.GeoDataFrame(ruptured_discretized_polygons_gdf, geometry='geometry')
        ruptured_rectangle_outlines_gdf = rectangle_outlines_gdf[rectangle_outlines_gdf.fault_id.isin(ruptured_fault_ids)]
        ruptured_NSHM_traces_gdf = NSHM_traces_gdf[NSHM_traces_gdf.FaultID.isin(ruptured_fault_ids)]

        # set plot bounds and ticks for slip figure (part a)
        plot_a_xmin, plot_a_ymin, plot_a_xmax, plot_a_ymax, xmin_a_tick, xmax_a_tick, ymin_a_tick, ymax_a_tick, \
            tick_a_separation = get_figure_bounds(ruptured_rectangle_outlines_gdf, extent=part_a_figure_extent)

        # extract needed data from displacement dictionary
        disps_scenario = all_ruptures_disp_dict[rupture_id]["v_disps_m"]
        patch_slips = all_ruptures_disp_dict[rupture_id]["polygon_slips_m"]
        plot_x_data = all_ruptures_disp_dict[rupture_id]["x_data"]
        plot_y_data = all_ruptures_disp_dict[rupture_id]["y_data"]

        max_vert_disp = np.max(np.abs(disps_scenario))
        ruptured_discretized_polygons_gdf["slip"] = patch_slips

        if slip_taper is False:
            max_slip_color_val = np.max(patch_slips) / 0.76276
        else:
            max_slip_color_val = np.max(patch_slips)

        # format fig
        fig, axs = plt.subplots(1, 2, figsize=(6.5, 5))

        # ## plot ruptured patches and slip on the discretized polygons
        ruptured_discretized_polygons_gdf.plot(ax=axs[0], column=ruptured_discretized_polygons_gdf["slip"],
                                               cmap="viridis", vmin=0, vmax=max_slip_color_val)
        ruptured_rectangle_outlines_gdf.boundary.plot(ax=axs[0], linewidth= 0.5, color="0.5")
        ruptured_discretized_polygons_gdf.boundary.plot(ax=axs[0], linewidth= 0.5, color="0.2")
        ruptured_NSHM_traces_gdf.plot(ax=axs[0], linewidth= 0.5, color="r")

        if grid is True:
            # plot displacement as a grid
            # get x and y dimensions
            length_unique_x = len(np.unique(plot_x_data))
            length_unique_y = len(np.unique(plot_y_data))
            # reshape list back into a grid (for plotting)
            disps_scenario = np.reshape(disps_scenario, (length_unique_y, length_unique_x))
            plot_x_data = np.reshape(plot_x_data, (length_unique_y, length_unique_x))
            plot_y_data = np.reshape(plot_y_data, (length_unique_y, length_unique_x))

            disp_map = axs[1].imshow(disps_scenario[-1::-1], cmap="seismic", vmin=-3.0, vmax=3.0,
                                  extent=[plot_x_data.min(), plot_x_data.max(), plot_y_data.min(), plot_y_data.max()])
            # add def contours, levels can be int or list/array
            if disps_scenario.min() < -0.5:
                axs[1].contour(plot_x_data, plot_y_data, disps_scenario, levels=[-3, -2.5, -2, -1.5, -1, -0.5],
                                             colors="steelblue", linewidths=1, linestyles='dashed')
            if disps_scenario.max() > 1.0:
                axs[1].contour(plot_x_data, plot_y_data, disps_scenario, levels=[1, 2, 3],
                                             colors="maroon", linewidths=1)
            axs[1].contour(plot_x_data, plot_y_data, disps_scenario, levels=[0],
                           colors="0.5", linewidths=1)

        # plot point data using scatter
        else:
            if len(plot_x_data) < 20:  # if there are less than 20 points, plot with black edges
                disp_map = axs[1].scatter(plot_x_data, plot_y_data, s=15, c=disps_scenario, cmap="seismic",
                                                    edgecolors='black', linewidth=0.5, zorder=2, vmin=-3.0, vmax=3.0)
            else:  # otherwise plot without black edges
                disp_map = axs[1].scatter(plot_x_data, plot_y_data, s=15, c=disps_scenario, cmap="seismic",
                                                    edgecolors=None, linewidth=0.5, zorder=2, vmin=-3.0, vmax=3.0)

        # Format subplots
        for ax in axs:
            coastline.plot(ax=ax, color="k", linewidth=0.5)
            plate_boundary.plot(ax=ax, color="0.75", linewidth=1.0)
            ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mE'))
            ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mN'))
            plt.setp(ax.get_yticklabels(), rotation=90, ha="center", rotation_mode="anchor")
            ax.tick_params(axis="both", which='major', labelsize=6)
            ax.set_aspect("equal")

        # set left plot limits (slip)
        axs[0].set_xlim(plot_a_xmin, plot_a_xmax)
        axs[0].set_ylim(plot_a_ymin, plot_a_ymax)
        axs[0].set_xticks(np.arange(xmin_a_tick, xmax_a_tick, tick_a_separation))
        axs[0].set_yticks(np.arange(ymin_a_tick, ymax_a_tick, tick_a_separation))
        axs[0].set_aspect("equal")

        #set right plot limits (displacement)
        axs[1].set_xlim(plot_b_xmin, plot_b_xmax)
        axs[1].set_ylim(plot_b_ymin, plot_b_ymax)
        axs[1].set_xticks(np.arange(xmin_b_tick, xmax_b_tick, tick_b_separation))
        axs[1].set_yticks(np.arange(ymin_b_tick, ymax_b_tick, tick_b_separation))
        axs[1].set_aspect("equal")

        # color bars, labels, and stuff
        divider = make_axes_locatable(axs[0])
        cax1 = divider.append_axes('top', size='6%', pad=0.05)
        cbar1 = fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=max_slip_color_val)),
                             cax=cax1, orientation='horizontal')
        cbar1.set_label("slip (m)", fontsize=8)
        cbar1.ax.tick_params(labelsize=6)
        cbar1.ax.xaxis.set_ticks_position('top')
        cbar1.ax.xaxis.set_label_position('top')

        divider2 = make_axes_locatable(axs[1])
        cax2 = divider2.append_axes('top', size='5%', pad=0.05)
        cbar2 = fig.colorbar(disp_map, cax=cax2, orientation='horizontal')
        cbar2.set_label("Vertical deformation (m)", fontsize=8)
        cbar2.ax.tick_params(labelsize=6)
        cbar2.ax.xaxis.set_ticks_position('top')
        cbar2.ax.xaxis.set_label_position('top')

        fig.suptitle(f"Rupture {rupture_id}{taper_extension}\n{mesh_version}")
        fig.tight_layout()

        # make folder for outputs
        if not os.path.exists(f"../{results_version_directory}/{naming_extension}"):
            os.mkdir(f"../{results_version_directory}/{naming_extension}")
        if not os.path.exists(f"../{results_version_directory}/{naming_extension}/scenario_displacements"):
            os.mkdir(f"../{results_version_directory}/{naming_extension}/scenario_displacements")
        if not os.path.exists(f"../{results_version_directory}/{naming_extension}/rupture_geojson"):
            os.mkdir(f"../{results_version_directory}/{naming_extension}/rupture_geojson")

        # save figure and geojson with slip data
        for file_type in file_type_list:
            fig.savefig(f"../{results_version_directory}/{naming_extension}/scenario_displacements/"
                        f"rupture_{rupture_id}{taper_extension}.{file_type}", dpi=300)
        ruptured_discretized_polygons_gdf.to_file(
            f"../{results_version_directory}/{naming_extension}/rupture_geojson/scenario_slip_{rupture_id}{taper_extension}.geojson",
            driver="GeoJSON")


