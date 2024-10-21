
import os
import geopandas as gpd
from pcdhm.shared_helper import get_figure_bounds
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import matplotlib.ticker as mticker
import pickle as pkl
from pcdhm.shared_helper import get_probability_color, get_plot_order_list

##############

def get_mean_disp_barchart_data(site_PPE_dictionary, probability, exceed_type, site_list):
    """returns displacements at the X% probabilities of exceedance for each site
    This is effectively " find the X value for the desired Y"
    :param exceed_type: Options are "total_abs", "up", "down"
    """

    # get disp at 10% exceedance
    disps = []
    errs_plus = []
    errs_minus = []
    for site in site_list:
        threshold_vals = site_PPE_dictionary[site]["threshold_vals"]

        # displacement thresholds are negative for "down" exceedances
        if exceed_type == "down":
            threshold_vals = -threshold_vals

        # probability values at each threshold
        site_mean_probs = site_PPE_dictionary[site][f"weighted_exceedance_probs_{exceed_type}"]
        max_probs = site_PPE_dictionary[site][f"{exceed_type}_max_vals"]
        min_probs = site_PPE_dictionary[site][f"{exceed_type}_min_vals"]

        # get first index (x-value) for the probability (y-value) that is *just* < prob% (ideally we would interpolate
        # for exact value but don't have a function for that)
        mean_exceedance_index = next((index for index, value in enumerate(site_mean_probs) if value <= probability), -1)
        max_exceedance_index = next((index for index, value in enumerate(max_probs) if value <= probability), -1)
        min_exceedance_index = next((index for index, value in enumerate(min_probs) if value <= probability), -1)

        # displacement value at the X% probability
        disp = threshold_vals[mean_exceedance_index]
        disps.append(disp)

        #minimum and maximum values at the same index
        max_disp = threshold_vals[max_exceedance_index]
        min_disp = threshold_vals[min_exceedance_index]
        if exceed_type == "down":
            err_plus = abs(min_disp - disp)
            err_minus = abs(disp - max_disp)
        if exceed_type == "up" or exceed_type == "total_abs":
            err_plus = abs(max_disp - disp)
            err_minus = abs(disp - min_disp)

        errs_plus.append(err_plus)
        errs_minus.append(err_minus)

    return disps, errs_plus, errs_minus

def get_mean_prob_plot_data(site_PPE_dictionary, threshold, exceed_type, site_list):
    """ function that finds the probability at each site for the specified displacement threshold on the hazard curve
        Inputs:
        :param: site_PPE_dictionary, dictionary of exceedance probabilities for each site (key = site)
        :param exceedance type: string; "total_abs", "up", or "down"
        :param: site_list: list of sites to get data for. If None, will get data for all sites in site_PPE_dictionary.
                I made this option so that you could skip the sites you didn't care about (e.g., use "plot_order")

        Outputs:
        :return    probs: list of probabilities of exceeding the specified threshold (one per site)
        :return    errs_plus: list of (+) errors for each probability (one per site)
        :return    errs_minus: list of (-) errors for each probability (one per site)
            """

    # get disp at 10% exceedance
    probs = []
    errs_plus = []
    errs_minus = []

    # get list of probabilities at defined displacement threshold (one for each site)
    for site in site_list:
        site_PPE = site_PPE_dictionary[site][f"weighted_exceedance_probs_{exceed_type}"]
        threshold_vals = list(site_PPE_dictionary[site]["threshold_vals"])

        # find index in threshold_vals where the value matches the parameter threshold
        probs_index = threshold_vals.index(threshold)
        probs.append(site_PPE[probs_index])

        mean_prob = site_PPE[probs_index]
        max_prob = site_PPE_dictionary[site][f"{exceed_type}_max_vals"][probs_index]
        min_prob = site_PPE_dictionary[site][f"{exceed_type}_min_vals"][probs_index]

        err_plus = max_prob - mean_prob
        err_minus = mean_prob - min_prob
        errs_plus.append(err_plus)
        errs_minus.append(err_minus)

    return probs, errs_plus, errs_minus

def plot_weighted_mean_haz_curves(weighted_mean_PPE_dictionary, PPE_dictionary, exceed_type_list,
                                  model_version_title, out_directory, file_type_list, slip_taper, skipped_sites):
    """
    Intending for a single fault model (e.g., crustal CFM)
    Plots the weighted mean hazard curve for each site, for each exceedance type (total_abs, up, down)
    :param weighted_mean_PPE_dictionary: dictionary containing the weighted mean exceedance probabilities for each site.
    :param PPE_dictionary: dictionary containing the weighted mean exceedance probabilities for each branch
    :param exceed_type_list: list of strings, either "total_abs", "up", or "down"
    :return:
    """

    plot_order = get_plot_order_list(dictionary=weighted_mean_PPE_dictionary, skipped_sites=skipped_sites)

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    unique_id_list = list(PPE_dictionary.keys())

    for exceed_type in exceed_type_list:
        plt.close("all")
        fig, axs = plt.subplots(figsize=(8, 10))
        plt.subplots_adjust(hspace=0.3, wspace=0.3)

        # shade the region between the max and min value of all the curves at each site
        for i, site in enumerate(plot_order):
            ax = plt.subplot(4, 3, i + 1)

            # plots all three types of exceedance (total_abs, up, down) on the same plot
            max_probs = weighted_mean_PPE_dictionary[site][f"{exceed_type}_max_vals"]
            min_probs = weighted_mean_PPE_dictionary[site][f"{exceed_type}_min_vals"]
            threshold_vals = weighted_mean_PPE_dictionary[site]["threshold_vals"]
            threshold_vals = threshold_vals[1:]
            max_probs = max_probs[1:]
            min_probs = min_probs[1:]

            ax.fill_between(threshold_vals, max_probs, min_probs, color='0.9')

        # plot all the branches as light grey lines
        # for each branch, plot the exceedance probabilities for each site
        for k, unique_id in enumerate(unique_id_list):
            # this loop isn't really needed, but it's useful if you calculate Green's functions
            # at more sites than you want to plot

            for i, site in enumerate(plot_order):
                threshold_vals = PPE_dictionary[unique_id]["cumu_PPE_dict"][site]["thresholds"]
                site_exceedance_probs = PPE_dictionary[unique_id]["cumu_PPE_dict"][site][f"exceedance_probs_{exceed_type}"]

                ax = plt.subplot(4, 3, i + 1)
                ax.plot(threshold_vals, site_exceedance_probs, color='0.7')

        # loop through sites and add the weighted mean lines
        for i, site in enumerate(plot_order):
            ax = plt.subplot(4, 3, i + 1)

            # plots all three types of exceedance (total_abs, up, down) on the same plot
            weighted_mean_exceedance_probs = weighted_mean_PPE_dictionary[site][f"weighted_exceedance_probs_{exceed_type}"]
            threshold_vals = weighted_mean_PPE_dictionary[site]["threshold_vals"]

            line_color = get_probability_color(exceed_type)
            ax.plot(threshold_vals, weighted_mean_exceedance_probs, color=line_color, linewidth=2)

            ax.axhline(y=0.02, color="g", linestyle='dashed')
            ax.axhline(y=0.1, color="g", linestyle='dotted')

            xmin, xmax = 0.01, 3
            ymin, ymax = 0.000005, 1
            ax.set_title(site)
            ax.set_yscale('log'), ax.set_xscale('log')
            ax.set_ylim([ymin, ymax])
            ax.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
            ax.get_xaxis().set_major_formatter(ScalarFormatter())
            ax.ticklabel_format(axis='x', style='plain')
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.set_xlim([xmin, xmax])

        fig.text(0.5, 0, 'Vertical displacement threshold (m)', ha='center')
        fig.text(0, 0.5, 'Probability of exceedance in 100 years', va='center', rotation='vertical')
        fig.suptitle(f"weighted mean hazard curves\n{model_version_title} {taper_extension}\n{exceed_type}")
        plt.tight_layout()

        if not os.path.exists(f"../{out_directory}/weighted_mean_figures"):
            os.mkdir(f"../{out_directory}/weighted_mean_figures")

        for file_type in file_type_list:
            plt.savefig(
                f"../{out_directory}/weighted_mean_figures/weighted_mean_hazcurve_{exceed_type}{taper_extension}.{file_type}", dpi=300)

    # make a second graph with just the shaded envelope and weighted mean lines
    if len(exceed_type_list) > 1:
        plt.close("all")
        fig, axs = plt.subplots(figsize=(8, 10))
        plt.subplots_adjust(hspace=0.3, wspace=0.3)

        for i, site in enumerate(plot_order):
            ax = plt.subplot(4, 3, i + 1)
            for exceed_type in exceed_type_list:
                fill_color = get_probability_color(exceed_type)

                weighted_mean_max_probs = weighted_mean_PPE_dictionary[site][f"{exceed_type}_max_vals"]
                weighted_mean_min_probs = weighted_mean_PPE_dictionary[site][f"{exceed_type}_min_vals"]
                threshold_vals = weighted_mean_PPE_dictionary[site]["threshold_vals"]
                ax.fill_between(threshold_vals, weighted_mean_max_probs, weighted_mean_min_probs, color=fill_color,
                                alpha=0.2)

            # plot solid lines on top of the shaded regions
            for exceed_type in exceed_type_list:
                line_color = get_probability_color(exceed_type)
                weighted_mean_exceedance_probs = weighted_mean_PPE_dictionary[site][
                    f"weighted_exceedance_probs_{exceed_type}"]
                ax.plot(threshold_vals, weighted_mean_exceedance_probs, color=line_color, linewidth=2)

            # add 10% and 2% lines
            ax.axhline(y=0.02, color="g", linestyle='dashed')
            ax.axhline(y=0.1, color="g", linestyle='dotted')

            # make axes pretty
            ax.set_title(site)
            ax.set_yscale('log'), ax.set_xscale('log')
            ax.set_ylim([0.000005, 1]), ax.set_xlim([0.01, 3])
            ax.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
            ax.get_xaxis().set_major_formatter(ScalarFormatter())
            ax.ticklabel_format(axis='x', style='plain')
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))

        fig.text(0.5, 0, 'Vertical displacement threshold (m)', ha='center')
        fig.text(0, 0.5, 'Probability of exceedance in 100 years', va='center', rotation='vertical')
        exceed_types_string = ", ".join(exceed_type_list)
        fig.suptitle(f"weighted mean hazard curves \n{model_version_title} {taper_extension} \n{exceed_types_string} ")
        plt.tight_layout()

        for file_type in file_type_list:
            plt.savefig(f"../{out_directory}/weighted_mean_figures/weighted_mean_hazcurves{taper_extension}"
                        f".{file_type}", dpi=300)

def plot_weighted_mean_haz_curves_colorful(weighted_mean_PPE_dictionary, PPE_dictionary, exceed_type_list,
                                           model_version_title, out_directory, file_type_list, slip_taper, file_name,
                                           string_list, skipped_sites):
    """
    Plots the weighted mean hazard curve for each site, for each exceedance type (total_abs, up, down)
    :param weighted_mean_PPE_dictionary: dictionary containing the weighted mean exceedance probabilities for each site.
    :param PPE_dictionary: dictionary containing the weighted mean exceedance probabilities for each branch
    :param exceed_type_list: list of strings, either "total_abs", "up", or "down"
    :return:
    """

    plot_order = get_plot_order_list(weighted_mean_PPE_dictionary, skipped_sites=skipped_sites)

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    unique_id_list = list(PPE_dictionary.keys())

    for exceed_type in exceed_type_list:
        plt.close("all")
        fig, axs = plt.subplots(sharex=True, sharey=True, figsize=(12, 10), )
        plt.subplots_adjust(hspace=0.3, wspace=0.3)
        subplot_indices = [1, 2, 3, 5, 6, 7, 9, 10, 11, 13, 14, 15]

        # shade the region between the max and min value of all the curves at each site
        for i, site in enumerate(plot_order):
            # ax = plt.subplot(4, 4, i + 1)
            ax = plt.subplot(4, 4, subplot_indices[i])

            # plots all three types of exceedance (total_abs, up, down) on the same plot
            max_probs = weighted_mean_PPE_dictionary[site][f"{exceed_type}_max_vals"]
            min_probs = weighted_mean_PPE_dictionary[site][f"{exceed_type}_min_vals"]
            threshold_vals = weighted_mean_PPE_dictionary[site]["threshold_vals"]

            threshold_vals = threshold_vals[1:]
            max_probs = max_probs[1:]
            min_probs = min_probs[1:]

            ax.fill_between(threshold_vals, max_probs, min_probs, color='0.9', label="_nolegend_")

        # plot all the branches as light grey lines
        # for each branch, plot the exceedance probabilities for each site
        # make a list of random colors in the length on unique_ids
        # random_colors = [
        #     (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1)) for i in range(len(unique_id_list))
        # ]
        special_colors = [(235 / 255, 97 / 255, 35 / 255), (64 / 255, 176 / 255, 166 / 255), (116 / 255, 43 / 255,
                                                                                         140 / 255)]

        for k, unique_id in enumerate(unique_id_list):
            # this loop isn't really needed, but it's useful if you calculate Green's functions
            # at more sites than you want to plot


            if string_list[0] in unique_id:
                line_color = special_colors[0]
                linewidth = 1
            elif string_list[1] in unique_id:
                line_color = special_colors[1]
                linewidth = 1
            else:
                line_color = special_colors[2]
                linewidth = 1

            for i, site in enumerate(plot_order):

                threshold_vals = PPE_dictionary[unique_id]["cumu_PPE_dict"][site]["thresholds"]
                site_exceedance_probs = PPE_dictionary[unique_id]["cumu_PPE_dict"][site][f"exceedance_probs_{exceed_type}"]

                # skip the 0 value in the list
                threshold_vals = threshold_vals[1:]
                site_exceedance_probs = site_exceedance_probs[1:]

                # ax = plt.subplot(4, 3, i + 1)
                ax = plt.subplot(4, 4, subplot_indices[i])

                #ax.plot(threshold_vals, site_exceedance_probs, color='0.7')
                ax.plot(threshold_vals, site_exceedance_probs, color=line_color, linewidth=linewidth)

        # loop through sites and add the weighted mean lines
        for i, site in enumerate(plot_order):
            # ax = plt.subplot(4, 3, i + 1)
            ax = plt.subplot(4, 4, subplot_indices[i])

            # plots all three types of exceedance (total_abs, up, down) on the same plot
            weighted_mean_exceedance_probs = weighted_mean_PPE_dictionary[site][f"weighted_exceedance_probs_{exceed_type}"]
            threshold_vals = weighted_mean_PPE_dictionary[site]["threshold_vals"]

            threshold_vals = threshold_vals[1:]
            weighted_mean_exceedance_probs = weighted_mean_exceedance_probs[1:]

            line_color = get_probability_color(exceed_type)
            ax.plot(threshold_vals, weighted_mean_exceedance_probs, color=line_color, linewidth=2)

            ax.axhline(y=0.02, color="0.3", linestyle='dashed')
            ax.axhline(y=0.1, color="0.3", linestyle='dotted')

            xmin, xmax = 0.01, 3
            ymin, ymax = 0.000005, 1
            ax.set_title(site)
            ax.set_yscale('log'), ax.set_xscale('log')
            ax.set_ylim([ymin, ymax])
            ax.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
            ax.get_xaxis().set_major_formatter(ScalarFormatter())
            ax.ticklabel_format(axis='x', style='plain')
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax.set_xlim([xmin, xmax])

        # fig.show()
        fig.legend(unique_id_list, loc='center right', fontsize="x-small")
        fig.subplots_adjust(left=0.0)

        fig.text(0.5, 0, 'Vertical displacement threshold (m)', ha='center')
        fig.text(0, 0.5, 'Probability of exceedance in 100 years', va='center', rotation='vertical')
        fig.suptitle(f"weighted mean hazard curves\n{model_version_title} {taper_extension}\n{exceed_type}")
        plt.tight_layout()

        if not os.path.exists(f"../{out_directory}/weighted_mean_figures"):
            os.mkdir(f"../{out_directory}/weighted_mean_figures")

        for file_type in file_type_list:
            plt.savefig(
                f"../{out_directory}/weighted_mean_figures/"
                f"{file_name}_weighted_mean_hazcurve_{exceed_type}{taper_extension}.{file_type}", dpi=300)



# def map_and_plot_probabilities(weighted_mean_PPE_dictionary, exceed_type, model_version_title, out_directory,
#                            file_type_list,
#                        slip_taper, file_name, skipped_sites=None, probability=0.1, grid=False):
#
#     """ plots the displacement threshold value (m) for x% probability of exceedance in 100 years
#     CAVEATS/choices:
#     - currently set up with two subplots, a 10% and 2% probability of exceedance
#     - exeed_type is a list of ["total_abs", "up", "down"]
#     - fault_type can be "crustal" or "sz"
#
#     use grid=True if you want to plot the gridded data as an image rather than a set of scatter points"""
#
#     plot_order = get_plot_order_list(weighted_mean_PPE_dictionary, skipped_sites=skipped_sites)
#
#     if slip_taper is True:
#         taper_extension = "_tapered"
#     else:
#         taper_extension = "_uniform"
#
#     # load data
#     site_coords = [weighted_mean_PPE_dictionary[site]["site_coords"] for site in plot_order]
#     x_data = [coords[0] for coords in site_coords]
#     y_data = [coords[1] for coords in site_coords]
#
#     map = gpd.read_file("../data/map/nz_coastline.geojson")
#     wellington_boundary = gpd.read_file("../data/wellington_region_boundary.geojson")
#     plate_boundary = gpd.read_file("../data/map/plate_boundary.geojson")
#
#     # get plot bounds (set for Wellington Region at the moment)
#     plot_xmin, plot_ymin, plot_xmax, plot_ymax, \
#         xmin_tick, xmax_tick, ymin_tick, ymax_tick, tick_separation \
#         = get_figure_bounds(extent="Wellington close")
#
#     plt.close("all")
#     # two part figure, plot on the left and map on the right
#     fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.5))
#
#     colors = make_qualitative_colormap("custom", len(plot_order))
#
#
#     ##### plot disps on the left
#     disps, errs_plus, errs_minus = \
#         get_mean_disp_barchart_data(site_PPE_dictionary=weighted_mean_PPE_dictionary, exceed_type=exceed_type,
#                                     site_list=plot_order, probability=probability)
#     bar_width = 0.6
#
#     # plot bars and error bars
#     x = np.arange(len(plot_order))
#     bars = axs[0].bar(x, disps, bar_width, color=colors, linewidth=0.5)
#     axs[0].errorbar(x, disps, yerr=[errs_minus, errs_plus], fmt='none', ecolor='0.6', capsize=3,
#                     linewidth=1, markeredgewidth=1)
#
#     # add zero line
#     axs[0].axhline(y=0, color="k", linewidth=0.5)
#
#     # add value labels to bars
#     label_offset = 3 * 0.05
#     label_size = 6
#     # add value labels to bars
#     for bar in bars:
#         bar_color = bar.get_facecolor()
#         axs[0].text(x=(bar.get_x() + bar.get_width() * 0.6), y=bar.get_height() + label_offset,
#                     s=round(bar.get_height(), 1), ha='left',
#                     va='center', color=bar_color, fontsize=label_size, fontweight='bold')
#     axs[0].set_ylim(0.0, 3.0)
#     axs[0].set_ylabel("Displacement (m)", fontsize=8)
#     axs[0].tick_params(axis='y', labelrotation=90, labelsize=6)
#     axs[0].tick_params(axis='x', labelrotation=90, labelsize=6)
#     axs[0].set_xticks(x, plot_order)
#
#     #### Format map subplot
#     map.plot(ax=axs[1], color="k", linewidth=0.5)
#     plate_boundary.plot(ax=axs[1], color="0.75", linewidth=1.0)
#     axs[1].set_xticks(np.arange(xmin_tick, xmax_tick, tick_separation))
#     axs[1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mE'))
#     axs[1].set_yticks(np.arange(ymin_tick, ymax_tick, tick_separation))
#     axs[1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mN'))
#     plt.setp(axs[1].get_yticklabels(), rotation=90, ha="center", rotation_mode="anchor")
#     axs[1].tick_params(axis="both", which='major', labelsize=6)
#     axs[1].set_xlim(plot_xmin, plot_xmax)
#     axs[1].set_ylim(plot_ymin, plot_ymax)
#     axs[1].set_aspect("equal")
#
#     # add site points
#     axs[1].scatter(x_data, y_data, s=20, c=colors, edgecolors='black', linewidth=0.5, zorder=2)
#
#     # set titles and stuff
#     probability_string = str(int(100 * probability))
#     fig.suptitle(f"Displacement at {probability_string}%\n{model_version_title} {exceed_type} (100 yrs)")
#     fig.tight_layout()
#
#     # make directory for hazard map if it doesn't exist
#     if not os.path.exists(f"../{out_directory}/weighted_mean_figures"):
#         os.mkdir(f"../{out_directory}/weighted_mean_figures")
#
#     for type in file_type_list:
#         fig.savefig(f"../{out_directory}/weighted_mean_figures/disp{probability_string}_plot_hazard_map_{file_name}."
#                     f"{type}", dpi=300)

def map_and_plot_probabilities(PPE_path, plot_name, title, outfile_directory, plot_order, threshold=0.2,
                                 labels_on=True, file_type_list=["png"], exceed_type="down", colorbar_max=None):
    """Makes a two-part plot with probability of exceeding a threshold on the left and a map of sites, colored by
    probability, on the right"""

    # Load saved PPE dictionary
    with open(f"{PPE_path}", "rb") as fid:
        PPE_dict = pkl.load(fid)

        # set x-axis plot order info for probability subplot
        probability_x_vals = np.arange(len(plot_order))  # the site label locations

    # load data from PPE dictionary
    site_coords = [PPE_dict[site]["site_coords"] for site in plot_order]
    map_x_data = [coords[0] for coords in site_coords]
    # Load other map data
    coastline = gpd.read_file("../data/map/nz_coastline.geojson")
    plate_boundary = gpd.read_file("../data/map/plate_boundary.geojson")

    # get plot bounds (set for Wellington Region at the moment)
    plot_xmin, plot_ymin, plot_xmax, plot_ymax, \
        xmin_tick, xmax_tick, ymin_tick, ymax_tick, tick_separation \
        = get_figure_bounds(extent="Wellington close")
    map_y_data = [coords[1] for coords in site_coords]

    # Set plot and map point size
    point_size = 30
    fontsize = 8
    labelsize = 6
    color_map = "plasma"

    # two part figure, plot on the left and map on the right
    plt.close("all")
    fig, axs = plt.subplots(1, 2, figsize=(6.5, 3.5))

    ##### plot disps on the left
    # find the maximum y value for error bars so that the plot can be scaled correctly
    mean_probs, errs_plus, errs_minus = \
        get_mean_prob_plot_data(site_PPE_dictionary=PPE_dict, exceed_type=exceed_type,
                                threshold=threshold, site_list=plot_order)
    max_errs_y = 0.02 + max([mean_probs[j] + errs_plus[j] for j in range(len(mean_probs))])
    if colorbar_max is None:
        max_prob_color_val = round(max(mean_probs), 2) + 0.01
    else:
        max_prob_color_val = colorbar_max

    # add point and error bars to plot
    axs[0].errorbar(probability_x_vals, mean_probs, yerr=[errs_minus, errs_plus], fmt='none', ecolor='0.4', capsize=4,
                    linewidth=1, markeredgewidth=1)
    axs[0].scatter(probability_x_vals, mean_probs, s=point_size, c=mean_probs, cmap=color_map, zorder=3,
                   edgecolors='k', linewidths=0.5, vmin=0, vmax=max_prob_color_val)

    ##### Plot probability sites on map
    if len(map_x_data) < 20:  # if there are less than 20 points, plot with black edges
        disp_map = axs[1].scatter(map_x_data, map_y_data, s=point_size, c=mean_probs, cmap=color_map,
                                  edgecolors='black', linewidth=0.5, zorder=2, vmin=0, vmax=max_prob_color_val)
    else:  # otherwise plot without black edges
        disp_map = axs[1].scatter(map_x_data, map_y_data, s=point_size, c=mean_probs, cmap=color_map,
                                  edgecolors=None, linewidth=0.5, zorder=2, vmin=0, vmax=max_prob_color_val)

    #### Format the left plot (probability graph)
    # add other plot data
    coastline.plot(ax=axs[1], color="k", linewidth=0.5)
    plate_boundary.plot(ax=axs[1], color="0.75", linewidth=1.0)

    # probability labels: multiply the probability values by 100 to get a percentage, round to nearest %
    labels = [f"{round(100 * prob)}%" for prob in mean_probs]
    # this adds a 0.03 offset to the y value of the label so that it doesn't overlap with the point
    label_y_vals = [prob + 0.03 for prob in mean_probs]
    if labels_on:
        for site, q in enumerate(probability_x_vals):
            axs[0].text(x=probability_x_vals[q], y=label_y_vals[q], s=labels[q],
                        horizontalalignment='center', fontsize=6, fontweight='bold')

    # format the left subplot axes, ticks, and labels
    if max_errs_y < 0.25:
        axs[0].set_ylim(0.0, 0.25)
    else:
        axs[0].set_ylim(0.0, max_errs_y)
    axs[0].tick_params(axis='x', labelrotation=90, labelsize=labelsize)
    axs[0].tick_params(axis='y', labelsize=labelsize)
    axs[0].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    axs[0].yaxis.set_major_locator(mticker.MultipleLocator(0.1))
    axs[0].set_xticks(probability_x_vals, plot_order, va='top', ha='center')
    axs[0].set_ylabel("Probabilty", fontsize=fontsize)
    axs[0].set_title(f"Probability of exceeding {threshold} m {exceed_type}", fontsize=fontsize)

    ### Format the right subplot (map)
    # axes, ticks, labels
    axs[1].set_title(f"Probability of exceeding {threshold} m {exceed_type}", fontsize=fontsize)
    axs[1].set_xticks(np.arange(xmin_tick, xmax_tick, tick_separation))
    axs[1].xaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mE'))
    axs[1].set_yticks(np.arange(ymin_tick, ymax_tick, tick_separation))
    axs[1].yaxis.set_major_formatter(mticker.FormatStrFormatter('%.f mN'))
    plt.setp(axs[1].get_yticklabels(), rotation=90, ha="center", rotation_mode="anchor")
    axs[1].tick_params(axis="both", which='major', labelsize=labelsize)
    axs[1].set_xlim(plot_xmin, plot_xmax)
    axs[1].set_ylim(plot_ymin, plot_ymax)
    axs[1].set_aspect("equal")

    # colorbar
    divider = make_axes_locatable(axs[1])
    cax1 = divider.append_axes('right', size='6%', pad=0.05)
    cbar1 = fig.colorbar(disp_map, cax=cax1, orientation='vertical')
    cbar1.set_label("Probability", fontsize=fontsize)
    cbar1.ax.tick_params(labelsize=labelsize)
    cbar1.ax.yaxis.set_ticks_position('right')
    cbar1.ax.yaxis.set_label_position('right')

    # make output directory folder if it doesn't already exist
    if not os.path.exists(f"../{outfile_directory}"):
        os.makedirs(f"../{outfile_directory}")

    # set titles and stuff
    fig.suptitle(f"Probability of exceeding {threshold}\n{title} {exceed_type} (100 yrs)")
    fig.tight_layout()

    # make directory for hazard map if it doesn't exist
    if not os.path.exists(f"../{outfile_directory}"):
        os.mkdir(f"../{outfile_directory}")

    for file_type in file_type_list:
        fig.savefig(f"../{outfile_directory}/map_prob_{plot_name}.{file_type}", dpi=300)

