import os
import pickle as pkl
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib.patches import Rectangle
from pcdhm.shared_helper import make_qualitative_colormap
from run_mean_bar_charts import get_mean_prob_barchart_data #get_mean_disp_barchart_data



plot_order = ["Paraparaumu", "Porirua CBD north", "South Coast", "Wellington Airport", "Wellington CBD", "Petone",
                  "Seaview", "Eastbourne", "Turakirae Head", "Lake Ferry", "Cape Palliser", "Flat Point"]

# names = ["CFM mean uniform", "CFM mean tapered"]
# mean_PPE1_path = "../results_jde/crustal_CFM/weighted_mean_PPE_dict_testing_uniform.pkl"
# mean_PPE2_path = "../results_jde/crustal_CFM/weighted_mean_PPE_dict_testing_tapered.pkl"
#plot_name = "CFM_mean_uniform_tapered"

names = ["CFM mean uniform", "Alt. Fault 1 mean uniform"]
mean_PPE_paths = ["../results_jde/crustal_CFM/weighted_mean_PPE_dict_testing_uniform.pkl",
                   "../results_jde/crustal_Model1/weighted_mean_PPE_dict_testing_uniform.pkl"]
plot_name = "CFM_mean_uniform_alt1_mean_uniform"


outfile_directory = "../results_jde/compare_fault_models"

title = "CFM mean uniform vs Alt. Fault 1 mean uniform"


def get_probability_plot_data(site_PPE_dictionary, exceed_type, threshold, site_list=None):
    # I can't remember why I put a site list parameter.
    # define exceedance type. Options are "total_abs", "up", "down"

    if site_list == None:
        site_list = list(site_PPE_dictionary.keys())

    # get probability at defined displacement threshold
    probs_threshold = []

    for site in site_list:
        site_PPE = site_PPE_dictionary[site][f"exceedance_probs_{exceed_type}"]
        threshold_vals = list(site_PPE_dictionary[site]["thresholds"])

        # find index in threshold_vals where the value matches the parameter threshold
        index = threshold_vals.index(threshold)

        probs_threshold.append(site_PPE[index])

    return probs_threshold


def get_exceedance_plot_data(site_PPE_dictionary, probability, exceed_type, site_list=None):
    """returns displacements at the X% probabilities of exceedance for each site

    define exceedance type. Options are "total_abs", "up", "down"
    """

    if site_list == None:
        site_list = list(site_PPE_dictionary.keys())

    # get disp at 10% exceedance
    disps = []
    for site in site_list:
        threshold_vals = site_PPE_dictionary[site]["thresholds"]

        # dispalcement threasholds are negative for "down" exceedances
        if exceed_type == "down":
            threshold_vals = -threshold_vals

        site_PPE = site_PPE_dictionary[site][f"exceedance_probs_{exceed_type}"]

        # get first index that is < 10% (ideally we would interpolate for exaxt value but don't have a function)
        exceedance_index = next((index for index, value in enumerate(site_PPE) if value <= probability), -1)
        disp = threshold_vals[exceedance_index]
        disps.append(disp)

    return disps


# Chris pointed out that using a bar chart of this isn't really a good idea
def make_branch_compare_prob_barchart(crustal_file_suffix, sz_file_suffix, slip_taper,
                              results_directory="results_files", threshold=0.2):
    """ What is the probability of exceeding 0.2 m subsidence, 0.2 m uplift at each site?
    pickle files are a string with the path
    """


    exceed_type_list = ["up", "down"]

    if slip_taper:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"


    with open(f"../{results_directory}/sites_{crustal_file_suffix}/crustal_cumu_exceed_prob_"
              f"{crustal_file_suffix}{taper_extension}.pkl", "rb") as fid:
        crustal_site_PPE_dictionary = pkl.load(fid)
    with open(f"../{results_directory}/sites_{sz_file_suffix}/sz_cumu_exceed_prob_"
              f"{sz_file_suffix}.pkl", "rb") as fid:
        sz_site_PPE_dictionary = pkl.load(fid)

    # set up custom color scheme
    colors = make_qualitative_colormap("custom", len(plot_order))

    fig, axs = plt.subplots(1, 2, figsize=(7, 4))
    x = np.arange(len(plot_order))  # the site label locations
    width = 0.25  # the width of the bars

    max_min_y_vals = []
    for i, exceed_type in enumerate(exceed_type_list):
        crustal_probs = \
            get_probability_plot_data(site_PPE_dictionary=crustal_site_PPE_dictionary, exceed_type=exceed_type,
                                      threshold=threshold, site_list=plot_order)

        sz_probs = \
            get_probability_plot_data(site_PPE_dictionary=sz_site_PPE_dictionary, exceed_type=exceed_type,
                                      threshold=threshold, site_list=plot_order)
        max_min_y_vals.append(max(crustal_probs + sz_probs))
        max_min_y_vals.append(min(crustal_probs + sz_probs))

        # add bars to plot
        bars_crustal = axs[i].bar(x, crustal_probs, width, color=colors, edgecolor="0.2", linewidth=0.5, alpha=0.5)
        bars_sz = axs[i].bar(x + width, sz_probs, width, color=colors, edgecolor='black', linewidth=0.5)


        for bar in bars_crustal:
            bar_color = bar.get_facecolor()
            # add value label to each bar
            axs[i].text(x=(-0.12 + bar.get_x() + bar.get_width() / 2), y=bar.get_height() + 0.01,
                        s=f"{int(100 * round(bar.get_height(), 2))}%", horizontalalignment='center', color='k',
                        fontsize=6, fontweight='black', alpha=0.9)
            axs[i].text(x=(-0.12 + bar.get_x() + bar.get_width() / 2), y=bar.get_height() + 0.01,
                        s=f"{int(100 * round(bar.get_height(), 2))}%", horizontalalignment='center', color='white',
                        fontsize=6, fontweight='bold')
            axs[i].text(x=(-0.12 + bar.get_x() + bar.get_width() / 2), y=bar.get_height() + 0.01,
                        s=f"{int(100 * round(bar.get_height(), 2))}%", horizontalalignment='center', color=bar_color,
                        fontsize=6, fontweight='bold')
        for bar in bars_sz:
            bar_color = bar.get_facecolor()
            # add value label to each bar
            axs[i].text(x=(0.12 + bar.get_x() + bar.get_width() / 2), y=bar.get_height() + 0.01,
                        s=f"{int(100 * round(bar.get_height(), 2))}%", horizontalalignment='center', color=bar_color,
                        fontsize=6, fontweight='bold')
    for i in range(len(exceed_type_list)):
        axs[i].set_ylim(0.0, 0.4)
        axs[i].tick_params(axis='x', labelrotation=90, labelsize=6)
        axs[i].tick_params(axis='y', labelsize=8)
        axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # set tick labels to be every 0.2
        axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.2))
        axs[i].set_xticks(x + width/2, plot_order)

    # set indidual subplot stuff
    axs[0].set_ylabel("Probabilty of exceedance", fontsize=8)
    axs[1].tick_params(axis='y', labelleft=False)

    axs[0].set_title(f"{threshold} m uplift", fontsize=8)
    axs[1].set_title(f"{threshold} m subsidence", fontsize=8)

    # manually make legend with rectangles and text
    swatch_width, swatch_height = width, max(max_min_y_vals) * 0.08
    swatch_minx, swatch_miny = -1 * (len(plot_order) / 30), max(max_min_y_vals)
    crustal_swatches, sz_swatches = [], []
    offset = 0
    for color in colors:
        crustal_swatches.append(Rectangle((swatch_minx + offset, swatch_miny), swatch_width, swatch_height,
                                          facecolor=color, edgecolor='black',linewidth=0.5,))
        sz_swatches.append(Rectangle((swatch_minx + offset, swatch_miny - 2 * swatch_height), swatch_width,
                                     swatch_height, facecolor=color, edgecolor="0.2", linewidth=0.5, alpha=0.5))
        offset += swatch_width

    for i in range(len(crustal_swatches)):
        axs[0].add_patch(crustal_swatches[i])
        axs[0].add_patch(sz_swatches[i])

    axs[0].text(len(crustal_swatches) * swatch_width + swatch_width, swatch_miny, "subduction zone", fontsize=6)
    axs[0].text(len(crustal_swatches) * swatch_width + swatch_width, swatch_miny - 2 * swatch_height,
                "crustal faults", fontsize=6)

    fig.suptitle(f"100 yr exceedance probabilities", fontsize=10)
    fig.tight_layout()

    # make a directory for the figures if it doesn't already exist
    if not os.path.exists(f"../{results_directory}/compare_{crustal_file_suffix}_{sz_file_suffix}"):
        os.makedirs(f"../{results_directory}/compare_{crustal_file_suffix}_{sz_file_suffix}")

    fig.savefig(f"../{results_directory}/compare_{crustal_file_suffix}_{sz_file_suffix}"
                f"/probs_bar_chart_{crustal_file_suffix}_{sz_file_suffix}.png",
                dpi=300)


# def compare_disp_bar_chart(crustal_file_suffix, sz_file_suffix,
#                                    slip_taper, results_directory):
#     """ makes bar charts of the displacement value at the 10% and 2% probability of exceence thresholds for each site
#         crustal extension list: one is extension1 that designates the branch, the other defines the taper (uniform or tapered)
#         """
#     plot_order = ["Paraparaumu", "Porirua CBD north", "South Coast", "Wellington Airport", "Wellington CBD", "Petone",
#                   "Seaview", "Eastbourne", "Turakirae Head", "Lake Ferry", "Cape Palliser", "Flat Point"]
#
#     probability_list = [0.1, 0.02]
#
#     if slip_taper:
#         taper_extension = "_tapered"
#     else:
#         taper_extension = "_uniform"
#
#     with open(f"../{results_directory}/crustal{crustal_model_extension}"
#               f"/sites{crustal_file_suffix}"
#               f"/crustal_cumu_exceed_prob_sites{crustal_file_suffix}{taper_extension}.pkl", "rb") as fid:
#         crustal_site_PPE_dictionary = pkl.load(fid)
#     with open(f"../{results_directory}/sz_{sz_model_version}"
#               f"/sizes{sz_file_suffix}/sz_cumu_exceed_prob_sites{sz_file_suffix}"
#               f".pkl", "rb") as fid:
#         sz_site_PPE_dictionary = pkl.load(fid)
#
#     fig, axs = plt.subplots(1, 2, figsize=(10, 6))
#     x = np.arange(len(plot_order))  # the site label locations
#     width = 0.3  # the width of the bars
#
#     max_min_y_vals = []
#     for i, probability in enumerate(probability_list):
#         crustal_disps_up = \
#             get_exceedance_plot_data(site_PPE_dictionary=crustal_site_PPE_dictionary, exceed_type="up",
#                                     site_list=plot_order, probability=probability)
#         crustal_disps_down = \
#             get_exceedance_plot_data(site_PPE_dictionary=crustal_site_PPE_dictionary, exceed_type="down",
#                                      site_list=plot_order, probability=probability)
#
#         sz_disps_up = \
#             get_exceedance_plot_data(site_PPE_dictionary=sz_site_PPE_dictionary, exceed_type="up",
#                                      site_list=plot_order, probability=probability)
#         sz_disps_down = \
#             get_exceedance_plot_data(site_PPE_dictionary=sz_site_PPE_dictionary, exceed_type="down",
#                                      site_list=plot_order, probability=probability)
#
#         # max value in either crustal or sz up
#         max_min_y_vals.append(max(crustal_disps_up + sz_disps_up))
#         max_min_y_vals.append(min(crustal_disps_down + sz_disps_down))
#
#         color_up, color_down = (189/255, 0, 0), (15/255, 72/255, 186/255)
#         label_size1, label_size2 = 8, 6
#         label_offset = label_size2 / 150
#
#         # add bars to plot
#         bars_crustal_up = axs[i].bar(x, crustal_disps_up, width, color=color_up, linewidth=0.5, alpha=0.5)
#         bars_crustal_down = axs[i].bar(x, crustal_disps_down, width, color=color_down, linewidth=0.5, alpha=0.5)
#         bars_sz_up = axs[i].bar(x + width, sz_disps_up, width, color=color_up, linewidth=0.5)
#         bars_sz_down = axs[i].bar(x + width, sz_disps_down, width, color=color_down, linewidth=0.5)
#
#         # add value labels to bars
#         for bar in bars_crustal_up:
#             bar_color = bar.get_facecolor()
#             axs[i].text(x=(bar.get_x()) + 0.8 * width, y=bar.get_height() + label_offset,
#                         s=round(bar.get_height(), 1), ha='right', va='center', color=bar_color,
#                         fontsize=label_size2, fontweight='bold')
#         for bar in bars_sz_up:
#             bar_color = bar.get_facecolor()
#             axs[i].text(x=bar.get_x() + 0.2 * width, y=bar.get_height() + label_offset,
#                         s=round(bar.get_height(), 1), ha='left', va='center', color=bar_color,
#                         fontsize=label_size2, fontweight='bold')
#         for bar in bars_crustal_down:
#             bar_color = bar.get_facecolor()
#             axs[i].text(x=bar.get_x() + 0.8 * width, y=bar.get_height() - label_offset,
#                         s=round(bar.get_height(), 1), ha='right', va='center',
#                         color=bar_color, fontsize=label_size2, fontweight='bold')
#         for bar in bars_sz_down:
#             bar_color = bar.get_facecolor()
#             axs[i].text(x=bar.get_x() + 0.2 * width, y=bar.get_height() - label_offset,
#                         s=round(bar.get_height(), 1), ha='left', va='center',
#                         color=bar_color, fontsize=label_size2, fontweight='bold')
#
#     for i in range(len(axs)):
#         axs[i].axhline(y=0, color="k", linewidth=0.5)
#         axs[i].set_ylim(min(max_min_y_vals) - 0.2, max(max_min_y_vals) + 0.2)
#         axs[i].tick_params(axis='x', labelrotation=90, labelsize=label_size1)
#         axs[i].tick_params(axis='y', labelsize=label_size1)
#         axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#         axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.5))
#         axs[i].set_xticks(x + width/2, plot_order)
#
#     # set indidual subplot stuff
#     axs[0].set_ylabel("Displacement (m)", fontsize=8)
#     axs[1].tick_params(axis='y', labelleft=False)
#
#     axs[0].set_title(f"10% probability of exceedance", fontsize=label_size1)
#     axs[1].set_title(f"2% probability of exceedance", fontsize=label_size1)
#
#     # manually make legend with rectangles and text
#     swatch_width, swatch_height = width, max(max_min_y_vals) * 0.08
#     swatch_minx, swatch_miny = -1 * (len(plot_order) / 30), max(max_min_y_vals)
#     axs[0].add_patch(Rectangle((swatch_minx, swatch_miny), swatch_width, swatch_height,
#                                facecolor=color_up, edgecolor=None))
#     axs[0].add_patch(Rectangle((swatch_minx + swatch_width, swatch_miny), swatch_width, swatch_height,
#                                facecolor=color_down, edgecolor=None))
#     axs[0].add_patch(Rectangle((swatch_minx, swatch_miny - 2 * swatch_height), swatch_width, swatch_height,
#                                facecolor=color_up, edgecolor=None, alpha=0.5))
#     axs[0].add_patch(Rectangle((swatch_minx + swatch_width, swatch_miny - 2 * swatch_height), swatch_width, swatch_height,
#                                facecolor=color_down, edgecolor=None, alpha=0.5))
#
#     axs[0].text(swatch_minx + 3 * swatch_width, swatch_miny + swatch_height/2, "subduction zone",
#                 fontsize=label_size1, va='center')
#     axs[0].text(swatch_minx + 3 * swatch_width, swatch_miny - 3/2 * swatch_height, "crustal faults",
#                 fontsize=label_size1,
#                 va='center')
#
#     fig.suptitle(f"100 yr exceedance displacements", fontsize=label_size1 + 2)
#     fig.tight_layout()
#
#     # make a directory for the figures if it doesn't already exist
#     if not os.path.exists(f"../{results_directory}/compare_{crustal_file_suffix}_{sz_file_suffix}"):
#         os.makedirs(f"../{results_directory}/compare_{crustal_file_suffix}_{sz_file_suffix}")
#
#     fig.savefig(f"../{results_directory}/compare_{crustal_file_suffix}_{sz_file_suffix}"
#                 f"/disps_bar_chart_{crustal_file_suffix}_{sz_file_suffix}.png",
#                 dpi=300)

# def compare_branches(branch_extension_list, model_version, fault_type, results_directory):
#     """ only works for calculations done at named sites
#     must be two branches of the same fualt type (crustal """
#
#     probability_list = [0.1, 0.02]
#
#     branch_dict = {}
#     for i, extension in enumerate(branch_extension_list):
#         with open(f"../{results_directory}/{fault_type}{model_version}", "rb") as fid:
#             dict_i = pkl.load(fid)
#         branch_dict[extension] = dict_i
#
#     fig, axs = plt.subplots(1, 2, figsize=(10, 6))
#     x = np.arange(len(plot_order))  # the site label locations
#     width = 0.3  # the width of the bars
#     label_size1, label_size2 = 8, 6
#
#     max_min_y_vals = []
#     colors = plt.get_cmap('tab20')(np.linspace(0, 1, len(branch_dict)))
#
#     for j, branch in enumerate(branch_dict):
#         for i, probability in enumerate(probability_list):
#             disps_up = \
#                 get_exceedance_plot_data(site_PPE_dictionary=branch_dict[branch], exceed_type="up",
#                                          site_list=plot_order, probability=probability)
#             disps_down = \
#                 get_exceedance_plot_data(site_PPE_dictionary=branch_dict[branch], exceed_type="down",
#                                          site_list=plot_order, probability=probability)
#
#             # max value in either crustal or sz up
#             max_min_y_vals.append(max(disps_up))
#             max_min_y_vals.append(min(disps_down))
#
#             # add points to plot
#             scatter_up = axs[i].scatter(x, disps_up, color=colors[j])
#             scatter_down = axs[i].scatter(x, disps_down, color=colors[j])
#
#     for i in range(len(axs)):
#         axs[i].axhline(y=0, color="k", linewidth=0.5)
#         axs[i].set_ylim(min(max_min_y_vals) - 0.2, max(max_min_y_vals) + 0.2)
#         axs[i].tick_params(axis='x', labelrotation=90, labelsize=label_size1)
#         axs[i].tick_params(axis='y', labelsize=label_size1)
#         axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#         axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.5))
#         axs[i].set_xticks(x, plot_order)
#
#     # set indidual subplot stuff
#     axs[0].set_ylabel("Displacement (m)", fontsize=8)
#     axs[1].tick_params(axis='y', labelleft=False)
#
#     axs[0].set_title(f"10% probability of exceedance", fontsize=label_size1)
#     axs[1].set_title(f"2% probability of exceedance", fontsize=label_size1)
#
#     #fig.legend()
#     legend_start_x, legend_start_y = -1 * (len(plot_order) / 30), max(max_min_y_vals)
#     multiplier = 0
#     for i, extension in enumerate(extension_list):
#         if i == round(len(extension_list) / 2):
#             multiplier = 0
#             legend_start_x += 1.4
#         axs[0].scatter(legend_start_x, legend_start_y - multiplier, color=colors[i], label=extension)
#         axs[0].text(legend_start_x + 0.2, legend_start_y - multiplier, extension_list[i], fontsize=label_size2, va='center')
#         multiplier += 0.1
#
#     fig.suptitle(f"100 yr exceedance displacements, {fault_type}", fontsize=label_size1 + 2)
#     fig.tight_layout()
#
#     # make a directory for the figures if it doesn't already exist
#     if not os.path.exists(f"../{results_directory}/compare_branches"):
#         os.makedirs(f"../{results_directory}/compare_branches")
#
#     fig.savefig(f"../{results_directory}/{fault_type}_disps_branches.png", dpi=300)


def compare_faultmodel_prob_plot(PPE_paths, plot_name, title, names, outfile_directory,
                                 file_type_list=["png"],
                                      threshold=0.2):
    """ """

    exceed_type_list = ["up", "down"]


    PPE_dicts = []
    for PPE_path in PPE_paths:
        with open(f"{PPE_path}",  "rb") as fid:
            PPE_dict = pkl.load(fid)
        PPE_dicts.append(PPE_dict)

    # set up custom color scheme
    colors = [(255/255, 190/255, 106/255), (64/255, 176/255, 166/255)]

    # set up figure and subplots
    fig, axs = plt.subplots(1, 2, figsize=(7, 5))
    x = np.arange(len(plot_order))  # the site label locations

    for p, PPE_dict in enumerate(PPE_dicts):
        all_max_errs_y = []
        for i, exceed_type in enumerate(exceed_type_list):
            mean_prob, errs_plus, errs_minus = \
                get_mean_prob_barchart_data(site_PPE_dictionary=PPE_dict, exceed_type=exceed_type,
                                            threshold=threshold, site_list=plot_order)

            # add point and error bars to plot
            axs[i].errorbar(x, mean_prob, yerr=[errs_plus, errs_minus], fmt='none', ecolor=colors[p], capsize=5,
                            linewidth=1, markeredgewidth=1)
            axs[i].scatter(x, mean_prob, s=24, color=colors[p], zorder=3, edgecolors='k', linewidths=0.5)

            max_errs_y = max([mean_prob[j] + errs_plus[j] for j in range(len(mean_prob))])
            all_max_errs_y.append(max_errs_y)

            ## add labels to each point
            # for j, value in enumerate(mean_prob):
            #     point_color = err_colors[p]
            #
            #     # add value label to each bar
            #     label_offset = 0.05 * max_errs_y
            #     axs[i].text(x=x[j] + label_offset, y=mean_prob[j] + label_offset,
            #                 s=f"{int(100 * round(mean_prob[j], 2))}%",
            #                 horizontalalignment='left', color=point_color,
            #                 fontsize=6, fontweight='demibold')

            if max(all_max_errs_y) < 0.25:
                axs[i].set_ylim(0.0, 0.25)
            else:
                axs[i].set_ylim(0.0, max(all_max_errs_y))
            axs[i].tick_params(axis='x', labelrotation=90, labelsize=6)
            axs[i].tick_params(axis='y', labelsize=8)
            axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            # set tick labels to be every 0.2
            axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.1))
            axs[i].set_xticks(x, plot_order)

    # set indidual subplot stuff
    fontsize = 8
    # legend1 = axs[0].legend(loc="upper left", title="Classes")
    # axs[0].add_artist(legend1)
    axs[0].legend(names, loc="upper left", title="Fault models", title_fontsize=fontsize, fontsize=fontsize)

    axs[0].set_ylabel("Probabilty", fontsize=8)
    axs[1].tick_params(axis='y', labelleft=False)

    axs[0].set_title(f"Probability of exceeding {threshold} m uplift", fontsize=fontsize)
    axs[1].set_title(f"Probability of exceeding {threshold} m subsidence", fontsize=fontsize)

    fig.suptitle(f"Compare {title} \n100 yrs")
    fig.tight_layout()

    if not os.path.exists(outfile_directory):
        os.makedirs(outfile_directory)

    for file_type in file_type_list:
        fig.savefig(f"{outfile_directory}/compare_probs_{plot_name}.{file_type}", dpi=300)



# for i in range(len(crustal_file_suffix_list)):
#     crustal_file_suffix = crustal_file_suffix_list[i]
#     sz_file_suffix = sz_file_suffix_list[i]
#
#     make_compare_disp_10_2_bar_chart(crustal_file_suffix=crustal_file_suffix, sz_file_suffix=sz_file_suffix,
#                                    slip_taper=slip_taper, results_directory=results_directory)
#     make_combo_prob_bar_chart(crustal_file_suffix=crustal_file_suffix, sz_file_suffix=sz_file_suffix,
#                               slip_taper=slip_taper, results_directory=results_directory,
#                               threshold=0.2)

# pkl_list = []
# for extension in sz_extension1_list:
#     pkl_list.append(f"../{results_directory}/{extension}/{fault_type}_cumu_exceed_prob_{extension}"
#                             f"_{extension3}.pkl")
#
# compare_branches(extension_list=sz_extension1_list, pkl_list=pkl_list, fault_type=fault_type,
#                         results_directory=results_directory)


compare_faultmodel_prob_plot(PPE_paths=mean_PPE_paths, plot_name=plot_name,
                             outfile_directory=outfile_directory, title=title, names=names, file_type_list=["png"],
                             threshold=0.2)