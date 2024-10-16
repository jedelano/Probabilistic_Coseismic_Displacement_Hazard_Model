from pcdhm.probabalistic_displacement import \
    make_sz_crustal_paired_PPE_dict, \
     make_fault_model_PPE_dict, \
     get_weighted_mean_PPE_dict
from pcdhm.weighted_mean_plotting import \
    plot_weighted_mean_haz_curves, \
    plot_weighted_mean_haz_curves_colorful
from pcdhm.shared_helper import \
    make_branch_weight_dict
import pickle as pkl

"""this script takes the displcements and probabilities from each fault model branch and aggregates them into a
weighed mean. It then plots the hazard curves/probability figures.

you can run this script for a single fault model or a paired crustal and subduction model. The paired model takes
the individual branch results and combines them, effectively creating n = X crustal * Y subduction new branches."""

#### USER INPUTS   #####
# set up file directories
results_directory = 'results_r1'
# Optional; something to tack on to the end of the aggregated PPE dictionary so you don't overwrite files
outfile_extension = ''

# choose fault type and model versions
slip_taper = False                           # True or False, only matters if crustal. Defaults to False for sz.

# Run PPE or figures for a single fault model or a paired crustal/subduction model?
paired_crustal_sz = False                  # True: combined sz and crustal; False: uses [single] fault type
single_fault_type = 'sz'                       # "crustal or "sz"; only matters for single fault model

# "_Model1", "_Model2", or "_CFM" for crustal; must match results subdirectory suffix
crustal_model_version = "_CFM"
sz_model_version = "_multi50_gentlerdip"

# Do you want to calculate PPEs (probabilities) and weighted mean PPE for the fault model?
# True: Only has to be done once because it is saved in a pickle file
# False: uses saved pickle file; saves time
calculate_fault_model_PPE = True               # calculate poissonian probabilities of exceedance

# choose which figures to make
figure_file_type_list = ["png", "pdf"]             # file types for figures
probability_plot = False                         # plots the prob at the 0.2 m uplift and subsidence thresholds
displacement_chart = False                       # plots disp at the 10% and 2% probability of exceedance
make_hazcurves = False                               # plots all branches and weight mean haz curves
make_colorful_hazcurves = True                     # plots branches colored by unique_id_keyphrase_list
make_map = False                                     # displacements on the left, map on the right

# sites to omit from plots
skipped_sites = ["Porirua CBD south"]       # ["site1", "site2", "site3"]

# identify which branches are colored for "colorful haz curves"
# more than one can be true, it just makes multiple different figures.
# only matters if make_colorful_hazcurves is True
show_b_n_variation = True                   # magnitude-frequency distribution
show_s_value_variation = True               # non-stationary moment rate scaling
show_def_model_variation = True            # deformation model
show_time_dependence = True                # time dependent or time independent

#choose which version of the figures to make
exceed_type_list = ["down"] # e.g. ["up", "down", "total_abs"]

########## script ############################################
# get branch weights from the saved Excel spreadsheet
branch_weight_file_path = f"../data/branch_weight_data.xlsx"
crustal_sheet_name = "crustal_weights_4_2"
sz_sheet_name = "sz_weights_4_0"

#######
# only works for sites greens functions due to computational intensity
gf_name = "sites"

if slip_taper is True:
    taper_extension = "_tapered"
else:
    taper_extension = "_uniform"

# plot_order = get_plot_order_list(skipped_sites=skipped_sites)

# get branch weights from the saved Excel spreadsheet
crustal_branch_weight_dict = make_branch_weight_dict(branch_weight_file_path=branch_weight_file_path,
                                                     sheet_name=crustal_sheet_name)
sz_branch_weight_dict = make_branch_weight_dict(branch_weight_file_path=branch_weight_file_path,
                                                sheet_name=sz_sheet_name)

# these directories should already be made from calculating displacements in a previous script
if paired_crustal_sz:
    crustal_model_version_results_directory = f"{results_directory}/crustal{crustal_model_version}"
    sz_model_version_results_directory = f"{results_directory}/sz{sz_model_version}"

# designate which branch weight dictionary to use based on the fault type
if not paired_crustal_sz and single_fault_type == "crustal":
    fault_model_branch_weight_dict = crustal_branch_weight_dict
    fault_model_version = crustal_model_version
if not paired_crustal_sz and single_fault_type == "sz":
    fault_model_branch_weight_dict = sz_branch_weight_dict
    fault_model_version = sz_model_version


### make a dictionary of all the branch probabilities, oranized by site within each branch
# option to skip this step if you've already run it once and saved to a pickle file
if not paired_crustal_sz:
    n_samples = 1000000

    model_version_results_directory = f"{results_directory}/{single_fault_type}{fault_model_version}"
    fault_model_PPE_filepath = f"../{model_version_results_directory}/allbranch_PPE_dict{outfile_extension}{taper_extension}.pkl"

    # Calculate the probabilities for each branch and save to a pickle file.
    if calculate_fault_model_PPE:
        make_fault_model_PPE_dict(
            branch_weight_dict=fault_model_branch_weight_dict,
            model_version_results_directory=model_version_results_directory, n_samples=n_samples,
            slip_taper=slip_taper, outfile_extension=outfile_extension)

        with open(fault_model_PPE_filepath, 'rb') as f:
            PPE_dict = pkl.load(f)

        weighted_mean_PPE_dict = get_weighted_mean_PPE_dict(fault_model_PPE_dict=PPE_dict,
                                                            out_directory=model_version_results_directory,
                                                            outfile_extension=outfile_extension, slip_taper=slip_taper)

    else:
        with open(fault_model_PPE_filepath, 'rb') as f:
            PPE_dict = pkl.load(f)

##### paired crustal and sz PPE
if paired_crustal_sz:
    n_samples = 100000
    model_version_results_directory = f"{results_directory}/paired_c{crustal_model_version}_sz{sz_model_version}"

    paired_PPE_pickle_name = f"sz_crustal_paired_allbranch_PPE_dict{outfile_extension}{taper_extension}.pkl"
    paired_PPE_filepath = f"../{model_version_results_directory}/{paired_PPE_pickle_name}"

    #### calculate the probabilities for each branch combination (n Crustal * m SZ) and saves to a pickle file
    # skip this part if you've already run it once and saved to a pickle file
    if calculate_fault_model_PPE:
        # the script automatically saves the paired PPE dictionary to a pickle file
        make_sz_crustal_paired_PPE_dict(
            crustal_branch_weight_dict=crustal_branch_weight_dict, sz_branch_weight_dict=sz_branch_weight_dict,
            crustal_model_version_results_directory=crustal_model_version_results_directory,
            sz_model_version_results_directory=sz_model_version_results_directory,
            slip_taper=slip_taper,
            n_samples=n_samples,
            out_directory=model_version_results_directory,
            paired_PPE_pickle_name=paired_PPE_pickle_name)

        with open(paired_PPE_filepath, 'rb') as f:
            PPE_dict = pkl.load(f)
        weighted_mean_PPE_dict = get_weighted_mean_PPE_dict(fault_model_PPE_dict=PPE_dict,
                                                            out_directory=model_version_results_directory,
                                                            outfile_extension=outfile_extension, slip_taper=slip_taper)

    else:
        with open(paired_PPE_filepath, 'rb') as f:
            PPE_dict = pkl.load(f)

# open the saved weighted mean PPE dictionary
weighted_mean_PPE_filepath = f"../{model_version_results_directory}/weighted_mean_PPE_dict{outfile_extension}" \
                             f"{taper_extension}.pkl"
with open(weighted_mean_PPE_filepath, 'rb') as f:
    weighted_mean_PPE_dict = pkl.load(f)

# plot hazard curves and save to file
if paired_crustal_sz:
    model_version_title = f"paired crustal{crustal_model_version} and sz{sz_model_version}"
else:
    model_version_title = f"{single_fault_type}{fault_model_version}"

if make_hazcurves:
    plot_weighted_mean_haz_curves(
        PPE_dictionary=PPE_dict,
        weighted_mean_PPE_dictionary=weighted_mean_PPE_dict,
        model_version_title=model_version_title,
        exceed_type_list=["up", "down", "total_abs"],
        out_directory=model_version_results_directory,
        file_type_list=figure_file_type_list,
        slip_taper=slip_taper,
        skipped_sites=skipped_sites)

if make_colorful_hazcurves:
    if show_b_n_variation is True:
        if single_fault_type == "sz":
            unique_id_keyphrase_list = ["N165", "N279"]
        elif single_fault_type == "crustal":
            unique_id_keyphrase_list = ["N27", "N46"]
        plot_weighted_mean_haz_curves_colorful(weighted_mean_PPE_dictionary=weighted_mean_PPE_dict,
                                               PPE_dictionary=PPE_dict,
                                               exceed_type_list=exceed_type_list,
                                               model_version_title=model_version_title,
                                               out_directory=model_version_results_directory,
                                               file_type_list=figure_file_type_list,
                                               slip_taper=slip_taper, file_name=f"bn_val_branches"
                                                                                f"{fault_model_version}",
                                               string_list=unique_id_keyphrase_list, skipped_sites=skipped_sites)

    if show_s_value_variation is True:
        if single_fault_type == "sz":
            unique_id_keyphrase_list = ["S042", "S158"]
        elif single_fault_type == "crustal":
            unique_id_keyphrase_list = ["S066", "S141"]
        plot_weighted_mean_haz_curves_colorful(weighted_mean_PPE_dictionary=weighted_mean_PPE_dict,
                                               PPE_dictionary=PPE_dict,
                                               exceed_type_list=exceed_type_list,
                                               model_version_title=model_version_title,
                                               out_directory=model_version_results_directory,
                                               file_type_list=figure_file_type_list,
                                               slip_taper=slip_taper, file_name=f"s_val_branches{fault_model_version}",
                                               string_list=unique_id_keyphrase_list, skipped_sites=skipped_sites)

    if show_def_model_variation is True:
        if single_fault_type =="crustal":
            unique_id_keyphrase_list = ["geologic", "geodetic"]
        plot_weighted_mean_haz_curves_colorful(weighted_mean_PPE_dictionary=weighted_mean_PPE_dict,
                                               PPE_dictionary=PPE_dict,
                                               exceed_type_list=exceed_type_list,
                                               model_version_title=model_version_title,
                                               out_directory=model_version_results_directory,
                                               file_type_list=figure_file_type_list,
                                               slip_taper=slip_taper,
                                               file_name=f"def_model_branches{fault_model_version}",
                                               string_list=unique_id_keyphrase_list,
                                               skipped_sites=skipped_sites)
    if show_time_dependence is True:
        if single_fault_type =="crustal":
            unique_id_keyphrase_list = ["_TD_", "_TI_"]
        plot_weighted_mean_haz_curves_colorful(weighted_mean_PPE_dictionary=weighted_mean_PPE_dict,
                                               PPE_dictionary=PPE_dict,
                                               exceed_type_list=exceed_type_list,
                                               model_version_title=model_version_title,
                                               out_directory=model_version_results_directory,
                                               file_type_list=figure_file_type_list,
                                               slip_taper=slip_taper,
                                               file_name=f"time_dependence_branches{fault_model_version}",
                                               string_list=unique_id_keyphrase_list,
                                               skipped_sites=skipped_sites)

######### This needs a lil TLC before it works correctly in this script
# if make_map:
#     for exceed_type in exceed_type_list:
#         for disp in [0.1, 0.02]:
#             map_and_plot_probabilities(PPE_path=weighted_mean_PPE_filepath,
#                                        plot_name=f"map_probs{fault_model_version}",
#                                        exceed_type=exceed_type,
#                                        title=f"probability of exceeding {disp} {exceed_type}\n{fault_model_version}",
#                                        outfile_directory=model_version_results_directory,
#                                        plot_order=plot_order,
#                                        labels_on=True,
#                                        file_type_list=file_type_list,
#                                        threshold=disp,
#                                        colorbar_max=0.3)

