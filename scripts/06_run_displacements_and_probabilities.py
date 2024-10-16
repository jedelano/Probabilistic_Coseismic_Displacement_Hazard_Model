import pickle as pkl
import random
import os
import matplotlib
from pcdhm.shared_helper import get_rupture_disp_dict
from pcdhm.rupture_scenario_plotting import vertical_disp_figure
from pcdhm.probabalistic_displacement import get_site_disp_dict, get_cumu_PPE, plot_branch_hazard_curve, \
     make_10_2_disp_plot, make_branch_prob_plot

##### USER INPUTS   #####
# must run crustal and subduction lists/loops separately
results_directory = "results_r1"

slip_taper = False                    # True or False, only matters if crustal otherwise it defaults to false later.
fault_type = "sz"                  # "crustal or "sz"

# How many branches do you want to run?
# True: picks the most central branch (i.e. geologic, time independent, mid b and N). Saves time, useful for
# generating displacement maps for scenarios (e.g., gf_type=grid) or troubleshooting.
# False: runs all branches in a fault model. All branches are needed for '07_run_aggregate_weighted_branches.py'.
# False is most useful for gf_name = "sites"
single_branch = True

# True: Use the predefined subset of rupture IDs
# False: Make a random sample of rupture IDs + the predefined subset
specific_rupture_ids = False

# can only run one type of GF and fault geometry at a time
gf_name = "sites"                               # "sites" or "grid" or "coastal"
crustal_mesh_version = "_CFM"           # e.g.,  "_Model1", "_Model2", or "_CFM"
sz_mesh_version = "_multi50"                    # must match suffix in the subduction directory
out_files_directory = "mesh_gf_outfiles_r1"     # used for grabbing geojsons

# Can run more than one type of deformation model at a time (only matters for crustal)
deformation_model = "geologic and geodetic"          # "geologic" or "geodetic" or "geologic and geodetic"

# can select one or both "time dependent" and "independent branches" (only matters for crustal)
time_dependent = True       # True or False
time_independent = True     # True or False

# True: calculates displacements at the greens function type/locations for all ruptures
# False: skips calculating displacements(assumed it's already calculated and saved all_rupture_disps pickle file)
calculate_displacements = False

# True: calculates PPE dictionaries for the "sites" only gf_name; saves to "cumu_exceed_prob_{extension1}.pkl"
# False: skips calculating PPE dictionaries; uses saved file for plotting
calculate_probabilities = False

# making figures
file_type_list = ["png", "pdf"]               # e.g. ["png", "pdf"] or ["png"]
skipped_sites = ["Porirua CBD south"]       # ["site1", "site2", "site3"]
# Scenario displacement maps
# True: makes displacement map figures of a sample of rupture scenarios. Mostly useful for grid gfs.
# False: skips making displacement map figures
make_scenario_displacement_maps = False

# True: uses saved dictionaries to make probability and displacement figures for each branch
# False: skips making probability and displacement figures for each branch
make_branch_probability_figures = True

################ script ###################
# this makes figure fonts editable in Illustrator
matplotlib.rcParams['pdf.fonttype'] = 42

# This is quite possibly the dumbest part of my coding, but it's going to be annoying to fix.
# Set up which branches you want to calculate displacements and probabilities for
# better way is to set up a csv or something and read in the paths and names
file_suffix_list = []
NSHM_directory_list = []
if fault_type == "crustal":
    mesh_version = crustal_mesh_version
    if not single_branch:
        if time_independent:
            if "geologic" in deformation_model:
                file_suffix_list_i = ["_c_MDA2", "_c_MDEz", "_c_MDE1"]
                NSHM_directory_list_i = ["crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDA2",
                                         "crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDEz",
                                         "crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDE1"]
                file_suffix_list.extend(file_suffix_list_i)
                NSHM_directory_list.extend(NSHM_directory_list_i)
            if "geodetic" in deformation_model:
                file_suffix_list_i = ["_c_MDE2", "_c_MDE5", "_c_MDI0"]
                NSHM_directory_list_i = ["crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDE2",
                                         "crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDE5",
                                         "crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDI0"]
                file_suffix_list.extend(file_suffix_list_i)
                NSHM_directory_list.extend(NSHM_directory_list_i)

        if time_dependent:
            if "geologic" in deformation_model:
                file_suffix_list_i = ["_c_NjE5", "_c_NjI2", "_c_NjI3"]
                NSHM_directory_list_i = ["crustal_solutions/NZSHM22_TimeDependentInversionSolution-QXV0b21hdGlvblRhc2s6MTExNjE5",
                                         "crustal_solutions/NZSHM22_TimeDependentInversionSolution-QXV0b21hdGlvblRhc2s6MTExNjI2",
                                         "crustal_solutions/NZSHM22_TimeDependentInversionSolution-QXV0b21hdGlvblRhc2s6MTExNjI3"]
                file_suffix_list.extend(file_suffix_list_i)
                NSHM_directory_list.extend(NSHM_directory_list_i)
            if "geodetic" in deformation_model:
                file_suffix_list_i = ["_c_NjI5", "_c_NjMy", "_c_NjM3"]
                NSHM_directory_list_i = ["crustal_solutions/NZSHM22_TimeDependentInversionSolution-QXV0b21hdGlvblRhc2s6MTExNjI5",
                                         "crustal_solutions/NZSHM22_TimeDependentInversionSolution-QXV0b21hdGlvblRhc2s6MTExNjMy",
                                         "crustal_solutions/NZSHM22_TimeDependentInversionSolution-QXV0b21hdGlvblRhc2s6MTExNjM3"]
                file_suffix_list.extend(file_suffix_list_i)
                NSHM_directory_list.extend(NSHM_directory_list_i)

    if single_branch:
        file_suffix_list_i = ["_c_MDEz"]
        NSHM_directory_list_i = ["crustal_solutions/NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDEz"]
        file_suffix_list.extend(file_suffix_list_i)
        NSHM_directory_list.extend(NSHM_directory_list_i)

elif fault_type == "sz":
    mesh_version = sz_mesh_version
    slip_taper = False
    if not single_branch:
        file_suffix_list_i = ["_sz_MzI2", "_sz_MzMx",  "_sz_MzMy"]
        NSHM_directory_list_i = ["sz_solutions/NZSHM22_AveragedInversionSolution-QXV0b21hdGlvblRhc2s6MTA3MzI2",
                                 "sz_solutions/NZSHM22_AveragedInversionSolution-QXV0b21hdGlvblRhc2s6MTA3MzMx",
                                 "sz_solutions/NZSHM22_AveragedInversionSolution-QXV0b21hdGlvblRhc2s6MTA3MzMy"]
        file_suffix_list.extend(file_suffix_list_i)
        NSHM_directory_list.extend(NSHM_directory_list_i)
    if single_branch:
        file_suffix_list_i = ["_sz_MzMx"]
        NSHM_directory_list_i = ["sz_solutions/NZSHM22_AveragedInversionSolution-QXV0b21hdGlvblRhc2s6MTA3MzMx"]
        file_suffix_list.extend(file_suffix_list_i)
        NSHM_directory_list.extend(NSHM_directory_list_i)

# raise an error if the number of file suffixes and NSHM directories are not equal
if len(file_suffix_list) != len(NSHM_directory_list):
    raise ValueError("Number of file suffixes and NSHM directories must be equal")
# raise an error if there are duplicate file suffixes. Ask me how I know.
if len(file_suffix_list) != len(set(file_suffix_list)):
    raise ValueError("Duplicate file suffixes found")
# raise an error if there are duplicate NSHM directories.
if len(NSHM_directory_list) != len(set(NSHM_directory_list)):
    raise ValueError("Duplicate NSHM directories found")


naming_extension_list = [gf_name + suffix for suffix in file_suffix_list]
model_version_results_directory = f"{results_directory}/{fault_type}{mesh_version}"
grid = True if gf_name == "grid" else False


if not os.path.exists(f"../{results_directory}"):
    os.mkdir(f"../{results_directory}")
if not os.path.exists(f"../{model_version_results_directory}"):
    os.mkdir(f"../{model_version_results_directory}")

crustal_outfiles_path = f"../{out_files_directory}/crustal{crustal_mesh_version}"
sz_outfiles_path = f"../{out_files_directory}/sz{sz_mesh_version}"

if calculate_displacements is True:
    # Calculate displacements and make displacement dictionary once per branch. Save to pickle file in branch directory.
    # file saved as "all_rupture_disps_{extension1}.pkl"
    for i in range(len(naming_extension_list)):
        print (f"Calculating displacements for branch {i} in {len(naming_extension_list)}")
        get_rupture_disp_dict(NSHM_directory=NSHM_directory_list[i], naming_extension=naming_extension_list[i],
                              slip_taper=slip_taper, fault_type=fault_type, gf_name=gf_name,
                              results_version_directory=model_version_results_directory,
                              crustal_outfiles_path=crustal_outfiles_path, sz_outfiles_path=sz_outfiles_path)

### make vertical displacement figures (random sample of ~10 ruptures per branch)
## will make a colored gridded map for "grid" and a map with colored points for "sites"
if make_scenario_displacement_maps is True:

    taper_extension = "_tapered" if slip_taper is True else "_uniform"

    for i in range(len(naming_extension_list)):
        with open(f"../{model_version_results_directory}/{naming_extension_list[i]}/all_rupture_disps_{naming_extension_list[i]}"
                  f"{taper_extension}.pkl", "rb") as fid:
            all_ruptures_disp_dict = pkl.load(fid)
        rupture_id = list(all_ruptures_disp_dict.keys())

        if specific_rupture_ids:
            target_rupture_ids = []
        else:
            if fault_type == "crustal": num_ruptures = 5
            if fault_type == "sz": num_ruptures = 10
            target_rupture_ids = random.sample(rupture_id, num_ruptures)

        # tack on a few rupture scenarios that we want to see results from. Only include if they are in the rupture list
        if fault_type == "sz":
            ids = [948, 3043, 4365, 4977]
            for id in ids:
                if id in rupture_id: target_rupture_ids.append(id)
        elif fault_type == "crustal":
            ids = [20890, 96084, 97010, 166970, 305270, 368024, 375389, 401491]
            for id in ids:
                if id in rupture_id: target_rupture_ids.append(id)

        # make figures
        print(f"\n*~ Making displacement figures for {naming_extension_list[i]} ~*")
        vertical_disp_figure(NSHM_directory=NSHM_directory_list[i], all_ruptures_disp_dict=all_ruptures_disp_dict,
                             target_rupture_ids=target_rupture_ids, naming_extension=naming_extension_list[i],
                             extent="Wellington", slip_taper=slip_taper, grid=grid, fault_type=fault_type,
                             results_version_directory=model_version_results_directory,
                             crustal_outfiles_path=crustal_outfiles_path,
                             sz_outfiles_path=sz_outfiles_path,
                             mesh_version=mesh_version,
                             file_type_list=file_type_list)

if calculate_probabilities is True and gf_name == "sites":
    ## calculate rupture branch probabilities and make plots
    for i in range(len(naming_extension_list)):
        print(f"*~ Calculating probabilities for {naming_extension_list[i]} ~*\n")
        # ### step 1: get site displacement dictionary
        branch_site_disp_dict = get_site_disp_dict(naming_extension_list[i], slip_taper=slip_taper,
                                                   model_version_results_directory=model_version_results_directory)
        #
        # ### step 2: get exceedance probability dictionary
        #    # file saved as "cumu_exceed_prob_{extension1}.pkl"
        get_cumu_PPE(extension1=naming_extension_list[i], branch_site_disp_dict=branch_site_disp_dict,
                     model_version_results_directory=model_version_results_directory, slip_taper=slip_taper,
                     time_interval=100, n_samples=1000000)
        #
# ## step 3 (optional): plot hazard curves
if make_branch_probability_figures is True and gf_name == "sites":
    for i in range(len(naming_extension_list)):
        print(f"*~ Making probability figures for {naming_extension_list[i]} ~*\n")
        plot_branch_hazard_curve(extension1=naming_extension_list[i],
                                 model_version_results_directory=model_version_results_directory,
                                 slip_taper=slip_taper, file_type_list=file_type_list, skipped_sites=skipped_sites)

        # step 4 (optional): plot hazard maps
        # plot_cumu_disp_hazard_map(extension1=extension1_list[i], slip_taper=slip_taper, grid=grid, fault_type=fault_type,
        #                           model_version_results_directory=model_version_results_directory,
        #                           crustal_directory=crustal_directory,
        #                           sz_directory=sz_directory, model_version=model_version)

        ## step 5: probability and displacement plots
        # graph with probability (%) at the 0.2 m uplift and subsidence thresholds (can change threshold)
        make_branch_prob_plot(naming_extension_list[i], slip_taper=slip_taper, threshold=0.2, skipped_sites=skipped_sites,
                              model_version_results_directory=model_version_results_directory,
                              model_version=mesh_version)
        # bar chart with displacements at 10% and 2% exceedance probabilities
        make_10_2_disp_plot(extension1=naming_extension_list[i], slip_taper=slip_taper,
                            model_version_results_directory=model_version_results_directory,
                            skipped_sites=skipped_sites)
