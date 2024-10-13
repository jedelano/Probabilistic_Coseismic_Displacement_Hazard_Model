import pandas as pd
import os
import itertools
import pickle as pkl
import matplotlib.ticker as mticker
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
from pcdhm.shared_helper import  make_qualitative_colormap,  get_probability_color, get_plot_order_list
from matplotlib.patches import Rectangle

# This makes fonts editable in exported PDFs
matplotlib.rcParams['pdf.fonttype'] = 42


def get_site_disp_dict(extension1, slip_taper, model_version_results_directory):
    """
        inputs: uses extension naming scheme to load displacement dictionary created with the
        get_rupture_disp_dict function. State slip taper (True or False).

        functions: reshapes the dictionary to be organized by site name (key = site name).

        outputs: a dictionary where each key is a location/site name. contains displacements (length = number of
        rupture ids), annual rate (length = number of rupture ids), site name list (should be same length as green's
        function), and site coordinates (same length as site name list)

        CAVEATS/choices:
        - a little clunky because most of the dictionary columns are repeated across all keys.
        """

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    # load saved displacement data
    # disp dictionary has keys for each rupture id and displacement data for each site ( disp length = # of sites)
    # for grids, the disp_dict is still in grid shape instead of list form.
    with open(f"../{model_version_results_directory}/{extension1}/all_rupture_disps_{extension1}{taper_extension}.pkl",
              "rb") as fid:
        rupture_disp_dictionary = pkl.load(fid)

    ###### reshape displacement data to be grouped by site location.
    first_key = list(rupture_disp_dictionary.keys())[0]
    site_names = rupture_disp_dictionary[first_key]["site_name_list"]
    site_coords = rupture_disp_dictionary[first_key]["site_coords"]

    # makes a list of lists. each item is a rupture scenario that contains a list of displacements at each site.
    # could probably simplify this because it's the same shape as the dictionary, but I'm not sure how to do that yet
    # and this is easier to understand I think.
    disps_by_scenario = []
    annual_rates_by_scenario = []
    for rupture_id in rupture_disp_dictionary.keys():
        disps = rupture_disp_dictionary[rupture_id]["v_disps_m"]
        disps_by_scenario.append(disps)
        annual_rate = rupture_disp_dictionary[rupture_id]["annual_rate"]
        annual_rates_by_scenario.append(annual_rate)

    # list of lists. each item is a site location that contains displacements from each scenario (disp list length =
    # number of rupture scenarios)
    disps_by_location = []
    annual_rates_by_location = []
    for site_num in range(len(site_names)):
        site_disp = [scenario[site_num] for scenario in disps_by_scenario]
        disps_by_location.append(site_disp)
        annual_rates_by_location.append(annual_rates_by_scenario)

    # make dictionary of displacements and other data. key is the site name.
    site_disp_dictionary = {}
    for i, site in enumerate(site_names):
        site_disp_dictionary[site] = {"disps": disps_by_location[i], "rates": annual_rates_by_location[i],
                                           "site_coords": site_coords[i]}

    if slip_taper is True:
        extension3 = "_tapered"
    else:
        extension3 = "_uniform"

    # with open(f"../{results_version_directory}/{extension1}/site_disp_dict_{extension1}{extension3}.pkl",
    #           "wb") as f:
    #     pkl.dump(site_disp_dictionary, f)
    return site_disp_dictionary


def get_cumu_PPE(slip_taper, model_version_results_directory, branch_site_disp_dict,  n_samples,
                 extension1, time_interval=100, sd=0.4):
    """calculates the poissonian probability of exceedance for each site for each displacement threshold value
    Must first run get_site_disp_dict to get the dictionary of displacements and rates

    Time_interval is in years

    outputs: pickle file with probability dictionary

    """

    # use random number generator to initial monte carlo sampling
    rng = np.random.default_rng()

    # Load the displacement/rate data for all sites
    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    ## loop through each site and generate a bunch of 100 yr interval scenarios
    site_PPE_dict = {}
    for i, site_of_interest in enumerate(branch_site_disp_dict.keys()):
        # if i == 0:
        #     if branch_key not in ["nan", ""]:
        #         print(f"calculating {branch_key} PPE for site {i} of {len(branch_site_disp_dict.keys())}")
        #     if extension1 not in ["nan", ""]:
        #         print(f"calculating {extension1} PPE for site {i} of {len(branch_site_disp_dict.keys())}")
        # print(f"calculating {branch_key} PPE ({i} of {len(branch_site_disp_dict.keys())} branches)")
        site_dict_i = branch_site_disp_dict[site_of_interest]

        ## Set up params for sampling
        investigation_time = time_interval

        if "scaled_rates" not in site_dict_i.keys():
            # if no scaled_rate column, assumes scaling of 1 (equal to "rates")
            scaled_rates = site_dict_i["rates"]
        else:
            scaled_rates = site_dict_i["scaled_rates"]

        # average number of events per time interval (effectively R*T from Ned's guide)
        lambdas = investigation_time * np.array(scaled_rates)

        # Generate n_samples of possible earthquake ruptures for random 100 year intervals (or investigation time interval)
        # returns boolean array where 0 means "no event" and 1 means "event". rows = 100 yr window, columns = earthquake
        # rupture
        scenarios = rng.poisson(lambdas, size=(n_samples, lambdas.size))

        # assigns a normal distribution with a mean of 1 and a standard deviation of sd
        # effectively a multiplier for the displacement value
        disp_uncertainty = rng.normal(1., sd, size=(n_samples, lambdas.size))

        # for each 100 yr scenario, get displacements from EQs that happened
        disp_scenarios = scenarios * site_dict_i["disps"]
        # multiplies displacement by the uncertainty multiplier
        disp_scenarios = disp_scenarios * disp_uncertainty
        # sum all displacement values at that site in that 100 yr interval.
        # End up with n-intervals-sized array, with cumulative displacement at that site.
        cumulative_disp_scenarios_total_abs = np.abs(disp_scenarios).sum(axis=1)

        # replace all negative values in disp_scenarios with 0
        disp_scenarios_up = np.where(disp_scenarios > 0, disp_scenarios, 0)
        cumulative_disp_scenarios_up = disp_scenarios_up.sum(axis=1)
        # replace all positive values in disp_scenarios with 0
        disp_scenarios_down = np.where(disp_scenarios < 0, disp_scenarios, 0)
        cumulative_disp_scenarios_down = disp_scenarios_down.sum(axis=1)

        # get displacement thresholds for calculating exceedance (hazard curve x axis)
        thresholds = np.arange(0, 3, 0.01)
        thresholds_neg = thresholds * -1
        # sum all the displacements in the 100 year window that exceed threshold
        n_exceedances_total_abs = np.zeros_like(thresholds)
        n_exceedances_up = np.zeros_like(thresholds)
        n_exceedances_down = np.zeros_like(thresholds)
        # for threshold value:
        for threshold in thresholds:
            # replaces index in zero array with the number of times the cumulative displacement exceeded the threshold
            # across all of the 100 yr scenarios

            # sums the absolute value of the disps if the abs value is greater than threshold. e.g., -0.5 + 0.5 = 1
            n_exceedances_total_abs[thresholds == threshold] = (
                        np.abs(cumulative_disp_scenarios_total_abs) > threshold).sum()
            n_exceedances_up[thresholds == threshold] = (cumulative_disp_scenarios_up > threshold).sum()
        for threshold in thresholds_neg:
            n_exceedances_down[thresholds_neg == threshold] = (cumulative_disp_scenarios_down < threshold).sum()

        # the probability is the number of times that threshold was exceeded divided by the number of samples. so,
        # quite high for low displacements (e.g. 25%). Means there's a ~25% change an earthquake will
        # exceed 0.01 m in next 100 years across all earthquakes in the catalogue (at that site).
        exceedance_probs_total_abs = n_exceedances_total_abs / n_samples
        exceedance_probs_up = n_exceedances_up / n_samples
        exceedance_probs_down = n_exceedances_down / n_samples

        # CAVEAT: at the moment only absolute value thresholds are stored, but for "down" the thresholds are
        # actually negative.
        site_PPE_dict[site_of_interest] = {"thresholds": thresholds,
                                           "exceedance_probs_total_abs": exceedance_probs_total_abs,
                                           "exceedance_probs_up": exceedance_probs_up,
                                           "exceedance_probs_down": exceedance_probs_down,
                                           "site_coords": site_dict_i["site_coords"],
                                           "standard_deviation": sd}

    if extension1 != "":
        with open(f"../{model_version_results_directory}/{extension1}/cumu_exceed_prob_{extension1}"
              f"{taper_extension}.pkl", "wb") as f:
            pkl.dump(site_PPE_dict, f)

    else:
        return site_PPE_dict

def make_fault_model_PPE_dict(branch_weight_dict, model_version_results_directory, slip_taper, n_samples,
                              outfile_extension):
    """ This function takes the branch dictionary and calculates the PPEs for each branch.
    It then combines the PPEs (key = unique branch ID).

    Must run this function with crustal, subduction, or a combination of two.

    :param crustal_branch_dict: from the function make_branch_weight_dict
    :param results_version_directory: string; path to the directory with the solution files
    :return mega_branch_PPE_dictionary and saves a pickle file.
    """

    gf_name = "sites"
    counter = 0

    if slip_taper:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"


    fault_model_allbranch_PPE_dict = {}
    print("Creating fault model PPE dictionary")
    for branch_id in branch_weight_dict.keys():

        print(f"calculating {branch_id} PPE ({counter} of {len(branch_weight_dict.keys())} branches)")
        counter += 1

        # get site displacement dictionary and branch weights
        extension1 = gf_name + branch_weight_dict[branch_id]["file_suffix"]
        branch_weight = branch_weight_dict[branch_id]["total_weight_RN"]

        # Extract rates from the NSHM solution directory, but it is not scaled by the rate scaling factor
        branch_site_disp_dict = get_site_disp_dict(extension1, slip_taper=slip_taper,
                           model_version_results_directory=model_version_results_directory)
        # multiply the rates by the rate scaling factor
        rate_scaling_factor = branch_weight_dict[branch_id]["S"]
        for site in branch_site_disp_dict.keys():
            # multiply each value in the rates array by the rate scaling factor
            branch_site_disp_dict[site]["scaled_rates"] = [rate * rate_scaling_factor for rate in branch_site_disp_dict[
                site]["rates"]]

        ### get exceedance probability dictionary
        branch_cumu_PPE_dict = get_cumu_PPE(branch_site_disp_dict=branch_site_disp_dict,
                    model_version_results_directory=model_version_results_directory, slip_taper=slip_taper,
                    time_interval=100, n_samples=n_samples, extension1="")

        fault_model_allbranch_PPE_dict[branch_id] = {"cumu_PPE_dict": branch_cumu_PPE_dict, "branch_weight":
            branch_weight}

    outfile_name = f"allbranch_PPE_dict{outfile_extension}{taper_extension}"

    with open(f"../{model_version_results_directory}/{outfile_name}.pkl", "wb") as f:
        pkl.dump(fault_model_allbranch_PPE_dict, f)
    return fault_model_allbranch_PPE_dict


def get_weighted_mean_PPE_dict(fault_model_PPE_dict, out_directory, outfile_extension, slip_taper):
    """takes all the branch PPEs and combines them based on the branch weights into a weighted mean PPE dictionary

    :param fault_model_PPE_dict: The dictionary has PPEs for each branch (or branch pairing).
    Each branch contains "branch_weight" and "cumu_PPE_dict".
    "cumu_PPE_dict" is organized by site. Nested in sites is "thresholds", "exceedance_probs_up",
    "exceedance_probs_down", and "exceedance_probs_total_abs"
    :return dictionary of sites, with lists of weighted mean PPEs and threshold values.
    """

    if slip_taper:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    unique_id_list = list(fault_model_PPE_dict.keys())
    site_list = fault_model_PPE_dict[unique_id_list[0]]["cumu_PPE_dict"].keys()

    # weight the probabilities by NSHM branch weights to get a weighted mean
    branch_weights = [fault_model_PPE_dict[unique_id]["branch_weight"] for unique_id in
                      fault_model_PPE_dict.keys()]

    # need a more elegant solution to this I think
    threshold_vals = np.arange(0, 3, 0.01)

    # extract site coordinates from fault model PPE dictionary
    site_coords_dict = {}
    for site in site_list:
        site_coords = fault_model_PPE_dict[unique_id_list[0]]["cumu_PPE_dict"][site]["site_coords"]
        site_coords_dict[site] = site_coords

    weighted_mean_site_probs_dictionary = {}
    for site in site_list:
        weighted_mean_site_probs_dictionary[site] = {}

    for exceed_type in ["total_abs", "up", "down"]:
        for i, site in enumerate(site_list):

            site_df = {}
            for unique_id in unique_id_list:
                probabilities_i_site = fault_model_PPE_dict[unique_id]["cumu_PPE_dict"][site][
                    f"exceedance_probs_{exceed_type}"]
                site_df[unique_id] = probabilities_i_site
            site_probabilities_df = pd.DataFrame(site_df)

            # collapse each row into a weighted mean value
            site_weighted_mean_probs = site_probabilities_df.apply(
                lambda x: np.average(x, weights=branch_weights), axis=1)
            site_max_probs = site_probabilities_df.max(axis=1)
            site_min_probs = site_probabilities_df.min(axis=1)

            weighted_mean_site_probs_dictionary[site]["threshold_vals"] = threshold_vals

            weighted_mean_site_probs_dictionary[site][f"weighted_exceedance_probs_{exceed_type}"] = site_weighted_mean_probs
            weighted_mean_site_probs_dictionary[site][f"{exceed_type}_max_vals"] = site_max_probs
            weighted_mean_site_probs_dictionary[site][f"{exceed_type}_min_vals"] = site_min_probs
            weighted_mean_site_probs_dictionary[site]["site_coords"] = site_coords_dict[site]

    with open(f"../{out_directory}/weighted_mean_PPE_dict{outfile_extension}{taper_extension}.pkl", "wb") as f:
        pkl.dump(weighted_mean_site_probs_dictionary, f)

    return weighted_mean_site_probs_dictionary

def make_sz_crustal_paired_PPE_dict(crustal_branch_weight_dict, sz_branch_weight_dict,
                                    crustal_model_version_results_directory, sz_model_version_results_directory,
                                    slip_taper, n_samples, out_directory, paired_PPE_pickle_name):
    """ This function takes the branch dictionary and calculates the PPEs for each branch.
    It then combines the PPEs (key = unique branch ID).

    Must run this function with crustal, subduction, or a combination of two.

    :param crustal_branch_dict: from the function make_branch_weight_dict
    :param results_version_directory: string; path to the directory with the solution files
    :return mega_branch_PPE_dictionary and saves a pickle file.
    """

    gf_name = "sites"

    if slip_taper:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"
    #paired_PPE_pickle_name = f"sz_crustal_paired_PPE_dict{outfile_extension}{taper_extension}.pkl"

    # make a dictionary of displacements at each site from all the crustal earthquake scenarios

    all_crustal_branches_site_disp_dict = {}
    for branch_id in crustal_branch_weight_dict.keys():

        extension1 = gf_name + crustal_branch_weight_dict[branch_id]["file_suffix"]
        # get site displacement dictionary
        # this extracts the rates from the solution directory, but it is not scaled by the rate scaling factor
        crustal_branch_site_disp_dict = get_site_disp_dict(
            extension1, slip_taper=slip_taper, model_version_results_directory=crustal_model_version_results_directory)

        # multiply the rates by the rate scaling factor
        rate_scaling_factor = crustal_branch_weight_dict[branch_id]["S"]
        for site in crustal_branch_site_disp_dict.keys():
            # multiply each value in the rates array by the rate scaling factor
            crustal_branch_site_disp_dict[site]["scaled_rates"] = \
                [rate * rate_scaling_factor for rate in crustal_branch_site_disp_dict[site]["rates"]]

        all_crustal_branches_site_disp_dict[branch_id] = {"site_disp_dict":crustal_branch_site_disp_dict,
                                                   "branch_weight":crustal_branch_weight_dict[branch_id][
                                                       "total_weight_RN"]}

    # make a dictionary of displacements at each site from all the crustal earthquake scenarios
    all_sz_branches_site_disp_dict = {}
    for branch_id in sz_branch_weight_dict.keys():
        sz_slip_taper = False

        extension1 = gf_name + sz_branch_weight_dict[branch_id]["file_suffix"]
        # get displacement dictionary
        # this extracts the rates from the solution directory, but it is not scaled by the rate scaling factor
        sz_branch_site_disp_dict = get_site_disp_dict(
            extension1, slip_taper=sz_slip_taper, model_version_results_directory=sz_model_version_results_directory)

        # multiply the rates by the rate scaling factor
        rate_scaling_factor = sz_branch_weight_dict[branch_id]["S"]
        for site in sz_branch_site_disp_dict.keys():
            # multiply each value in the rates array by the rate scaling factor
            sz_branch_site_disp_dict[site]["scaled_rates"] = \
                [rate * rate_scaling_factor for rate in sz_branch_site_disp_dict[site]["rates"]]

        all_sz_branches_site_disp_dict[branch_id] = {"site_disp_dict":sz_branch_site_disp_dict,
                                                   "branch_weight":sz_branch_weight_dict[branch_id]["total_weight_RN"]}

    # make all the combinations of crustal and subduction zone branch pairs
    crustal_sz_branch_pairs = list(itertools.product(crustal_branch_weight_dict.keys(),
                                                      sz_branch_weight_dict.keys()))

    counter = 0
    paired_crustal_sz_PPE_dict = {}
    print("Making the sz_crustal_paired dictionary")
    for pair in crustal_sz_branch_pairs:
        # get the branch unique ID for the crustal and sz combos
        crustal_unique_id, sz_unique_id = pair[0], pair[1]
        pair_unique_id = crustal_unique_id + "_" + sz_unique_id

        print(f"calculating {pair_unique_id} PPE ({counter} of {len(crustal_sz_branch_pairs)} branches)")
        counter += 1

        site_names = list(all_crustal_branches_site_disp_dict[crustal_unique_id]["site_disp_dict"].keys())

        pair_weight = all_crustal_branches_site_disp_dict[crustal_unique_id]["branch_weight"] * \
                       all_sz_branches_site_disp_dict[sz_unique_id]["branch_weight"]

        # loop over all the sites for the crustal and sz branches of interest
        # make one long list of displacements and corresponding scaled rates per site
        pair_site_disp_dict = {}
        for j, site in enumerate(site_names):
            site_coords = all_crustal_branches_site_disp_dict[crustal_unique_id]["site_disp_dict"][site]["site_coords"]

            crustal_site_disps = all_crustal_branches_site_disp_dict[crustal_unique_id]["site_disp_dict"][site]["disps"]
            sz_site_disps = all_sz_branches_site_disp_dict[sz_unique_id]["site_disp_dict"][site]["disps"]

            crustal_site_scaled_rates = all_crustal_branches_site_disp_dict[crustal_unique_id]["site_disp_dict"][
                site]["scaled_rates"]
            sz_site_scaled_rates = all_sz_branches_site_disp_dict[sz_unique_id]["site_disp_dict"][site]["scaled_rates"]

            pair_site_disps = crustal_site_disps + sz_site_disps
            pair_scaled_rates = crustal_site_scaled_rates + sz_site_scaled_rates

            pair_site_disp_dict[site] = {"disps": pair_site_disps, "scaled_rates": pair_scaled_rates,
                                                 "site_coords": site_coords}


        # get exceedence probabilities for each crustal/sz pair

        if not os.path.exists(f"../{out_directory}"):
            os.mkdir(f"../{out_directory}")
        pair_cumu_PPE_dict = get_cumu_PPE(branch_site_disp_dict=pair_site_disp_dict,
                                            model_version_results_directory=out_directory,
                                            slip_taper=slip_taper, time_interval=100,
                                            n_samples=n_samples, extension1="")

        paired_crustal_sz_PPE_dict[pair_unique_id] = {"cumu_PPE_dict": pair_cumu_PPE_dict, "branch_weight": pair_weight}

    with open(f"../{out_directory}/{paired_PPE_pickle_name}", "wb") as f:
        pkl.dump(paired_crustal_sz_PPE_dict, f)
    return paired_crustal_sz_PPE_dict

def get_exceedance_bar_chart_data(site_PPE_dictionary, probability, exceed_type, site_list):
    """returns displacements at the X% probabilities of exceedance for each site

    define exceedance type. Options are "total_abs", "up", "down"
    """

    # get disp threshold (x-value) at defined probability (y-value)
    disps = []
    for site in site_list:
        threshold_vals = site_PPE_dictionary[site]["thresholds"]

        # dispalcement threasholds are negative for "down" exceedances
        if exceed_type == "down":
            threshold_vals = -threshold_vals

        site_PPE = site_PPE_dictionary[site][f"exceedance_probs_{exceed_type}"]

        # get first index that is < 10% (ideally we would interpolate for exact value but don't have a function)
        exceedance_index = next((index for index, value in enumerate(site_PPE) if value <= probability), -1)
        disp = threshold_vals[exceedance_index]
        disps.append(disp)

    return disps

def get_probability_bar_chart_data(site_PPE_dictionary, exceed_type, threshold, site_list):
    """ function that finds the probability at each site for the specified displacement threshold on the hazard curve
        Inputs:
        :param: dictionary of exceedance probabilities for each site (key = site)
        :param exceedance type: string; "total_abs", "up", or "down"
        :param: list of sites to get data for. If None, will get data for all sites in site_PPE_dictionary.
                I made this option so that you could skip the sites you didn't care about (e.g., use "plot_order")

        Outputs:
        :return    probs_threshold: list of probabilities of exceeding the specified threshold (one per site)
            """


    # get list of probabilities at defined displacement threshold (one for each site)
    probs_threshold = []
    for site in site_list:
        site_PPE = site_PPE_dictionary[site][f"exceedance_probs_{exceed_type}"]
        threshold_vals = list(site_PPE_dictionary[site]["thresholds"])

        # find index in threshold_vals where the value matches the parameter threshold
        index = threshold_vals.index(threshold)
        probs_threshold.append(site_PPE[index])

    return probs_threshold


def plot_branch_hazard_curve(extension1, slip_taper, model_version_results_directory, file_type_list,
                             skipped_sites):
    """makes hazard curves for each site. includes the probability of cumulative displacement from multiple
    earthquakes exceeding a threshold in 100 years."""

    exceed_type_list = ["total_abs", "up", "down"]

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    with open(f"../{model_version_results_directory}/{extension1}/cumu_exceed_prob_{extension1}"
              f"{taper_extension}.pkl",
              "rb") as fid:
        PPE_dictionary = pkl.load(fid)

    # make a list of sites to plot in a specific order. Can skip sites you don't want to plot.
    plot_order = get_plot_order_list(dictionary=PPE_dictionary, skipped_sites=skipped_sites)

    plt.close("all")
    fig, axs = plt.subplots(figsize=(8, 10))
    plt.subplots_adjust(hspace=0.3, wspace=0.3)
    fig.suptitle("Sites", fontsize=18, y=0.95)

    #loop through sites and make a subplot for each one
    for i, site in enumerate(plot_order):
        ax = plt.subplot(4, 3, i + 1)

        # plots all three types of exceedance (total_abs, up, down) on the same plot
        for j, exceed_type in enumerate(exceed_type_list):
            curve_color = get_probability_color(exceed_type)

            exceedance_probs = PPE_dictionary[site][f"exceedance_probs_{exceed_type}"]
            threshold_vals = PPE_dictionary[site]["thresholds"]

            ax.plot(threshold_vals, exceedance_probs, color=curve_color)
            ax.axhline(y=0.02, color="0.7", linestyle='dashed')
            ax.axhline(y=0.1, color="0.7", linestyle='dotted')

        ax.set_title(site)
        ax.set_yscale('log'), ax.set_xscale('log')
        ax.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1])
        xmin, xmax = 0.01, 3
        ymin, ymax = 0.000005, 1
        ax.set_ylim([ymin, ymax])
        ax.set_xlim([xmin, xmax])
        ax.get_xaxis().set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(axis='x', style='plain')
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    fig.text(0.5, 0, 'Vertical displacement threshold (m)', ha='center')
    fig.text(0, 0.5, 'Probability of exceedance in 100 years', va='center', rotation='vertical')
    fig.suptitle(f"cumulative exceedance hazard curves \n{taper_extension}")
    plt.tight_layout()

    #save hazard curve figure
    # make directory for hazard curve if it doesn't exist
    if not os.path.exists(f"../{model_version_results_directory}/{extension1}/probability_figures"):
        os.mkdir(f"../{model_version_results_directory}/{extension1}/probability_figures")

    for file_type in file_type_list:
        plt.savefig(f"../{model_version_results_directory}/{extension1}/probability_figures/hazard_curve_{extension1}"
                    f"{taper_extension}.{file_type}", dpi=300)

def make_10_2_disp_plot(extension1, slip_taper, model_version_results_directory, skipped_sites,
                        file_type_list=["png", "pdf"]):
    """ makes bar charts of the displacement value at the 10% and 2% probability of exceence thresholds for each site
        extension1 = "sites_c_MDEz" or whatever
        fault_type = "crustal" or "sz"
        slip_taper = True or False
        model_version_results_directory = "{results_directory}/{fault_type}{fault_model}"
    """

    probability_list = [0.1, 0.02]

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"


    with open(f"../{model_version_results_directory}/{extension1}/cumu_exceed_prob_{extension1}"
              f"{taper_extension}.pkl", "rb") as fid:
        site_PPE_dictionary = pkl.load(fid)

    plot_order = get_plot_order_list(dictionary=site_PPE_dictionary, skipped_sites=skipped_sites)

    plt.close("all")
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.4))
    x = np.arange(len(plot_order))  # the site label locations
    width = 0.4  # the width of the bars
    # find maximum value in all the "up" columns in PPE dictionary

    max_min_y_vals = []
    for i, probability in enumerate(probability_list):
        disps_up = \
            get_exceedance_bar_chart_data(site_PPE_dictionary=site_PPE_dictionary, exceed_type="up",
                                    site_list=plot_order, probability=probability)
        disps_down= \
            get_exceedance_bar_chart_data(site_PPE_dictionary=site_PPE_dictionary, exceed_type="down",
                                     site_list=plot_order, probability=probability)

        max_min_y_vals.append(max(disps_up))
        max_min_y_vals.append(min(disps_down))

        color_up = (189/255, 0, 0)
        color_down = (15/255, 72/255, 186/255)
        label_size = 6
        label_offset = label_size / 60

        # add bars to plot, add black horizontal line at zero.
        bars_up = axs[i].bar(x, disps_up, width, color=color_up, linewidth=0.5)
        bars_down = axs[i].bar(x, disps_down, width, color=color_down, linewidth=0.5)
        axs[i].axhline(y=0, color="k", linewidth=0.5)

        # add value labels to bars
        for bar in bars_up:
            bar_color = bar.get_facecolor()
            axs[i].text(x=bar.get_x(), y=bar.get_height() + label_offset, s=round(bar.get_height(), 1), ha='left',
                        va='center', color=bar_color, fontsize=label_size, fontweight='bold')

        for bar in bars_down:
            bar_color = bar.get_facecolor()
            axs[i].text(x=bar.get_x(), y=bar.get_height() - label_offset, s=round(bar.get_height(), 1), ha='left',
                        va='center', color=bar_color, fontsize=label_size, fontweight='bold')

    for i in range(len(probability_list)):
        axs[i].set_ylim(min(max_min_y_vals) - 0.25, max(max_min_y_vals) + 0.25)
        axs[i].tick_params(axis='x', labelrotation=90, labelsize=label_size)
        axs[i].set_xticks(x, plot_order)
        axs[i].tick_params(axis='y', labelsize=8)
        axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # set tick labels to be every 0.2
        axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.5))

    # set indidual subplot stuff
    axs[0].set_ylabel("Minimum displacement (m)", fontsize=8)
    axs[1].tick_params(axis='y', labelleft=False)

    axs[0].set_title(f"10% probability of exceedance", fontsize=8)
    axs[1].set_title(f"2% probability of exceedance", fontsize=8)

    # manually make legend with rectangles and text
    swatch_width, swatch_height = width, max(max_min_y_vals) * 0.08
    swatch_minx, swatch_miny = -1 * (len(plot_order) / 30), max(max_min_y_vals)
    axs[0].add_patch(Rectangle((swatch_minx, swatch_miny), swatch_width, swatch_height,
                               facecolor=color_up, edgecolor=None))
    axs[0].add_patch(Rectangle((swatch_minx, swatch_miny - 2 * swatch_height), swatch_width, swatch_height,
                               facecolor=color_down, edgecolor=None))


    axs[0].text(swatch_minx + 2 * swatch_width, swatch_miny, "uplift", fontsize=8)
    axs[0].text(swatch_minx + 2 * swatch_width, swatch_miny - 2 * swatch_height, "subsidence", fontsize=8)

    fig.suptitle(f"100 yr exceedance minimum displacements\n{extension1}{taper_extension}", fontsize=10)
    fig.tight_layout()

    # make a directory for the figures if it doesn't already exist
    outfile_directory = f"../{model_version_results_directory}/{extension1}/probability_figures"
    if not os.path.exists(f"{outfile_directory}"):
        os.makedirs(f"{outfile_directory}")
    for file_type in file_type_list:
        fig.savefig(f"{outfile_directory}/10_2_disps_{extension1}{taper_extension}.{file_type}", dpi=300)


# What is the probability of exceeding 0.2 m subsidence, 0.2 m uplift at each site?
def make_prob_bar_chart(extension1,  slip_taper, model_version, model_version_results_directory, skipped_sites,
                        threshold=0.2):
    # What is the probability of exceeding 0.2 m subsidence, 0.2 m uplift at each site?
    """ determines the probability of exceeding a defined displacement threshold at each site and plots as a bar chart
        two-part plot, one for up and one for down. y axis is probability, x axis is site name
        :param extension1: string, name of the NSHM branch suffix etc.
        :param slip_taper: boolean, True if slip tapers, False if uniform slip
        :param fault_type: string, "crustal" or sz"
        :param threshold: float, displacement threshold to determine exceedance probability
        :param results_directory: string, name of directory where results are stored
    """

    exceed_type_list = ["up", "down"]

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    with open(f"../{model_version_results_directory}/{extension1}/cumu_exceed_prob_{extension1}"
              f"{taper_extension}.pkl",
              "rb") as fid:
        site_PPE_dictionary = pkl.load(fid)

    # make list of sites to plot, skipping sites if necessary
    plot_order = get_plot_order_list(dictionary=site_PPE_dictionary, skipped_sites=skipped_sites)

    # set up custom color scheme
    colors = make_qualitative_colormap("custom", len(plot_order))

    # set up figure and subplots
    fig, axs = plt.subplots(1, 2, figsize=(7, 5))
    x = np.arange(len(plot_order))  # the site label locations
    width = 0.6  # the width of the bars

    for i, exceed_type in enumerate(exceed_type_list):
        probs_threshold_exceed_type = \
            get_probability_bar_chart_data(site_PPE_dictionary=site_PPE_dictionary, exceed_type=exceed_type,
                                      threshold=threshold, site_list=plot_order)

        # add bars to plot
        bars_10cm = axs[i].bar(x, probs_threshold_exceed_type, width, color=colors)

        for bar in bars_10cm:
            bar_color = bar.get_facecolor()
            # add value label to each bar
            axs[i].text(x=bar.get_x() + bar.get_width() / 2, y=bar.get_height() + 0.03,
                        s=f"{int(100 * round(bar.get_height(), 2))}%", horizontalalignment='center', color=bar_color,
                        fontsize=6, fontweight='bold')

        axs[i].set_ylim(0.0, 0.5)
        axs[i].tick_params(axis='x', labelrotation=90, labelsize=6)
        axs[i].tick_params(axis='y', labelsize=8)
        axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # set tick labels to be every 0.2
        axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.2))
        axs[i].set_xticks(x, plot_order)

    # set indidual subplot stuff
    axs[0].set_ylabel("Probabilty", fontsize=8)
    axs[1].tick_params(axis='y', labelleft=False)

    axs[0].set_title(f"Probability of exceeding {threshold} m uplift", fontsize=8)
    axs[1].set_title(f"Probability of exceeding {threshold} m subsidence", fontsize=8)

    fig.suptitle(f"{model_version} faults (100 yrs)")
    fig.tight_layout()

    # make directory for hazard curve if it doesn't exist
    if not os.path.exists(f"../{model_version_results_directory}/{extension1}/probability_figures"):
        os.mkdir(f"../{model_version_results_directory}/{extension1}/probability_figures")
    fig.savefig(f"../{model_version_results_directory}/{extension1}/probability_figures/prob_bar_chart_{extension1}"
                f"{taper_extension}.pdf", dpi=300)
    fig.savefig(f"../{model_version_results_directory}/{extension1}/probability_figures/prob_bar_chart_{extension1}"
                f"{taper_extension}.png", dpi=300)

def make_branch_prob_plot(extension1,  slip_taper, model_version, model_version_results_directory, skipped_sites,
                      file_type_list=["png", "pdf"], threshold=0.2):
    """ """

    exceed_type_list = ["up", "down"]

    if slip_taper is True:
        taper_extension = "_tapered"
    else:
        taper_extension = "_uniform"

    with open(f"../{model_version_results_directory}/{extension1}/cumu_exceed_prob_{extension1}"
              f"{taper_extension}.pkl",
              "rb") as fid:
        PPE_dict = pkl.load(fid)

    plot_order = get_plot_order_list(PPE_dict, skipped_sites=skipped_sites)

    # set up custom color scheme
    # set up custom color scheme
    colors = make_qualitative_colormap("custom", len(plot_order))
    point_size = [35]

    # set up figure and subplots
    fig, axs = plt.subplots(1, 2, figsize=(7, 3.5))
    x = np.arange(len(plot_order))  # the site label locations

    for i, exceed_type in enumerate(exceed_type_list):
        probs = \
            get_probability_bar_chart_data(site_PPE_dictionary=PPE_dict, exceed_type=exceed_type,
                                           threshold=threshold, site_list=plot_order)

        # add point and error bars to plot
        axs[i].scatter(x, probs, s=point_size, color=colors, zorder=3, edgecolors='k', linewidths=0.5)

        labels = [f"{int(100 * round(prob, 2))}%" for prob in probs]
        label_y_vals = [prob + 0.03 for prob in probs]
        for site, q in enumerate(x):
            axs[i].text(x=x[q], y=label_y_vals[q], s=labels[q],
                        horizontalalignment='center', fontsize=6, fontweight='bold')

        ymin, ymax  = 0.0, 0.3
        axs[i].set_ylim([ymin, ymax])
        axs[i].tick_params(axis='x', labelrotation=90, labelsize=6)
        axs[i].tick_params(axis='y', labelsize=8)
        axs[i].yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        # set tick labels to be every 0.2
        axs[i].yaxis.set_major_locator(mticker.MultipleLocator(0.1))
        axs[i].set_xticks(x, plot_order, va='top', ha='center')

    # set indidual subplot stuff
    fontsize = 8
    # I'm doing it this way instead of just using "names" because I want to make sure the legend is in the correct
    # order.

    axs[0].set_ylabel("Probabilty", fontsize=8)
    axs[1].tick_params(axis='y', labelleft=False)

    axs[0].set_title(f"Probability of exceeding {threshold} m uplift", fontsize=fontsize)
    axs[1].set_title(f"Probability of exceeding {threshold} m subsidence", fontsize=fontsize)

    fig.suptitle(f"{model_version} {extension1} (100 yrs)")
    fig.tight_layout()

    outfile_directory = f"../{model_version_results_directory}/{extension1}/probability_figures"
    if not os.path.exists(f"{outfile_directory}"):
        os.mkdir(f"{outfile_directory}")

    for file_type in file_type_list:
        fig.savefig(f"{outfile_directory}/probs_chart_{extension1}{taper_extension}.{file_type}", dpi=300)

