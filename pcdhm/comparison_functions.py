
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

