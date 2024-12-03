
import os
import matplotlib
from glob import glob
from pcdhm.weighted_mean_plotting import map_and_plot_probabilities
from pcdhm.compare_fault_model import compare_faultmodel_prob_plot, compare_disps_chart, compare_mean_hazcurves, \
    compare_disps_with_net

########## USER INPUTS #######################
plot_order_name = "12 sites"                 # "12 sites" or "pared" or "porirua"
results_directory = "results_EXAMPLE2"
exceed_type = "down"                     # "down", "up", or "total_abs"

# Choose what models to compare. These names should be in the results folder already.
model_subdirectory_names = ['sz_multi50', 'crustal_CFM']

# used for plot labels/titles. must be in same order as model_subdirectory_names
pretty_names = model_subdirectory_names

# plot crustal logic tree branch parameters.
# Must already have calculated the specific "weighted_mean_PPE_dict" in step 07
# set model subdirectory name (above) to appropriate crustal directory
plot_branch_sensitivity = False
sensitivity_param = 'time_dependency'  #time_dependency, 'deformation_model' 's_value', 'b_n_value'

file_type_list = ["png"]     # generally png and/or pdf
probability_plot = True             # plots the probability of exceedance at the 0.2 m uplift and subsidence thresholds
displacement_chart = False           # plots the displacement at the 10% and 2% probability of exceedance thresholds
compare_hazcurves = True        # plots the different hazard curves on the same plot
make_map = False
disps_net = False
labels_on = False                # displacement number labels for bar charts and probability plots


#### script ###################
if plot_branch_sensitivity is True:
    if sensitivity_param == 'time_dependency':
        branch_params = ['TI', 'TD']
    elif sensitivity_param == 'deformation_model':
        branch_params = ['geologic', 'geodetic']
    elif sensitivity_param == 's_value':
        branch_params = ['S066', 'S10', 'S141']
    elif sensitivity_param == 'b_n_value':
        branch_params = ['N46_b1089', 'N34_b0959', 'N27_b0823']
    else:
        # print an exit statement with a string error
        raise ValueError('the sensitivity parameter is wrong')
        # exit the script

    pretty_names = branch_params

# makes the text edible upon export to a pdf
matplotlib.rcParams['pdf.fonttype'] = 42

displacement_threshold_list = [0.2]

if not plot_branch_sensitivity:
    title = " vs ".join(pretty_names)
    file_name = "_".join(pretty_names)
    file_name = file_name.replace(" ", "_")
else:
    title = f"{model_subdirectory_names[0]} {sensitivity_param}"
    file_name = f'{sensitivity_param}_{exceed_type}'


mean_PPE_path_list = []
if plot_branch_sensitivity is False:
    compare_results_directory = f"{results_directory}/compare_fault_models"
    outfile_directory = f"{compare_results_directory}/{file_name}"
    for name in model_subdirectory_names:
        mean_PPE_path_i = glob(f"../{results_directory}/{name}/weighted_mean_PPE_dict*.pkl")
        mean_PPE_path_list.append(mean_PPE_path_i[0])
else:   # can only be one object in model_subdirectory_names
    compare_results_directory = f"{results_directory}/compare_branch_sensitivity"
    outfile_directory = f"{compare_results_directory}/{model_subdirectory_names[0]}_{sensitivity_param}"
    for branch_param in branch_params:
        mean_PPE_path_i = glob(f"../{results_directory}/{model_subdirectory_names[0]}/weighted_mean_PPE_dict_{branch_param}_*.pkl")
        mean_PPE_path_list.append(mean_PPE_path_i[0])


if not os.path.exists(f"../{compare_results_directory}"):
        os.makedirs(f"../{compare_results_directory}")

if not os.path.exists(f"../{outfile_directory}"):
        os.makedirs(f"../{outfile_directory}")

if plot_order_name == "12 sites":
    plot_order = ["Paraparaumu", "Porirua CBD north", "South Coast", "Wellington Airport", "Wellington CBD", "Petone",
                   "Seaview", "Eastbourne", "Turakirae Head", "Lake Ferry", "Cape Palliser", "Flat Point"]

if plot_order_name == "pared":
    plot_order = ["South Coast", "Porirua CBD north",  "Petone", "Eastbourne", "Turakirae Head", "Lake Ferry",
                    "Cape Palliser", "Flat Point"]

if plot_order_name == "porirua":
    plot_order = ["Porirua CBD north", "Porirua CBD south"]


if probability_plot:
    if exceed_type == "total_abs":
        exceed_type_list = ["total_abs"]
    else:
        exceed_type_list = ["up", "down"]
    compare_faultmodel_prob_plot(PPE_paths=mean_PPE_path_list,
                                 plot_name=file_name,
                                 outfile_directory=outfile_directory,
                                 title=title,
                                 pretty_names=pretty_names,
                                 exceed_type_list=exceed_type_list,
                                 plot_order=plot_order,
                                 labels_on=labels_on,
                                 file_type_list=file_type_list,
                                 threshold=0.2)


if displacement_chart:
    compare_disps_chart(PPE_paths=mean_PPE_path_list,
                        plot_name=file_name,
                        outfile_directory=outfile_directory,
                        title=title,
                        pretty_names=pretty_names,
                        plot_order=plot_order,
                        labels_on=labels_on,
                        file_type_list=file_type_list)

if compare_hazcurves:
    compare_mean_hazcurves(PPE_paths=mean_PPE_path_list, plot_name=file_name, outfile_directory=outfile_directory,
                           title=title, pretty_names=pretty_names, exceed_type=exceed_type,
                           plot_order=plot_order,
                           file_type_list=file_type_list)

if disps_net:
    compare_disps_with_net(PPE_paths=mean_PPE_path_list, plot_name=file_name, outfile_directory=outfile_directory,
                           title=title, pretty_names=pretty_names,
                           file_type_list=file_type_list)
if make_map:
    PPE_dicts = []
    for PPE_path in mean_PPE_path_list:
        for disp in displacement_threshold_list:
            map_and_plot_probabilities(PPE_path=PPE_path,
                                       plot_name=file_name,
                                       exceed_type=exceed_type,
                                       title=title,
                                       outfile_directory=outfile_directory,
                                       plot_order=plot_order,
                                       labels_on=True,
                                       file_type_list=file_type_list,
                                       threshold=disp,
                                       colorbar_max=0.35)
