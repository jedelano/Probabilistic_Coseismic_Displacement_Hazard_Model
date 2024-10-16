
import os
import matplotlib
from pcdhm.weighted_mean_plotting import map_and_plot_probabilities
from pcdhm.compare_fault_model import compare_faultmodel_prob_plot, compare_disps_chart, compare_mean_hazcurves, \
    compare_disps_with_net

########## USER INPUTS #######################
plot_order_name = "porirua"                 # "12 sites" or "pared" or "porirua"
results_directory = "results_r1"
exceed_type = "down"                     # "down", "up", or "total_abs"
slip_taper = False

# Choose what models to compare. These names should be in the results folder already.
model_subdirectory_names = ['crustal_CFM', 'crustal_Model1', 'crustal_Model2']
#model_subdirectory_names = ['sz_multi50', 'sz_multi50_steeperdip', 'sz_multi50_gentlerdip']

# used for plot labels/titles. must be in same order as model_subdirectory_names
pretty_names = model_subdirectory_names

file_type_list = ["png", "pdf"]     # generally png and/or pdf
probability_plot = True             # plots the probability of exceedance at the 0.2 m uplift and subsidence thresholds
displacement_chart = True           # plots the displacement at the 10% and 2% probability of exceedance thresholds
compare_hazcurves = True        # plots the different hazard curves on the same plot
make_map = False
disps_net = False
labels_on = False                # displacement number labels for bar charts and probability plots


#### script ###################
# makes the text edible upon export to a pdf
matplotlib.rcParams['pdf.fonttype'] = 42

displacement_threshold_list = [0.2]

title = " vs ".join(pretty_names)
file_name = "_".join(pretty_names)
file_name = file_name.replace(" ", "_")

if slip_taper: slip_taper_extension = "_tapered"
else: slip_taper_extension = "_uniform"

mean_PPE_path_list = []
for name in model_subdirectory_names:
    mean_PPE_path_i = f"../{results_directory}/{name}/weighted_mean_PPE_dict{slip_taper_extension}.pkl"
    mean_PPE_path_list.append(mean_PPE_path_i)

compare_results_directory = f"{results_directory}/compare_fault_models"
if not os.path.exists(f"../{compare_results_directory}"):
        os.makedirs(f"../{compare_results_directory}")

outfile_directory = f"{compare_results_directory}/{file_name}"
if not os.path.exists(f"../{outfile_directory}"):
        os.makedirs(f"../{outfile_directory}")

if plot_order_name == "12 sites":
    plot_order = ["Paraparaumu", "Porirua CBD north", "South Coast", "Wellington Airport", "Wellington CBD", "Petone",
                   "Seaview", "Eastbourne", "Turakirae Head", "Lake Ferry", "Cape Palliser", "Flat Point"]

if plot_order_name == "pared":
    plot_order = ["South Coast", "Porirua CBD north",  "Petone", "Eastbourne", "Turakirae Head", "Lake Ferry",
                    "Cape Palliser", "Flat Point"]

if plot_order_name == "porirua":
    plot_order = ["Porirua CBD north","Porirua CBD south"]


if probability_plot:
    compare_faultmodel_prob_plot(PPE_paths=mean_PPE_path_list, plot_name=file_name,
                                 outfile_directory=outfile_directory, title=title, pretty_names=pretty_names,
                                 plot_order=plot_order,
                                 labels_on=labels_on,
                                 file_type_list=file_type_list,
                                 threshold=0.2)

if displacement_chart:
    compare_disps_chart(PPE_paths=mean_PPE_path_list, plot_name=file_name, outfile_directory=outfile_directory,
                        title=title, pretty_names=pretty_names,
                        plot_order=plot_order,
                        labels_on=labels_on, file_type_list=file_type_list)

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
                                       colorbar_max=0.3)
