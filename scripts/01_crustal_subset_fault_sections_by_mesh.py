import geopandas as gpd
import os

# This script extracts the fault sections from the NSHM that are located along meshes in the Wellington area.
# Uses keywords that you have to input. Outputs a geojson with filtered traces.
# only have to do once no matter the NSHM branch since the faults don't change

#######INPUTS
# find only the rupture scenarios that use faults we have meshes for
NSHM_directory = "NZSHM22_InversionSolution-QXV0b21hdGlvblRhc2s6MTA3MDEz"  # geologic, mid B and N, C 4.2
crustal_mesh_version = "CFM"         # "_Model1" or "_Model2" or "_CFM
out_files_directory = "mesh_gf_outfiles_EXAMPLE" # name the outfiles location, same for subsequent scripts.

#############

# This system isn't great, but it's fine for now. Choose a unique part of each name in fault sections,
# avoid using hyphens and stuff to remove ambiguity. It's ok to "capture" more faults than you have meshes for at this
# stage.
target_fault_names = ["Aotea", "Huangarua", "Fisherman", "Honeycomb", "Moonshine", "Otaki", "Ohariu", "Okupe",
                      "Opouawe", "Uruti", "Otaraia", "Pahaua", "Palliser", "Pukerua", "Riversdale", "Shepherds Gully",
                      "Mana", "Otaheke",  "Wairarapa", "Wellington Hutt", "Whareama", "Wharekauhau", "Whitemans"]

# load all fault sections
traces = gpd.read_file(f"../data/crustal_solutions/{NSHM_directory}/ruptures/fault_sections.geojson"
                       "").to_crs(epsg=2193)

# find all fault ids that have a name match to the target name list
filtered_trace_ids = []
for i, row in traces.iterrows():
    for fault_name in target_fault_names:
        if fault_name in row.FaultName:
            filtered_trace_ids.append(row.FaultID)

# remove duplicates and sort
fault_set = set(filtered_trace_ids)
filtered_trace_ids = list(fault_set)
filtered_trace_ids.sort()

# subset fault sections by traget rupture id
filtered_traces_gdf = traces[traces.FaultID.isin(filtered_trace_ids)]

# make directory for outputs if it doesn't already exist
if not os.path.exists(f"../{out_files_directory}"):
    os.makedirs(f"../{out_files_directory}")
if not os.path.exists(f"../{out_files_directory}/crustal_{crustal_mesh_version}"):
    os.mkdir(f"../{out_files_directory}/crustal_{crustal_mesh_version}")

filtered_traces_gdf.to_file(f"../{out_files_directory}/crustal_{crustal_mesh_version}"
                            f"/name_filtered_fault_sections.geojson",
                            driver="GeoJSON")

