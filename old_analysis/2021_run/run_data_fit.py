import os

#Parameters from efficiecies, known branching fractions are taken as constants, or are floating as gaussian constraints
# gc_onflag = 0 or 1

#Fit to the data using mc as a starting point.
# dp indicates we fit to the combined peaks
# mcp indicates we fit to the individual mc peaks
#_c indicates we constraine to the MC fit shapes

# fit_strat = "super, dp_f, df_c, mcp_f, mcp_c "

# "zz": ["9", "7", "10", "4", "8", "12"]
# "zz": ["DG", "BGEP",  "BGEP", "BGEP", "cb2", "cb2"]}

gc_onflag = 0
fit_strat = "dp_f_nn"

#spec is spectrum, tuples are (mc_id, fit shape)
specs = ["z","zz","p","m","st"]

run_name = f"z_{fit_strat}_{gc_onflag}"
# plot_data(run_name, specs)

sl = " ".join(str(x) for x in specs)

# os.system(f"python3 data_fit/data_ws_build_script.py\
#             --run_name {run_name}\
#             --gc_onflag {gc_onflag}\
#             --fit_strat {fit_strat}\
#             --specs {sl}\
#             2>&1 | tee data_fit/log_files/{run_name}_base_data_output.txt")
#
# os.system(f"python3 data_fit/data_ws_fit_script.py\
#             --run_name {run_name}\
#             --gc_onflag {gc_onflag}\
#             --specs {sl}\
#             2>&1 | tee data_fit/log_files/{run_name}_fit_data_output.txt")

os.system(f"python3 data_fit/data_plot_script.py\
            --run_name {run_name}\
            --fit_strat {fit_strat}\
            --specs {sl}\
            2>&1 | tee data_fit/log_files/{run_name}_plot_data_output.txt")
