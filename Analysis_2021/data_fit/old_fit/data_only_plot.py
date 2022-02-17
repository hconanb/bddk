import os

gc_onflag = 0
fit_strat = "dp_f"

#spec is spectrum, tuples are (mc_id, fit shape)
specs = ["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]
#,"zz","p","m","st"]

run_name = f"zo_{fit_strat}_{gc_onflag}"
# plot_data(run_name, specs)

sl = " ".join(str(x) for x in specs)

os.system(f"python3 data_plot_script.py\
            --run_name {run_name}\
            --fit_strat {fit_strat}\
            --specs {sl}\
            --data_only_flag 1\
            ")
