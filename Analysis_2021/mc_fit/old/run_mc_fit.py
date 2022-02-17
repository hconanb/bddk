import os


for nms in mc_event_shape_list:
    el = " ".join(str(x) for x in nms[1])
    sl = " ".join(str(x) for x in nms[2])
    os.system(f"python3 mc_fit/mc_ws_build_script.py\
                --ws_name {nms[0]}\
                --event_list {el}\
                --shape_list {sl}\
                --mean_guess {nms[3]}\
                --mean_window {nms[4]}\
                --fit_window {nms[5]}\
                2>&1 | tee mc_fit/log_files/{nms[0]}_base_mc_output.txt")
#
for nms in mc_event_shape_list:
    for shape in nms[2]:
        os.system(f"python3 mc_fit/mc_ws_fit_script.py\
                    --ws_name {nms[0]}\
                    --shape {shape}\
                    2>&1 | tee mc_fit/log_files/{nms[0]}_fit_mc_{shape}_output.txt")
        os.system(f"python3 mc_fit/mc_plot_script.py\
                    --ws_name {nms[0]}\
                    --shape {shape}\
                    2>&1 | tee mc_fit/log_files/{nms[0]}_plot_mc_{shape}_output.txt")
