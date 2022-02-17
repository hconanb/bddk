import os

bwindow = 20

p0s = 5280
p0sw = 60

p1s = 5140
p1sw = 60

p2s = 4990
p2sw = 60

mc_event_shape_list = [

    #spectrum_id, [List of MC files], [List of shapes to fit], "B Mean Guess", "B_Mean_Guess_Window", "B Fit Window"
    ["01_Z_m_p", ["01_Z_m_p_11198006"], ["DG"], 5280, 10, 30],
    ["02_Z_m_p", ["02_Z_m_p_11198400"], ["BGEP"], 5130, 10, 50],
    # ["z_03", ["02_z_11198400"], ["GEP","BGEP"], 5130, 10, 50],
    ["04_Z_m_p", ["04_Z_m_p_11198401"], ["BG"], 4990, 10, 50],

    # ["zz_09",["09_zz_11196019"], ["GEP","BG"], 5280, 10, 30],
    # # ["zz_07",["07_zz_12197024"], ["BG"],  5125, 15, 47],
    # # ["zz_10",["10_zz_11196413"], ["GAddBGEP"], 5130, 15, 60],
    #
    # ["zz_0710", ["07_zz_12197024","10_zz_11196413"], ["BGEP"], 5130, 20, 50],
    #
    # # ["zz_04",["04_zz_11198022"], ["GEP","BGEP","BG"], 4970, 10, 50],
    # # ["zz_08",["08_zz_12197422"], ["GEP","BGEP","BG"], 4980, 10, 65],
    # # ["zz_12",["12_zz_11196414"], ["GEP","BGEP","BG"], 4980, 10, 80],

    # ["zz_040812", ["04_zz_11198022","08_zz_12197422", "12_zz_11196414"], ["BGEP"], 4975, 10, 50]
    #
    # ["p_05", ["05_p_12197023"], ["G"], 5280, 10, 30],
    # #
    # # ["p_02", ["02_p_11198005"], "cb1L", 1],
    # # ["p_06", ["06_p_12197410"], ["GEP","BGEP"], 5130, 15, 80],
    # # ["p_07", ["07_p_12197400"], "cb1L", 1],
    # #
    # ["p_020607", ["02_p_11198005", "06_p_12197410", "07_p_12197400"], ["GAddBGEP"], 5130, 20, 50],
    # #
    # # ["p_04", ["04_p_11198410"], "cb1L", 2],
    # # ["p_08", ["08_p_12197401"], "cb1L", 2],
    # #
    # ["p_0408", ["04_p_11198410", "08_p_12197401"], ["GEP","BGEP","BG"], 4990, 20, 50],
    # #

    # ["m_03", ["02_p_11198005"], ["G", "DG"], 5130, 20, 50],
    # ["m_04", ["04_p_11198410"], ["BG","GEP","BGEP","G","DG","GAddBGEP"], 4990, 20, 50],


    # ["st_07", ["07_st_12197024"], ["G","DG","BG","GEP","BGEP"],  5280, 10, 30],
    # # # ["st_4", ["04_st_11198022"], "cb1L", 1],
    # # # ["st_8", ["08_st_12197422"], "cb1L", 1],
    # # #
    # ["st_0408", ["04_st_11198022","08_st_12197422"], ["BG","GEP","BGEP"], 5130, 10, 50],

    # ("13_s_13198040","DG"),
    # ("14_s_13198200", "cb1L"),
    # ("15_s_13198400", "cb1L"),
    # ("16_s_13198600", "cb1L"),
    # ["norm7", ["norm7_norm7_12197008"], "DG", 0],
    # ["norm8", ["norm8_norm8_11198007"], "DG", 0],
]

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
