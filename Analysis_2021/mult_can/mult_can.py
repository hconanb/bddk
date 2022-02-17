#!/usr/bin/python
import sys
import os
import ROOT
import pandas as pd
import uproot
import glob
import datetime

def getBestCand(file, spec, year, type):

    treeIn = uproot.open(f"{file}:DecayTreeTuple")

    mydf0 = treeIn.arrays(treeIn.keys(), library="pd")

    df = pd.DataFrame(mydf0[0], columns=treeIn.keys())

    df_test = df.duplicated(subset=['runNumber' ,'eventNumber'])

    if (len(df) == df_test.value_counts()[0]):
        cue = 1
        cc = 0
    else:
        cue = (df_test.value_counts()[0] + df_test.value_counts()[1])/df_test.value_counts()[0]
        cc = df_test.value_counts()[1]


    df_bestCands = df.sort_values('B_dtf_chi2', ascending=False).drop_duplicates(subset=['runNumber' ,'eventNumber']).sort_index().reset_index(level=1, drop=True)

    if type == "Data":
        if "norm" not in spec:
            fileOut = uproot.recreate(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/{year}/final_sample/{spec}.root")
            fileOut["DecayTreeTuple"] = df_bestCands
        else:
            df_T = df_bestCands[(df_bestCands['B_L0HadronDecision_TOS'] == 1) | (df_bestCands['B_L0MuonDecision_TOS'] == 1) | (df_bestCands['B_L0PhotonDecision_TOS'] == 1) | (df_bestCands['B_L0ElectronDecision_TOS'] == 1)]
            df_nTaT = df_bestCands[((df_bestCands['B_L0HadronDecision_TOS'] == 0) & (df_bestCands['B_L0MuonDecision_TOS'] == 0) & (df_bestCands['B_L0PhotonDecision_TOS'] == 0) & (df_bestCands['B_L0ElectronDecision_TOS'] == 0)) & ((df_bestCands['B_L0HadronDecision_TIS'] == 1) | (df_bestCands['B_L0MuonDecision_TIS'] == 1) | (df_bestCands['B_L0PhotonDecision_TIS'] == 1) | (df_bestCands['B_L0ElectronDecision_TIS'] == 1))]
            print(len(df_bestCands.index), len(df_T.index), len(df_nTaT.index),)
            fileOut = uproot.recreate(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/final_sample/{spec}.root")
            fileOut["DecayTreeTuple_T"] = df_T
            fileOut["DecayTreeTuple_nTaT"] = df_nTaT

    if type == "MC":
        # if "norm" not in spec:
        #     fileOut = uproot.recreate(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/final_sample/{spec}.root")
        #     fileOut["DecayTreeTuple"] = df_bestCands
        # else:
        df_T = df_bestCands[(df_bestCands['B_L0HadronDecision_TOS'] == 1) | (df_bestCands['B_L0MuonDecision_TOS'] == 1) | (df_bestCands['B_L0PhotonDecision_TOS'] == 1) | (df_bestCands['B_L0ElectronDecision_TOS'] == 1)]
        df_nTaT = df_bestCands[((df_bestCands['B_L0HadronDecision_TOS'] == 0) & (df_bestCands['B_L0MuonDecision_TOS'] == 0) & (df_bestCands['B_L0PhotonDecision_TOS'] == 0) & (df_bestCands['B_L0ElectronDecision_TOS'] == 0)) & ((df_bestCands['B_L0HadronDecision_TIS'] == 1) | (df_bestCands['B_L0MuonDecision_TIS'] == 1) | (df_bestCands['B_L0PhotonDecision_TIS'] == 1) | (df_bestCands['B_L0ElectronDecision_TIS'] == 1))]
        print(len(df_bestCands.index), len(df_T.index), len(df_nTaT.index),)


        dict_temp = {
                    'Spectrum' : spec,
                    'Year': year,
                    'Canidates/Unique Event': cue,
                    'Canidates Cut': cc,
                    'Final TOS Events': len(df_T.index),
                    'Final TIS Events': len(df_nTaT.index)
                    }

        print(dict_temp)

        # now = datetime.datetime.now()
        # if not os.path.exists(f'mult_can_txt_files/{now.month}_{now.day}/'):
        #     os.makedirs(f'mult_can_txt_files/{now.month}_{now.day}/')
        folder_path = '/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_efficiencies/multi_can_effs'

        df_dict = pd.DataFrame(dict_temp, index=[0])
        outefftxt = f"{folder_path}/{spec}_{year}_can"
        df_dict.to_csv(outefftxt)

        fileOut = uproot.recreate(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/final_sample/{spec}.root")
        fileOut["DecayTreeTuple_T"] = df_T
        fileOut["DecayTreeTuple_nTaT"] = df_nTaT
#"Z_m_p", "P_z_p", "Z_z_z", "M_m_z", "P_z_pst", "Zs_sm_p",

spec_list = ["norm7", "norm8"]
years = ["2016","2017","2018"]
#
for spec in spec_list:
    for year in years:
        files = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/post_d/{spec}.root")
        getBestCand(files[0], spec, year, "Data")

# mc_spec_list = [
#     # "01_Z_m_p_11198006",
#     # "02_Z_m_p_11198400",
#     # "02_P_z_p_11198005",
#     # "04_Z_m_p_11198401",
#     # "04_P_z_p_11198410",
#     # "04_Z_z_z_11198023",
#     # "04_P_z_pst_11198023",
#     # "05_P_z_p_12197023",
#     # "06_P_z_p_12197410",
#     # "07_P_z_p_12197400",
#     # "07_Z_z_z_12197045",
#     # "07_P_z_pst_12197045",
#     # "08_P_z_p_12197401",
#     # "08_Z_z_z_12197423",
#     # "08_P_z_pst_12197423",
#     # "09_Z_z_z_11196019",
#     # "10_Z_z_z_11196413",
#     # "12_Z_z_z_11196414",
#     # "13_Zs_sm_p_13198040",
#     # "14_Zs_sm_p_13198200",
#     # "15_Zs_sm_p_13198400",
#     # "16_Zs_sm_p_13198600",
#     "norm7_norm7_12197008",
#     "norm8_norm8_11198007",
# ]
#
# for file_id in mc_spec_list:
#     for year in ["2016","2017","2018"]:
#         if "Z_m_p" in file_id or "P_z_pst" in file_id or "Zs_sm_p" in file_id or "norm" in file_id:
#             file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/{year}/post_d/{file_id}.root"
#         elif ("Z_z_z" in file_id or "P_z_p" in file_id):
#             file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/{year}/post_veto/{file_id}.root"
#         getBestCand(file_path, file_id, year, "MC")
