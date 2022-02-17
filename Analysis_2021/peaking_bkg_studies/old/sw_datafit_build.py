import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

run_name = "sw_test"
fit_strat = "dp_f_nn"

specs = ["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]
ids_and_shapes = {
                     "Z_m_p" : [["01", "G", 5280, 15],
                               ["02", "BG", 5130, 15],
                               ["04", "G", 4980, 15]],

                     "Z_z_z": [["09", "DG", 5280, 10],
                           ["0710", "BG", 5130, 15],
                           ["040812", "BG", 4975, 20]],

                     "P_z_p": [["05", "G", 5280, 10],
                           ["020607","BGEP", 5130, 15],
                           ["0408", "BG", 4975, 20]],

                     "M_m_z": [["03", "BG", 5130, 15],
                           ["04", "BG", 4975, 20]],

                     "P_z_pst": [["07", "BG", 5280, 10],
                           ["0408","DG", 5130, 15]],

                     "Zs_sm_p" : [["13","G", 5370, 10],
                                 ["14","G", 5225, 15],
                                 ["15","G", 5225, 15],
                                 ["16","BG", 5075, 15]]}

#build fit for data -> Save
#run fit for data -> Save S weighted Data
#plot sweighted data for a peak vs relevant mc

dws = ROOT.RooWorkspace(run_name)
dws.factory(f"B_DTF_M[{bmin},{bmax}]")

for spec in specs:
    get_shapes_bkg(spec, "Exponential", dws)
    for s_list in ids_and_shapes[spec]:
        get_free_shapes(dws, spec, fit_strat, s_list)

for i in range(1, 18):
    dws.factory(f"nny_{i}[500,0,10000]")

dws.factory("SUM::Z_m_p_spectrum_all_fit(nny_1*Z_m_p_01_fit,nny_2*Z_m_p_02_fit,nny_3*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
dws.factory("SUM::Z_z_z_spectrum_all_fit(nny_4*Z_z_z_09_fit,nny_5*Z_z_z_0710_fit,nny_6*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
dws.factory("SUM::P_z_p_spectrum_all_fit(nny_7*P_z_p_05_fit,nny_8*P_z_p_020607_fit,nny_9*P_z_p_0408_fit,n_P_z_P_z_p_bkg[100,0,10000]*P_z_p_spectrum_bkg)")
dws.factory("SUM::M_m_z_spectrum_all_fit(nny_10*M_m_z_03_fit,nny_11*M_m_z_04_fit,n_M_m_z_bkg[100,0,10000]*M_m_z_spectrum_bkg)")
dws.factory("SUM::P_z_pst_spectrum_all_fit(nny_12*P_z_pst_07_fit,nny_13*P_z_pst_0408_fit,n_P_z_st_bkg[100,0,1000]*P_z_pst_spectrum_bkg)")
dws.factory("SUM::s_spectrum_all_fit(nny_14*Zs_sm_p_13_fit,nny_15*Zs_sm_p_14_fit,nny_16*Zs_sm_p_15_fit,nny_17*Zs_sm_p_16_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")
dws.Print()

b_dtf_m = dws.var("B_DTF_M")
data_args = ROOT.RooArgSet(b_dtf_m)
#
# for spec in specs:
    # file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/{spec}*.root")
    # tchain = ROOT.TChain("DecayTreeTuple")
    # for file_name in file_list:
    #     tchain.Add(file_name)
    #
    # data_set = ROOT.RooDataSet(f"{spec}_data", f"{spec}_data", tchain, data_args)
    # dws.Import(data_set)

# vars = dws.allVars()
# for i in vars:
#     if "Z_z_z" in i.GetName() and "mean" in i.GetName():
#         print(i)
#
# dws.writeToFile(f"root_files/{run_name}.root")
