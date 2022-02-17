import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

specs = ["Z_m_p"] #"Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]

run_name = f"swbkg"

dws = ROOT.RooWorkspace(run_name)
dws.factory(f"B_DTF_M[{bmin},{bmax}]")

ids_and_shapes = {
                    "Z_m_p" : [["01", "G", 5280, 15],
                               ["02", "BG", 5130, 15],
                               ["04", "G", 4980, 15]]}

ids_and_shapes

for spec in specs:
    get_shapes_bkg(spec, "Exponential", dws)
    for s_list in ids_and_shapes[spec]:
        get_free_shapes(dws, spec, fit_strat, s_list)

for i in range(1, 4):
    dws.factory(f"nny_{i}[500,0,10000]")

spec = "Z_m_p"

dws.factory("SUM::Z_m_p_spectrum_all_fit(nny_1*Z_m_p_01_fit,nny_2*Z_m_p_02_fit,nny_3*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
# dws.factory("SUM::Z_z_z_spectrum_all_fit(nny_4*Z_z_z_09_fit,nny_5*Z_z_z_0710_fit,nny_6*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
# dws.factory("SUM::P_z_p_spectrum_all_fit(nny_7*P_z_p_05_fit,nny_8*P_z_p_020607_fit,nny_9*P_z_p_0408_fit,n_P_z_P_z_p_bkg[100,0,10000]*P_z_p_spectrum_bkg)")
# dws.factory("SUM::M_m_z_spectrum_all_fit(nny_10*M_m_z_03_fit,nny_11*M_m_z_04_fit,n_M_m_z_bkg[100,0,10000]*M_m_z_spectrum_bkg)")
# dws.factory("SUM::P_z_pst_spectrum_all_fit(nny_12*P_z_pst_07_fit,nny_13*P_z_pst_0408_fit,n_P_z_st_bkg[100,0,1000]*P_z_pst_spectrum_bkg)")
# dws.factory("SUM::s_spectrum_all_fit(nny_14*s0_fit,nny_15*sL1_fit,nny_16*sR1_fit,nny_17*s2_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")
dws.Print()

file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/Z_m_p*.root")
tchain = ROOT.TChain("DecayTreeTuple")
for file_name in file_list:
    tchain.Add(file_name)

b_dtf_m = dws.var("B_DTF_M")
data_args = ROOT.RooArgSet(b_dtf_m)
data_set = ROOT.RooDataSet(f"Z_m_p_data", f"Z_m_p_data", tchain, data_args)

all_fit = dws.pdf("Z_m_p_spectrum_all_fit")
nall_fit = all_fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

cs = ROOT.TCanvas("cs","cs")
frame = b_dtf_m.frame()
data_set.plotOn(frame)
all_fit.plotOn(frame)
all_fit.plotOn(frame, ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
all_fit.plotOn(frame, ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
all_fit.plotOn(frame, ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
frame.Draw()
save_png(cs, f"fit_tests", f"{run_name}", rpflag = 0)

nyield_1 = dws.obj(f"nny_1")
nyield_2 = dws.obj(f"nny_2")
nyield_3 = dws.obj(f"nny_3")
n_Z_m_p_bkg = dws.obj(f"n_Z_m_p_bkg")

yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, n_Z_m_p_bkg)

    # print(f"using data '{data_test.GetName()}'")
    # # print(f"using model '{fit_test.GetName()}'")
    # # print(f"using yields '{[x.GetName() for x in yields]}'")

sData = ROOT.RooStats.SPlot("sData","An SPlot", data_set, all_fit, yields)
dws.Import(data_set, ROOT.RooFit.Rename("data_sw"))

data_set_sw = dws.data("data_sw")
b_dtf_m = dws.var("B_DTF_M")
sw_var = dws.var("nny_1_sw")
data_args = ROOT.RooArgSet(b_dtf_m, sw_var)

dataset_test = ROOT.RooDataSet(data_set_sw.GetName(), data_set_sw.GetTitle(), data_set_sw, data_args, "", "nny_1_sw")

cs2 = ROOT.TCanvas("cs2","cs2")
frame_sigsw = b_dtf_m.frame()
dataset_test.plotOn(frame_sigsw)
frame_sigsw.Draw()
save_png(cs2, f"fit_tests", f"{run_name}_sw", rpflag = 0)
# cs2 = ROOT.TCanvas("cs2","cs2")
# frame_test = b_dtf_m.frame()
# data_set.plotOn(frame_test)
# all_fit.plotOn(frame_test)
# all_fit.plotOn(frame_test, ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
# all_fit.plotOn(frame_test, ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
# all_fit.plotOn(frame_test, ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
# frame_test.Draw()
# save_png(cs2, f"fit_tests", f"{run_name}_sw", rpflag = 0)


dws.Print()
