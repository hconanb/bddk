import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
run_name = "sw_test"

specs = ["Z_z_z"]#,"P_z_p","M_m_z","P_z_pst","Zs_sm_p"]

dws_base_file = ROOT.TFile(f"root_files/{run_name}.root")
dws = dws_base_file.Get(f"{run_name}")
b_dtf_m = dws.var("B_DTF_M")

for spec in specs:
    presw_data_set = dws.data(f"{spec}_data")
    fit_model = dws.pdf(f"{spec}_spectrum_all_fit")
    gofit = fit_model.fitTo(presw_data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    if spec == "Z_m_p":
        nyield_1 = dws.obj(f"nny_1")
        nyield_2 = dws.obj(f"nny_2")
        nyield_3 = dws.obj(f"nny_3")
        n_bkg = dws.obj(f"n_{spec}_bkg")
    if spec == "Z_z_z":
        nyield_1 = dws.obj(f"nny_4")
        nyield_2 = dws.obj(f"nny_5")
        nyield_3 = dws.obj(f"nny_6")
        n_bkg = dws.obj(f"n_{spec}_bkg")
    if spec == "P_z_p":
        nyield_1 = dws.obj(f"nny_7")
        nyield_2 = dws.obj(f"nny_8")
        nyield_3 = dws.obj(f"nny_9")
        n_bkg = dws.obj(f"n_{spec}_bkg")
    if spec == "M_m_z":
        nyield_1 = dws.obj(f"nny_10")
        nyield_2 = dws.obj(f"nny_11")
        n_bkg = dws.obj(f"n_{spec}_bkg")
    if spec == "P_z_pst":
        nyield_1 = dws.obj(f"nny_12")
        nyield_2 = dws.obj(f"nny_13")
        n_bkg = dws.obj(f"n_{spec}_bkg")
    if spec == "Zs_sm_p":
        nyield_1 = dws.obj(f"nny_14")
        nyield_2 = dws.obj(f"nny_15")
        nyield_3 = dws.obj(f"nny_16")
        nyield_4 = dws.obj(f"nny_17")
        n_bkg = dws.obj(f"n_{spec}_bkg")
    cs = ROOT.TCanvas("cs","cs")
    frame = b_dtf_m.frame()
    presw_data_set.plotOn(frame)
    fit_model.plotOn(frame)
    if spec == "Z_m_p":
        fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
        fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
        fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
    if spec == "Z_z_z":
        fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_09_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
        fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_0710_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
        fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_040812_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))

    # vars = dws.allVars()
    # mean_vars = ROOT.RooArgSet()
    # for i in vars:
    #     if f"{spec}" in i.GetName() and "mean" in i.GetName():
    #         mean_vars.add(i)
    #
    # dtpave = ROOT.TPaveText(0.50, 0.65, 0.70, 0.85, "NB NDC")
    # dtpave.SetFillStyle(0)
    # dtpave.AddText(f"Mean 0: {round(d_chi2, 3)}")
    # dtpave.AddText(f"Mean 1: {round(d_chi2, 3)}")
    # dtpave.AddText(f"Mean 2: {round(d_chi2, 3)}")



    frame.Draw()
    save_png(cs, f"fit_tests", f"{run_name}", rpflag = 0)
    if spec == "Z_m_p" or spec == "P_z_p" or spec == "Z_z_z":
        yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, n_bkg)

    sData = ROOT.RooStats.SPlot("sData","An SPlot", presw_data_set, fit_model, yields)
    # dws.Import(presw_data_set, ROOT.RooFit.Rename(f"swdata_{spec}"))
    # dws.writeToFile(f"root_files/sw_{run_name}.root")

    swnames = sorted([f"{x.GetName()}_sw" for x in yields])
    print(f"weights: {swnames}")
    # create output file
    nf = ROOT.TFile.Open(f"root_files/sw_{run_name}_{spec}.root", "recreate")
    outtreename = "testtree"
    # create directory hierarchy
    nd = nf
    if len(outtreename.split("/")) > 1:
        for d in outtreename.split("/")[:-1]:
            nd = nd.mkdir(d)

    nd.cd()
    # create output TTree
    nt = ROOT.TTree(outtreename.split("/")[-1], outtreename.split("/")[-1])
    # declare output branches
    swvals = [array("f", [0]) for x in swnames]
    for val, nm in zip(swvals, swnames):
        nt.Branch(nm, val, f"{nm}/F")
    # loop data
    for i in range(presw_data_set.numEntries()):
        # get vars
        swvars = sorted(
            [x for x in presw_data_set.get(i) if x.GetName() in swnames],
            key=lambda x: x.GetName(),
        )
        assert [x.GetName() for x in swvars] == swnames  # check sorting worked
        # set values
        for val, var in zip(swvals, swvars):
            val[0] = var.getVal()
        # fill values
        nt.Fill()
    nt.Write()
    nf.Close()
# data_set_sw = dws.data("data_sw")
# b_dtf_m = dws.var("B_DTF_M")
# sw_var = dws.var("nny_1_sw")
# data_args = ROOT.RooArgSet(b_dtf_m, sw_var)
#
# dataset_test = ROOT.RooDataSet(data_set_sw.GetName(), data_set_sw.GetTitle(), data_set_sw, data_args, "", "nny_1_sw")
#
# cs2 = ROOT.TCanvas("cs2","cs2")
# frame_sigsw = b_dtf_m.frame()
# dataset_test.plotOn(frame_sigsw)
# frame_sigsw.Draw()
# save_png(cs2, f"fit_tests", f"{run_name}_sw", rpflag = 0)
# # cs2 = ROOT.TCanvas("cs2","cs2")
# # frame_test = b_dtf_m.frame()
# # data_set.plotOn(frame_test)
# # fit_model.plotOn(frame_test)
# # all_fit.plotOn(frame_test, ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
# # all_fit.plotOn(frame_test, ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
# # all_fit.plotOn(frame_test, ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
# # frame_test.Draw()
# # save_png(cs2, f"fit_tests", f"{run_name}_sw", rpflag = 0)
#
#
# dws.Print()
