import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021")
from essentials import *

def run_data_fit(run_name, gc_onflag, specs):

    lp = 0
    if "Z_m_p" in specs:
        lp = lp + 4
    if "Z_z_z" in specs:
        lp = lp + 6
    if "P_z_p" in specs:
        lp = lp + 6
    if "M_m_z" in specs:
        lp = lp + 2
    if "P_z_pst" in specs:
        lp = lp + 3
    if "Zs_sm_p" in specs:
        lp = lp + 4

    dws_base_file = ROOT.TFile(f"{analysis_path}/data_fit/base_data_files/{run_name}.root")
    dws = dws_base_file.Get(f"{run_name}")

    if gc_onflag == 1:
        glist = ROOT.RooArgSet()
        for i in range(0, lp):
            gvar = dws.pdf(f"{run_name}_nu{i}_g")
            glist.add(gvar)

    all_data_sets = dws.data("all_data_sets")
    all_cats = dws.cat("all_cats")
    for i in all_cats:
        print(i)
    b_dtf_m = dws.var("B_DTF_M")
    all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)

    dws.Print()

    if "Z_m_p" in specs:
        z_model = dws.pdf("Z_m_p_spectrum_all_fit")
        all_fit.addPdf(z_model, "Z_m_p_spectrum")
    if "Z_z_z" in specs:
        zz_model = dws.pdf("Z_z_z_spectrum_all_fit")
        all_fit.addPdf(zz_model, "Z_z_z_spectrum")
    if "P_z_p" in specs:
        p_model = dws.pdf("P_z_p_spectrum_all_fit")
        all_fit.addPdf(p_model, "P_z_p_spectrum")
    if "M_m_z" in specs:
        m_model = dws.pdf("M_m_z_spectrum_all_fit")
        all_fit.addPdf(m_model, "M_m_z_spectrum")
    if "P_z_pst" in specs:
        st_model = dws.pdf("P_z_pst_spectrum_all_fit")
        all_fit.addPdf(st_model, "P_z_pst_spectrum")
    if "Zs_sm_p" in specs:
        s_model = dws.pdf("Zs_sm_p_spectrum_all_fit")
        all_fit.addPdf(s_model, "Zs_sm_p_spectrum")
    if gc_onflag == 1:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(glist), ROOT.RooFit.Save())
    if gc_onflag == 0:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    fitws = ROOT.RooWorkspace(f"fit_{run_name}")
    fitws.Import(all_data_sets)
    fitws.Import(all_cats)
    fitws.Import(all_fit)
    fitws.Import(nall_fit)
    fitws.Print()
    fitws.writeToFile(f"{analysis_path}/data_fit/fit_data_files/fit_{run_name}.root")

    """Determine s-weights from fit
    arguments:
    outfilename -- name of .root file to create with `outtreename`
    outtreename -- name of TTree with s-weights to save in `outfilename`
    data -- RooDataSet to which `model` was fitted
    model -- fitted RooAbsPdf
    yields -- RooArgList of RooRealVars extracted from fitting `model` to `data`
    """
    fit_test = fitws.pdf("Z_m_p_spectrum_all_fit")
    data_test = fitws.data("all_data_sets")
    # data_test2 = data_test.Get("Z_m_p_spectrum")
    nyield_1 = fitws.obj(f"nny_1")
    nyield_2 = fitws.obj(f"nny_2")
    nyield_3 = fitws.obj(f"nny_3")
    n_Z_m_p_bkg = fitws.obj(f"n_Z_m_p_bkg")

    yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, n_Z_m_p_bkg)

    # print(f"using data '{data_test.GetName()}'")
    # # print(f"using model '{fit_test.GetName()}'")
    # # print(f"using yields '{[x.GetName() for x in yields]}'")

    sData = ROOT.RooStats.SPlot("sData","An SPlot", data_test, fit_test, yields)
    # fitws.Import(data_test, Rename("dataWithSWeights"))
    # fitws.Import(fit_test, ROOT.RooCmdArg.Rename("modelWithSWeights"))

    cs = ROOT.TCanvas("cs","cs")
    frame_test = b_dtf_m.frame()
    data_test.plotOn(frame_test,ROOT.RooFit.Name("Hist"))
    fit_test.plotOn(frame_test,ROOT.RooFit.Name("curvetot"))
    frame_test.Draw()
    save_png(cs, f"fit_tests", f"{run_name}", rpflag = 0)

    # model->plotOn(frame,Components(*mBModel),LineStyle(kDashed), LineColor(kRed)) ;
    # model->plotOn(frame,Components(*bkg),LineStyle(kDashed),LineColor(kGreen)) ;




parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--run_name')
parser.add_argument('--gc_onflag', type=int)
parser.add_argument('--specs', nargs='+')

args = parser.parse_args()
run_name = args.run_name
gc_onflag = args.gc_onflag
specs  = args.specs

run_data_fit(run_name, gc_onflag, specs)
