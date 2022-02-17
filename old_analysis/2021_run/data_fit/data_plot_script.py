import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *


def plot_data(run_name, fit_strat, specs):

    fws_base_plot_file = ROOT.TFile(f"{analysis_path}/data_fit/fit_data_files/fit_{run_name}.root")
    fws = fws_base_plot_file.Get(f"fit_{run_name}")

    b_dtf_m = fws.var("B_DTF_M")
    all_cats = fws.cat("all_cats")
    all_data_sets = fws.data("all_data_sets")
    all_fit = fws.pdf("super_fit_Pdf")
    all_cats_Set = ROOT.RooArgSet(all_cats)
    fitresult = fws.obj("fitresult_super_fit_Pdf_all_data_sets")

    for spec in specs:

        spectrum = f"{spec}_spectrum"
        frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}_spectrum"), ROOT.RooFit.Bins(bbins))

        all_data_sets.plotOn(frame, ROOT.RooFit.Cut(f"all_cats==all_cats::{spec}_spectrum"), ROOT.RooFit.Name("data"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_all_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))

        if spec == "z" :
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
            ylist = ["01","02","04"]
            title = "D^{-} D^{+} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{+}"
            ylist = ["nny_1", "nny_2", "nny_3"]
        if spec == "zz" :
            title = "#bar{D^{0}} D^{0} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{0}"
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_09_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0710_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_040812_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
            ylist = ["nny_4", "nny_5", "nny_6"]

        if spec == "p" :
            title = "#bar{D^{0}} D^{+} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{+}"
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_05_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_020607_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0408_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
            ylist = ["nny_7", "nny_8", "nny_9"]

        if spec == "m" :
            title = "D^{-} D^{0} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{0}"
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_03_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            ylist = ["nny_10", "nny_11"]

        if spec == "st" :
            title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
            d1 = "D^{-}"
            d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_07_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0408_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            ylist = ["nny_12", "nny_13"]

        if spec == "s":
            title = "D_{s}^{-} D^{+} K^{*0}"
            d1 = "D_{s}^{-}"
            d2 = "D^{+} "
        d3 = "K^{*0}"


        d_nyield_bkg = fws.obj(f"n_{spec}_bkg")
        d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

        d_nyield_0 = fws.obj(f"{ylist[0]}")
        d_nyield_0_err = d_nyield_0.getPropagatedError(fitresult)

        d_nyield_1 = fws.obj(f"{ylist[1]}")
        d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)

        if spec != "m":
            if spec != "st":
                d_nyield_2 = fws.obj(f"{ylist[2]}")
                d_nyield_2_err = d_nyield_2.getPropagatedError(fitresult)



        #     d_nbkg = fws.var(f"n_z_bkg")
        # if spec == "zz":

        #
        #     d_nyield_0 = fws.obj(f"n_09_zz")
        #     d_nyield_0_err = d_nyield_0.getPropagatedError(fitresult)
        #
        #     d_nyield_1 = fws.obj(f"n_07_zz")
        #     d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)
        #
        #     d_nyield_2 = fws.obj(f"n_10_zz")
        #     d_nyield_2_err = d_nyield_2.getPropagatedError(fitresult)
        #
        #     d_nyield_3 = fws.obj(f"n_04_zz")
        #     d_nyield_3_err = d_nyield_3.getPropagatedError(fitresult)
        #
        #     d_nyield_4 = fws.obj(f"n_08_zz")
        #     d_nyield_4_err = d_nyield_4.getPropagatedError(fitresult)
        #
        #     d_nyield_5 = fws.obj(f"n_12_zz")
        #     d_nyield_5_err = d_nyield_5.getPropagatedError(fitresult)
        #
        #     d_nbkg = fws.var(f"n_zz_bkg")
        # if spec == "p":

        #
        #     d_nyield_0 = fws.obj(f"n_05_p")
        #     d_nyield_0_err = d_nyield_0.getPropagatedError(fitresult)
        #
        #     d_nyield_1 = fws.obj(f"n_02_p")
        #     d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)
        #
        #     d_nyield_2 = fws.obj(f"n_06_p")
        #     d_nyield_2_err = d_nyield_2.getPropagatedError(fitresult)
        #
        #     d_nyield_3 = fws.obj(f"n_07_p")
        #     d_nyield_3_err = d_nyield_3.getPropagatedError(fitresult)
        #
        #     d_nyield_4 = fws.obj(f"n_04_p")
        #     d_nyield_4_err = d_nyield_4.getPropagatedError(fitresult)
        #
        #     d_nyield_5 = fws.obj(f"n_08_p")
        #     d_nyield_5_err = d_nyield_5.getPropagatedError(fitresult)
        #
        #     d_nbkg = fws.var(f"n_p_bkg")
        # if spec == "m":

        #
        #     d_nyield_0 = fws.obj(f"n_03_m")
        #     d_nyield_0_err = d_nyield_0.getPropagatedError(fitresult)
        #
        #     d_nyield_1 = fws.obj(f"n_04_m")
        #     d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)
        #
        #     d_nbkg = fws.var(f"n_m_bkg")
        # if spec == "st":

        #     d_nyield_0 = fws.obj(f"n_07_st")
        #     d_nyield_0_err = d_nyield_0.getPropagatedError(fitresult)
        #
        #     d_nyield_1 = fws.obj(f"n_04_st")
        #     d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)
        #
        #     d_nyield_2 = fws.obj(f"n_08_st")
        #     d_nyield_2_err = d_nyield_2.getPropagatedError(fitresult)
        #
        #     d_nbkg = fws.var(f"n_st_bkg")

        # if spec == "s":

        p = ROOT.TCanvas("p1","p1")
        p.cd()

        frame.GetXaxis().SetTitle(f" m({title}) [MeV]")
        frame.Draw()

        legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
        #legend.SetHeader(title,"C")
        legend.AddEntry("pdf","Total Fit","lp")
        legend.AddEntry("bkg","Background","l")
        legend.AddEntry("fit0", "Fit 0 Missing Particles","l")
        legend.AddEntry("fit1", "Fit 1 Missing Particles","l")
        if spec != "m":
            if spec != "st":
                legend.AddEntry("fit2","Fit 2 Missing Particles","l")
        legend.SetTextSize(0.030)
        legend.Draw()

        d_chi2 = frame.chiSquare(f"pdf", f"data")
        if spec != "zz":
            dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
        if spec == "zz":
            dtpave = ROOT.TPaveText(0.50, 0.65, 0.70, 0.85, "NB NDC")

        dtpave.SetFillStyle(0)
        dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
        dtpave.AddText(f"Yield Right: {round(d_nyield_0.getValV(),3)} #pm {round(d_nyield_0_err,3)}")
        dtpave.AddText(f"Yield Middle: {round(d_nyield_1.getValV(),3)} #pm {round(d_nyield_1_err,3)}")
        if spec != "m":
            if spec != "st":
                dtpave.AddText(f"Yield Left: {round(d_nyield_2.getValV(),3)}#pm {round(d_nyield_2_err,3)}")
        dtpave.AddText(f"Yield BKG: {round(d_nyield_bkg.getValV(),3)} #pm {round(d_nyield_bkg_err,3)}")
        dtpave.Draw()
        frame.addObject(dtpave)

        save_png(p, f"fit_tests", f"{spectrum}_{run_name}", rpflag = 0)

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--run_name')
parser.add_argument('--fit_strat')
parser.add_argument('--specs', nargs='+')

args = parser.parse_args()
run_name = args.run_name
fit_strat = args.fit_strat
specs  = args.specs

plot_data(run_name, fit_strat, specs)
