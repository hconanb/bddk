import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021")
from essentials import *

def plot_data(run_name, fit_strat, specs, data_only_flag):

    print(run_name)
    fws_base_plot_file = ROOT.TFile(f"{analysis_path}/data_fit/fit_data_files/fit_{run_name}.root")
    fws = fws_base_plot_file.Get(f"fit_{run_name}")
    fws.Print()

    b_dtf_m = fws.var("B_DTF_M")
    all_cats = fws.cat("all_cats")

    all_data_sets = fws.data("all_data_sets")
    all_fit = fws.pdf("super_fit_Pdf")

    all_cats_Set = ROOT.RooArgSet(all_cats)
    fitresult = fws.obj("fitresult_super_fit_Pdf_all_data_sets")

    bf_string_dict = {
        "1": "B^{0} -> D- D+ K*0",
        "2": "B^{0} -> D*- D+ K*0",
        "3": "B^{0} -> D- D*+ K*0",
        "4": "B^{0} -> D*- D*+ K*0",
        "5": "B^{+} -> #bar{D0} D+ K*0",
        "6": "B^{+} -> #bar{D*0} D+ K*0",
        "7": "B^{+} -> #bar{D0} D*+ K*0",
        "8": "B^{+} -> #bar{D*0} D*+ K*0",
        "9": "B^{0} -> D0 D0 K*0",
        "10": "B^{0} -> D*0 D0 K*0 + B^{0} -> D0 D*0 K*0",
        "12": "B^{0} -> D*0 D*0 K*0",
        "13": "B_{s}^{0} -> D_{s}- D+ K*0",
        "14": "B_{s}^{0} -> D_{s}*- D+ K*0",
        "15": "B_{s}^{0} -> D_{s}- D*+ K*0",
        "16": "B_{s}^{0} -> D_{s}*- D*+ K*0",
    }

    file = open("bf.txt","w")
    for i in range(1,13):
        if i != 11:
            bf = fws.var(f"bf_{i}")
            bf_val = bf.getValV()
            bf_err = bf.getPropagatedError(fitresult)
            bf_print = ufloat(bf_val, bf_err)
            string = bf_string_dict[str(i)]
            file.write(f"{string} : {bf_print:.2E} \n")
    file.close()


    for spec in specs:

        spectrum = f"{spec}_spectrum"
        frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}_spectrum"), ROOT.RooFit.Bins(bbins))
        all_data_sets.plotOn(frame, ROOT.RooFit.Cut(f"all_cats==all_cats::{spec}_spectrum"), ROOT.RooFit.Name("data"))

        if spec == "Z_m_p":
            ylist = ["01","0203","04"]
            title = "D^{-} D^{+} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{+}"

        if spec == "Z_z_z":
            title = "#bar{D^{0}} D^{0} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{0}"
            ylist = ["09", "0710", "040812"]

        if spec == "P_z_p":
            title = "#bar{D^{0}} D^{+} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{+}"
            ylist = ["05", "020607", "0408"]

        if spec == "M_m_z" :
            title = "D^{-} D^{0} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{0}"
            ylist = ["03", "04"]

        if spec == "P_z_pst" :
            title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
            d1 = "D^{-}"
            d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
            ylist = ["07", "0408"]

        if spec == "Zs_sm_p" :
            ylist = ["13","14","15","16"]
            title = "D_{s}^{-} D^{+} K^{*0}"
            d1 = "D_{s}^{-}"
            d2 = "D^{+}"

        if data_only_flag == 0:

            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("data"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_all_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))

            if spec == "Z_m_p":
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                ylist = ["01","0203","04"]
                title = "D^{-} D^{+} K^{*0}"
                d1 = "D^{-}"
                d2 = "D^{+}"

            if spec == "Z_z_z":
                title = "#bar{D^{0}} D^{0} K^{*0}"
                d1 = "#bar{D^{0}}"
                d2 = "D^{0}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_09_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0710_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_040812_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                ylist = ["09", "0710", "040812"]

            if spec == "P_z_p":
                title = "#bar{D^{0}} D^{+} K^{*0}"
                d1 = "#bar{D^{0}}"
                d2 = "D^{+}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_05_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_020607_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0408_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                ylist = ["05", "020607", "0408"]

            if spec == "M_m_z" :
                title = "D^{-} D^{0} K^{*0}"
                d1 = "D^{-}"
                d2 = "D^{0}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_03_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                ylist = ["03", "04"]

            if spec == "P_z_pst" :
                title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
                d1 = "D^{-}"
                d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_07_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0408_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                ylist = ["07", "0408"]

            if spec == "Zs_sm_p" :
                ylist = ["13","14","15","16"]
                title = "D_{s}^{-} D^{+} K^{*0}"
                d1 = "D_{s}^{-}"
                d2 = "D^{+}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_13_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_14_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_15_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_16_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet), ROOT.RooFit.Name("fit3"))


            d_nyield_bkg = fws.obj(f"n_{spec}_bkg")
            d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

        p = ROOT.TCanvas("p1","p1")
        p.cd()

        frame.GetXaxis().SetTitle(f" m({title}) [MeV]")
        frame.Draw()

        if data_only_flag == 0:

            legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
            legend.AddEntry("pdf","Total Fit","lp")
            legend.AddEntry("bkg","Background","l")
            legend.AddEntry("fit0", "Fit 0 Missing Particles","l")
            legend.AddEntry("fit1", "Fit 1 Missing Particles","l")
            if spec != "M_m_z":
                if spec != "P_z_pst":
                    legend.AddEntry("fit2","Fit 2 Missing Particles","l")

            legend.SetTextSize(0.030)
            legend.Draw()

            d_chi2 = frame.chiSquare(f"pdf", f"data")
            if spec != "Z_z_z":
                dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
            if spec == "Z_z_z":
                dtpave = ROOT.TPaveText(0.50, 0.65, 0.70, 0.85, "NB NDC")

            dtpave.SetFillStyle(0)
            dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
            dtpave.AddText(f"Yield BKG: {round(d_nyield_bkg.getValV(),3)} #pm {round(d_nyield_bkg_err,3)}")

            if len(ylist) == 3:
                for y,name in zip(ylist,["Yield Right", "Yield Middle", "Yield Left"]):
                    d_nyield = fws.obj(f"n_{y}_{spec}")
                    d_nyield_err = d_nyield.getPropagatedError(fitresult)
                    dtpave.AddText(f"{name}: {round(d_nyield.getValV(),3)} #pm {round(d_nyield_err,3)}")
            if len(ylist) == 2:
                for y,name in zip(ylist,["Yield Right", "Yield Left"]):
                    d_nyield = fws.obj(f"n_{y}_{spec}")
                    d_nyield_err = d_nyield.getPropagatedError(fitresult)
                    dtpave.AddText(f"{name}: {round(d_nyield.getValV(),3)} #pm {round(d_nyield_err,3)}")
            if len(ylist) == 4:
                for y,name in zip(ylist,["Yield Right", "Yield Middle 1","Yield Middle 2", "Yield Left"]):
                    d_nyield = fws.obj(f"n_{y}_{spec}")
                    d_nyield_err = d_nyield.getPropagatedError(fitresult)
                    dtpave.AddText(f"{name}: {round(d_nyield.getValV(),3)} #pm {round(d_nyield_err,3)}")

            dtpave.Draw()
            frame.addObject(dtpave)

        save_png(p, f"fit_tests", f"{spectrum}_{run_name}", rpflag = 0)

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--run_name')
parser.add_argument('--fit_strat')
parser.add_argument('--specs', nargs='+')
parser.add_argument('--data_only_flag', default=0)

args = parser.parse_args()
run_name = args.run_name
fit_strat = args.fit_strat
specs = args.specs
data_only_flag = args.data_only_flag

plot_data(run_name, fit_strat, specs, data_only_flag)
