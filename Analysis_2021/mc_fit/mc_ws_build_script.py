import ROOT as ROOT
from Analysis_2021.essentials import get_free_shapes, save_png
from Analysis_2021 import rootutils
import glob as glob
analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"

def grab_file_list(id_list):
    file_list = []
    for file_id in id_list:
        file_list = file_list + glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/*/final_sample/{file_id}.root")
    return file_list

def build_mc_ws(spec, id_list, shape_list, mean_guess, mean_window, fit_window, mc_truth_flag = False):

    for shape in shape_list:

        mcws = ROOT.RooWorkspace(f"MC_{spec}")
        b_min = mean_guess - fit_window
        b_max = mean_guess + fit_window

        mcws.factory(f"B_DTF_M[{b_min}, {b_max}]")

        get_free_shapes(mcws, spec, "MC", shape, mean_guess, mean_window)

        b_dtf_m = mcws.var("B_DTF_M")

        tc_T = ROOT.TChain(f"DecayTreeTuple_T")
        tc_nTaT = ROOT.TChain(f"DecayTreeTuple_nTaT")

        fl = grab_file_list(id_list)

        for file_name in fl:
            tc_T.Add(file_name)
            tc_nTaT.Add(file_name)

        data_args = ROOT.RooArgSet(b_dtf_m)
        data_set_T = ROOT.RooDataSet(f"{spec}_events_T", f"{spec}_events_T", tc_T, data_args)
        data_set_nTaT = ROOT.RooDataSet(f"{spec}_events_nTaT", f"{spec}_events_nTat", tc_nTaT, data_args)

        data_set_T.append(data_set_nTaT)
        data_set_T.SetNameTitle(f"{spec}_events", f"{spec}_events")
        mcws.Import(data_set_T)

        base_output_file = f"{analysis_path}/mc_fit/base_mc_files/{spec}_{shape}.root"
        mcws.writeToFile(base_output_file)
        print(f"Saved: {base_output_file}")
        mcws.Print()

def build_mc_pizy_ws(spec, piz_spec_list, y_spec_list, mean_start, mean_window, fit_window):

        mcws = ROOT.RooWorkspace(f"MC_{spec}")
        b_min = mean_start - fit_window
        b_max = mean_start + fit_window
        b_low = mean_start - mean_window
        b_high = mean_start + mean_window

        mcws.factory(f"B_DTF_M[{b_min}, {b_max}]")
        b_dtf_m = mcws.var("B_DTF_M")
        data_args = ROOT.RooArgSet(b_dtf_m)
        data_set_list = []
        all_cats = ROOT.RooCategory("all_cats", "all_cats")

        for dspec, id_list in zip(["piz_spectrum", "y_spectrum"], [piz_spec_list, y_spec_list]):
            all_cats.defineType(dspec)
            tc_T = ROOT.TChain(f"DecayTreeTuple_T")
            tc_nTaT = ROOT.TChain(f"DecayTreeTuple_nTaT")
            fl = grab_file_list(id_list)
            for file_name in fl:
                tc_T.Add(file_name)
                tc_nTaT.Add(file_name)
            data_set_T = ROOT.RooDataSet(f"{dspec}_events_T", f"{dspec}_events_T", tc_T, data_args)
            data_set_nTaT = ROOT.RooDataSet(f"{dspec}_events_nTaT", f"{dspec}_events_nTat", tc_nTaT, data_args)
            data_set_T.append(data_set_nTaT)
            data_set_T.SetNameTitle(f"{dspec}_events", f"{dspec}_events")
            data_set_list.append(data_set_T)

        mc_data_sets = ROOT.RooDataSet(
            "mc_data_sets",
            "mc_data_sets",
            data_args,
            ROOT.RooFit.Index(all_cats),
            ROOT.RooFit.Import("piz_spectrum", data_set_list[0]),
            ROOT.RooFit.Import("y_spectrum", data_set_list[1]),
        )
        # mcws.factory(f"BifurGauss:{spec}_y_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}],width_L_{spec}_a[10,0.01,40.0],width_R_{spec}_a[20,0.01,40.0])")
        mcws.factory(f"RooBifurGaussExp:{spec}_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}],width_L_{spec}_a[10,0.01,40.0],width_R_{spec}_a[20,0.01,40.0],alpha_1_{spec}_a[2,0.001,10.0],alpha_2_{spec}_a[3,0.001,10.0])")
        mcws.factory(f"Gaussian::{spec}_fit_b(B_DTF_M, mean_{spec}_a, width_{spec}_b[5.0,0.1,50.0])")
        # mcws.factory(f"RooGaussExp:{spec}_y_fit_b(B_DTF_M, mean_{spec}_a, width_{spec}_b[5.0,0.1,40.0], alpha_{spec}_b[0.1, 0, 10])")
        # mcws.factory(f"Gaussian::{spec}_y_fit_b(B_DTF_M, mean_{spec}_a, width_{spec}_b[5.0,0.1,30.0])")
        mcws.factory(f"SUM:{spec}_fit({spec}_a_frac[0.647, 0.2, 0.8]*{spec}_fit_a, {spec}_fit_b)")


        # mcws.factory(f"RooBifurGaussExp:{spec}_piz_fit(B_DTF_M, mean_{spec}_c[{mean_start},{b_low},{b_high}],width_L_{spec}_c[10,0.01,30.0],width_R_{spec}_c[20,0.01,30.0],alpha_1_{spec}_c[2,0.001,10.0],alpha_2_{spec}_c[3,0.001,10.0])")
        mcws.factory(f"RooBifurGaussExp:{spec}_piz_fit(B_DTF_M, mean_{spec}_c[{mean_start},{b_low},{b_high}], width_L_{spec}_a, width_R_{spec}_a,alpha_1_{spec}_a,alpha_2_{spec}_a)")
        # mcws.factory(f"BifurGauss:{spec}_piz_fit(B_DTF_M, mean_{spec}_c[{mean_start},{b_low},{b_high}],width_L_{spec}_a,width_R_{spec}_a)")

        # mcws.factory(f"RooGaussExp:{spec}_y_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}], width_{spec}_a[5.0,0.1,40.0], alpha_{spec}_a[1, 0, 10])")
        # mcws.factory(f"Gaussian::{spec}_y_fit_b(B_DTF_M, mean_{spec}_b[{mean_start},{b_low},{b_high}], width_{spec}_b[5.0,0.1,50.0])")
        # mcws.factory(f"SUM:{spec}_y_fit({spec}_a_frac[0.647]*{spec}_y_fit_a, {spec}_y_fit_b)")

        # mcws.factory(f"RooGaussExp:{spec}_piz_fit(B_DTF_M, mean_{spec}_c[{mean_start},{b_low},{b_high}], width_{spec}_a, alpha_{spec}_a)")

        all_mc_fit = ROOT.RooSimultaneous("mc_super_fit_Pdf", "mc_super_fit_Pdf", all_cats)

        piz_model = mcws.pdf(f"{spec}_piz_fit")
        all_mc_fit.addPdf(piz_model, "piz_spectrum")

        y_model = mcws.pdf(f"{spec}_fit")
        all_mc_fit.addPdf(y_model, "y_spectrum")

        fit = all_mc_fit.fitTo(mc_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

        all_cats_Set = ROOT.RooArgSet(all_cats)

        for spectrum, comp in zip(["piz_spectrum", "y_spectrum"], [f"{spec}_piz_fit", f"{spec}_fit"]):
            frame = b_dtf_m.frame(ROOT.RooFit.Title(spectrum))
            mc_data_sets.plotOn(frame, ROOT.RooFit.Cut(f"all_cats==all_cats::{spectrum}"), ROOT.RooFit.Name("datda"))
            all_mc_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(comp), ROOT.RooFit.ProjWData(all_cats_Set, mc_data_sets ), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("data"))
            if spectrum == "y_spectrum":
                all_mc_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_fit_a"), ROOT.RooFit.ProjWData(all_cats_Set, mc_data_sets ), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("data_2"))
                all_mc_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_fit_b"), ROOT.RooFit.ProjWData(all_cats_Set, mc_data_sets ), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("data_3"))
            p = ROOT.TCanvas("p1","p1")
            p.cd()
            frame.Draw()
            save_png(p, f"comp_tests", f"{comp}_mfinal", rpflag = 0)

        fitws_y = ROOT.RooWorkspace(f"fit_ws")

        # vars = mcws.allVars()
        # for i in vars:
        #     if i.GetName() in [f"mean_{spec}_a",f"width_L_{spec}_a",f"width_R_{spec}_a",f"alpha_1_{spec}_a",f"alpha_2_{spec}_a"]:
        #         temp = mcws.var(i.GetName())
        #         temp.setConstant(True)
        #         print(f"{temp} is constant")

        y_real_model = mcws.pdf(f"{spec}_fit")
        fitws_y.Import(y_real_model)
        fitws_y.writeToFile(f"{analysis_path}/mc_fit/fit_mc_files/{spec}_fit_pizy.root")

def run_mc_fit(spec, shape_list):

    for shape in shape_list:
        base_input_file = f"{analysis_path}/mc_fit/base_mc_files/{spec}_{shape}.root"
        mcws_base_file = ROOT.TFile(base_input_file)

        mcws = mcws_base_file.Get(f"MC_{spec}")
        b_dtf_m = mcws.var("B_DTF_M")
        data_set = mcws.data(f"{spec}_events")
        fit = mcws.pdf(f"{spec}_fit")
        fit_pdf = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
        fitws = ROOT.RooWorkspace(f"fit_ws")
        fitws.Import(data_set)
        fitws.Import(fit_pdf)
        fitws.Import(fit)
        fitws.Print()
        base_output_file = f"{analysis_path}/mc_fit/fit_mc_files/{spec}_fit_{shape}.root"
        fitws.writeToFile(base_output_file)

def plot_mc(spec, shape_list, pullflag = 1):

    for shape in shape_list:

        base_input_file = f"{analysis_path}/mc_fit/fit_mc_files/{spec}_fit_{shape}.root"
        fws_base_plot_file = ROOT.TFile(base_input_file)
        fws = fws_base_plot_file.Get(f"fit_ws")
        b_dtf_m = fws.var("B_DTF_M")


        #     sw_data_file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/data_fit/sw_files/sw_{data_name}.root")
        #     sw_data_tchain = ROOT.TChain("SW_tree")
        #     for file_name in sw_data_file_list:
        #         print(file_name)
        #         sw_data_tchain.Add(file_name)
        #
        #     for file_name in data_file_list:
        #         tchain.Add(file_name)
        #
        #     tchain.AddFriend(sw_data_tchain, "SW_tree")
        #     print(tchain.GetEntries())
        #     # tchain.SetAlias("nn1", "SW_tree.nny_1_sw")
        #     # print(spec)
        #
        #     if "_01_" in spec:
        #         fws.factory("nny_1_sw[-10, 10]")
        #         sw_var = fws.var("nny_1_sw")
        #
        #     real_data_set = ROOT.RooDataSet("real_data", "real_data", tchain, ROOT.RooArgSet(b_dtf_m, sw_var), cuts = "", wgtVarName = "nny_1_sw")
        # #, ROOT.RooFit.Bins(p_bbins)
        if "P_z_pst" in spec:
            bbb = 50
        else:
            bbb = 100
        frame = b_dtf_m.frame(ROOT.RooFit.Title(spec), ROOT.RooFit.Bins(bbb))

        spectrum = frame.GetTitle()
        sspectrum = spectrum.split("_")[0]

        data_set = fws.data(f"{spec}_events")
        data_set.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.MarkerColor(ROOT.kBlue))

        fit_pdf = fws.pdf(f"{spec}_fit")
        fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("total_fit") )

        if shape in ["DG", "GAddBGEP","GEPAddBG","GEPAddBGEP"] or "_fr" in shape:
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fa") )
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fb") )

        d_chi2 = frame.chiSquare("total_fit", "data")
        print(d_chi2)

        if spec == "Z_m_p_01":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow D^{-} D^{+} K^{*0}"
        if spec == "Z_m_p_02" or spec == "Z_m_p_0203":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow D^{*-} D^{+} K^{*0} / D^{-} D^{*+} K^{*0}"
        if spec == "Z_m_p_04":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow D^{*-} D^{*+} K^{*0}"
        if spec == "Z_z_z_04":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) (D^{*+} #rightarrow D^{0} #pi^{+}) K^{*0}"
        if spec == "Z_z_z_07" or spec == "Z_z_z_0710" or spec == "Z_z_z_040812":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow #bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi^{+}) K^{*0}"
        if spec == "Z_z_z_08":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow (#bar{D^{0}} #rightarrow #bar{D^{0}} #pi^{0}/#gamma)) (D^{*+} #rightarrow D^{0} #pi^{+}) K^{*0}"
        if spec == "Z_z_z_09":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow #bar{D^{0}} D^{0} K^{*0}"
        if spec == "Z_z_z_10":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow(#bar{D^{*0}} #rightarrow #bar{D^{0}} #pi^{0}/#gamma) D^{0} K^{*0}"
        if spec == "Z_z_z_12":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow (#bar{D^{*0}} #rightarrow #bar{D^{0}} #pi^{0}/#gamma) (D^{*0} #rightarrow D^{0} #pi^{0}/#gamma) K^{*0}"
        if spec == "P_z_p_02" or spec == "P_z_p_0408" or spec == "P_z_p_020607":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) D^{+} K^{*0}"
        if spec == "P_z_p_04":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) (D^{*+} #rightarrow D^{+} #pi^{0}/#gamma) K^{*0}"
        if spec == "P_z_p_05":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow #bar{D^{0}} D^{+} K^{*0}"
        if spec == "P_z_p_06":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow (#bar{D^{*0}} #rightarrow #bar{D^{0}} #pi^{0}/#gamma) D^{+} K^{*0}"
        if spec == "P_z_p_07":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow #bar{D^{0}} (D^{*+} #rightarrow D^{+} #pi^{0}/#gamma) K^{*0}"
        if spec == "P_z_p_08":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow (#bar{D^{*0}} #rightarrow #bar{D^{0}} #pi^{0}/#gamma) (D^{*+} #rightarrow D^{+} #pi^{0}/#gamma) K^{*0}"
        if spec == "M_m_z_03":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow D^{-} (D^{*+} #rightarrow D^{0} #pi^{+}) K^{*0}"
        if spec == "M_m_z_04":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow (D^{*-} #rightarrow D^{-} #pi^{0}/#gamma) (D^{*+} #rightarrow D^{0} #pi^{+}) K^{*0}"
        if spec == "P_z_pst_04":
            bsp = "B^{0}"
            le = "B^{0} #rightarrow D^{*-} D^{*+} K^{*0}"
        if spec == "P_z_pst_07" or spec == "P_z_pst_0408":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow #bar{D^{0}} (D^{*+} #rightarrow D^{+} #pi^{0}/#gamma) K^{*0}"
        if spec == "P_z_pst_08":
            bsp = "B^{+}"
            le = "B^{+} #rightarrow #bar{D^{*0}} (D^{*+} #rightarrow D^{+} #pi^{0}/#gamma) K^{*0}"
        if "norm7" in spec:
            bsp = "B^{+}"
            le = "B^{+} #rightarrow #bar{D^{0}} (D^{0} #rightarrow K#pi#pi#pi) K+"
        if "norm8" in spec:
            bsp = "B^{+}"
            le = "B^{0} #rightarrow D^{-} (D^{0} #rightarrow K#pi#pi#pi) K+"

        dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
        dtpave.SetFillStyle(0)
        dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
        dtpave.Draw()
        frame.addObject(dtpave)

        if pullflag == 1:

            p = rootutils.residualPlot()
            p.pt.cd()
            xaxis = frame.GetXaxis()
            xaxis.SetTickLength(0)
            xaxis.SetNdivisions(0)
            xaxis.SetLabelSize(0)
            # frame.addObject(txt)
            frame.Draw()

            legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
            #legend.SetHeader(title,"C")
            legend.AddEntry("total_fit","total_fit","l")
            if shape in ["DG", "GAddBGEP","GEPAddBG","GEPAddBGEP"] or "_fr" in shape:
                legend.AddEntry("fa","fa","l")
                legend.AddEntry("fb","fb","l")

            legend.AddEntry("data","data","lep")
            legend.SetTextSize(0.050)
            legend.Draw()
            # frame.addObject(legend)
            p.pb.cd()

            hpull = frame.pullHist("data","total_fit")
            pull = fws.var("B_DTF_M").frame();
            pull.addPlotable(hpull, "P");
            pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
            pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()*p.blabelratio)
            pull.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Y")*(1+p.padratio)/(2*p.padratio))
            pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
            # print(ROOT.gStyle.GetTitleSize()*p.blabelratio)
            # pull.GetXaxis().SetTitleSize(0.1)
            pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
            pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset()/p.blabelratio)
            pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength()*p.blabelratio)
            pull.GetYaxis().SetTitle("Pull")
            pull.GetYaxis().SetNdivisions(202,False)
            pull.GetYaxis().SetRangeUser(-3,3)
            pull.GetYaxis().CenterTitle()
            pull.Draw("AP")
            pull.GetXaxis().SetTitle(f"{le} [MeV]")

        if pullflag == 0:

            p = ROOT.TCanvas("p1","p1")
            p.cd()

            # frame.GetXaxis().SetTitle(f" {le} : m({bsp}) [MeV]")
            frame.Draw()

        save_png(p, f"mc_fits", f"{spec}_{shape}", rpflag = 1)
