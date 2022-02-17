import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *

p0_mc_list = [
    "01_z_11198006",
    "05_p_12197023",
    "09_zz_11196019",
    "13_s_13198040",
    "norm7_norm7_12197008",
    "norm8_norm8_11198007",
]

p1_mc_list = ["02_p_11198005", "06_p_12197410", "10_zz_11196413", "14_s_13198200","07_p_12197400","07_zz_12197024","07_st_12197024","15_s_13198400","02_z_11198400"]

p2_mc_list = [
    "04_z_11198401",
    "04_p_11198410",
    "04_zz_11198022",
    "04_st_11198022",
    "08_p_12197401",
    "08_zz_12197422",
    "08_st_12197422",
    "12_zz_11196414",
    "16_s_13198600",
]

def get_free_shapes(ws, spec, mean_start, window, shape_flag):

    b_low = mean_start - window
    b_high = mean_start + window

    #Gaussian
    if shape_flag == "G":
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_fit(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}], width_{spec}_{shape_flag}[5.0,0.0,20.0])"
        )

    #Double_Gaussian
    if shape_flag == "DG":
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_a(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}], width_{spec}_{shape_flag}[12, 0.0, 50.0])"
        )
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_b(B_DTF_M, mean_{spec}_{shape_flag}, width_b_{spec}_{shape_flag}[21.0, 0.0, 50.0])"
        )
        ws.factory(f"SUM::{spec}_{shape_flag}_fit({spec}_{shape_flag}_a_frac[0.8,0,1]*{spec}_{shape_flag}_a, {spec}_{shape_flag}_b)")

    # if shape_flag == "GEP":
    #     ws.factory(
    #         f"RooGaussExp::{spec}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[20,1.0,50.0],alpha_{spec}_{shape_flag}[1,0.00,10])"
    #     )
    # if shape_flag == "BGEP":
    #     ws.factory(
    #         f"RooBifurGaussExp::{spec}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_L_{spec}_{shape_flag}[20,0.00,50.0],width_R_{spec}_{shape_flag}[5,0.00,50.0],alpha_1_{spec}_{shape_flag}[1,0.00,10],alpha_2_{spec}_{shape_flag}[2,0.00,10])"
    #     )
    # if shape_flag == "BG":
    #     ws.factory(
    #         f"BifurGauss::{spec}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_1_{spec}_{shape_flag}[20,0.00,50.0],width_2_{spec}_{shape_flag}[5,0.00,50.0])"
    #     )
    # if shape_flag == "cb1R":
    #     ws.factory(
    #         f"CBShape::{spec}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,0,50],alpha_{spec}_{shape_flag}[-1,-5,-0.0], n_{spec}_{shape_flag}[1,0,100])"
    #     )
    # if shape_flag == "cb1L":
    #     ws.factory(
    #         f"CBShape::{spec}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,0,50], alpha_{spec}_{shape_flag}[1,0.0,5.0], n_{spec}_{shape_flag}[1,0,100])"
    #     )
    # if shape_flag == "cb2":
    #     ws.factory(
    #         f"CBShape::{spec}_a_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,0,100],alpha_l_{spec}_{shape_flag}[-1,-20,-0.0],  n1_{spec}_{shape_flag}[20,0,1000])"
    #     )
    #     ws.factory(
    #         f"CBShape::{spec}_b_fit(B_DTF_M,mean_{spec}_{shape_flag},                                                  width_{spec}_{shape_flag},         alpha_r_{spec}_{shape_flag}[1,0.0,20.0],   n2_{spec}_{shape_flag}[20,0,1000])"
    #     )
    #     ws.factory(
    #         f"SUM::{spec}_fit({spec}_a_frac[0.5,0,1]*{spec}_a_fit, {spec}_b_fit)"
    #     )
    #

def build_mc_ws(ws_name, event_list, shape_list, mean_window, b_start, b_min, b_max, wflag = 0):

    mcws = ROOT.RooWorkspace(ws_name)
    mcws.factory(f"B_DTF_M[{b_min}, {b_max}]")
    for shape_flag in shape_list:
        get_free_shapes(mcws, ws_name, b_start, mean_window, shape_flag)

    b_dtf_m = mcws.var("B_DTF_M")

    if wflag == 0:
        tc = ROOT.TChain(f"DecayTreeTuple")
        fl = grab_file_list("MC", event_list)
        for file_name in fl:
            tc.Add(file_name)
            print(file_name)
        data_args = ROOT.RooArgSet(b_dtf_m)
        data_set = ROOT.RooDataSet(f"{ws_name}_events", f"{ws_name}_events", tc, data_args)
        mcws.Import(data_set)
        mcws.writeToFile(f"{analysis_path}/mc_fit/base_mc_files/{ws_name}.root")
        print(f"Saved: {analysis_path}/mc_fit/base_mc_files/{ws_name}.root")
        mcws.Print()

    # if wflag > 0:
    #
    #     data_set_list = []
    #
    #     tc_1 = ROOT.TChain(f"DecayTreeTuple")
    #     tc_2 = ROOT.TChain(f"DecayTreeTuple")
    #
    #     fl_1 = grab_file_list("MC", [event_list[0]])
    #     for file_name in fl_1:
    #         tc_1.Add(file_name)
    #
    #     fl_2 = grab_file_list("MC", [event_list[1]])
    #     for file_name in fl_2:
    #         tc_2.Add(file_name)
    #
    #     if len(event_list) == 3:
    #         fl_3 = grab_file_list("MC", [event_list[2]])
    #         for file_name in fl_3:
    #             tc_3.Add(file_name)
    #
    #
    #     vallist = np.linspace(0, 1, wflag, endpoint = False)
    #
    #     name_valist = list(range(0, wflag))
    #
    #     if len(event_list) == 2:
    #         for wval_1, nid in zip(vallist, name_valist):
    #
    #             wval_2 = 1 - wval_1
    #
    #             data_args_1 = ROOT.RooArgSet(b_dtf_m)
    #             data_args_2 = ROOT.RooArgSet(b_dtf_m)
    #
    #             data_set_1 = ROOT.RooDataSet(f"{name_mc_ws}_events_{nid}", f"{name_mc_ws}_events_{nid}", tc_1, data_args_1)
    #             data_set_2 = ROOT.RooDataSet(f"{name_mc_ws}_events_2_{nid}_temp", f"{name_mc_ws}_events_2_{nid}_temp", tc_2, data_args_2)
    #
    #             wvar_1 = ROOT.RooRealVar("wvar_1","wvar_1", wval_1)
    #             wvar_2 = ROOT.RooRealVar("wvar_2","wvar_2", wval_2)
    #
    #             data_set_1.addColumn(wvar_1)
    #             data_set_2.addColumn(wvar_2)
    #
    #             data_args_1.add(wvar_1)
    #             data_args_2.add(wvar_2)
    #
    #             new_data_set_1 = ROOT.RooDataSet(data_set_1.GetName(), data_set_1.GetTitle(), data_set_1, data_args_1, "B_DTF_M > 0", "wvar_1");
    #             new_data_set_2 = ROOT.RooDataSet(data_set_2.GetName(), data_set_2.GetTitle(), data_set_2, data_args_2, "B_DTF_M > 0", "wvar_2");
    #
    #             new_data_set_1.append(new_data_set_2)
    #             data_set_list.append(new_data_set_1)



    else:
        for ds in data_set_list:
            mcws.Import(ds)
        mcws.writeToFile(f"../mc_ws_root_files/{name_mc_ws}.root")
        mcws.Print()

def run_mc_fit(ws_name, shape_list, wflag = 0):

        mcws_base_file = ROOT.TFile(f"{analysis_path}/mc_fit/base_mc_files/{ws_name}.root")
        mcws = mcws_base_file.Get(f"{ws_name}")
        b_dtf_m = mcws.var("B_DTF_M")
        data_set = mcws.data(f"{ws_name}_events")

        fitws = ROOT.RooWorkspace(f"fit_{ws_name}")
        fitws.Import(data_set)

        for shape in shape_list:
            fit = mcws.pdf(f"{ws_name}_{shape}_fit")
            if wflag == 0:
                fit_pdf = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
                fitws.Import(fit_pdf)
                fitws.Import(fit)
                fitws.Print()
                fitws.writeToFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_{ws_name}.root")

        #
        # if wflag !=0:
        #     name_valist = list(range(0,wflag))
        #     for nid in name_valist:
        #
        #         fitws = ROOT.RooWorkspace(f"fit_{name}_{nid}")
        #
        #         data_set = mcws.data(f"{name}_events_{nid}")
        #         fit_pdf = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.SumW2Error(True))
        #
        #
        #         fitws.Import(fit)
        #         fitws.Import(data_set)
        #
        #         fitws.Print()
        #         fitws.writeToFile(f"../mc_ws_root_files/fit_{name}_{nid}.root")



def plot_mc(ws_name, shape, pullflag = 0):

    fws_base_plot_file = ROOT.TFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_{ws_name}.root")
    fws = fws_base_plot_file.Get(f"fit_{ws_name}")

    b_dtf_m = fws.var("B_DTF_M")
    frame = b_dtf_m.frame(ROOT.RooFit.Title(ws_name), ROOT.RooFit.Bins(p_bbins))

    spectrum = frame.GetTitle()
    sspectrum = spectrum.split("_")[0]

    data_set = fws.data(f"{ws_name}_events")
    data_set.plotOn(frame, ROOT.RooFit.Name("data"))



    fit_pdf = fws.pdf(f"{ws_name}_{shape}_fit")
    fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{ws_name}_{shape}_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("total_fit") )

    if shape == "DG":
        fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{ws_name}_{shape}_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("dg_a") )
        fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{ws_name}_{shape}_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("dg_b") )

    d_chi2 = frame.chiSquare("total_fit", "data")
    print(d_chi2)

    dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.Draw()
    frame.addObject(dtpave)

    if pullflag == 1:

        p = residualPlot()
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
        if shape == "DG":
            legend.AddEntry("dg_a","g_a","l")
            legend.AddEntry("dg_b","g_b","l")

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
        pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
        pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset()/p.blabelratio)
        pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength()*p.blabelratio)
        pull.GetYaxis().SetTitle("Pull")
        pull.GetYaxis().SetNdivisions(202,False)
        pull.GetYaxis().SetRangeUser(-3,3)
        pull.GetYaxis().CenterTitle()
        pull.Draw("AP")
        pull.GetXaxis().SetTitle(f"m_{{B}} [MeV]")

    if pullflag == 0:

        p = ROOT.TCanvas("p1","p1")
        p.cd()

        title = "test"
        frame.GetXaxis().SetTitle(f" m({title}) [MeV]")
        frame.Draw()

    save_png(p, f"mc_fits", f"{ws_name}_{shape}", rpflag = 1)

def plot_mc_fit_test(name, shape, nf, wflag):

    name_valist = list(range(0,wflag))
    framelist = []

    p = ROOT.TCanvas("p1","p1")
    p.Divide(4,4)
    p.cd(0)
    wslist = []
    dslist = []
    framelist = []
    bvarlist = []
    filelist = []
    fit_pdflist = []

    chi2_list = []

    for nid in name_valist:

        filelist.append(ROOT.TFile(f"fit_ws_root_files/fit_{name}_{nid}.root"))
        wslist.append(filelist[nid].Get(f"fit_{name}_{nid}"))
        bvarlist.append(wslist[nid].var("B_DTF_M"))
        framelist.append(bvarlist[nid].frame(ROOT.RooFit.Title(f"{name}_{nid}"), ROOT.RooFit.Bins(p_bbins)))

        spectrum = framelist[nid].GetTitle()
        sspectrum = spectrum.split("_")[0]

        dslist.append(wslist[nid].data(f"{name}_events_{nid}"))
        dslist[nid].plotOn(framelist[nid], ROOT.RooFit.Name(f"data_{nid}"))

        p.cd(nid+1)

        fit_pdflist.append(wslist[nid].pdf(f"{name}_{nf}_fit"))
        fit_pdflist[nid].plotOn(framelist[nid], ROOT.RooFit.Components(f"{name}_{nf}_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name(f"total_fit_{nid}") )

        d_chi2 = framelist[nid].chiSquare(f"total_fit_{nid}", f"data_{nid}")
        chi2_list.append((nid, d_chi2))

        print(nid, d_chi2)

        dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
        dtpave.SetFillStyle(0)
        dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
        dtpave.Draw()
        framelist[nid].addObject(dtpave)
        framelist[nid].Draw()

    save_png(p, f"mc_fits_tests", f"{name}_{shape}", rpflag = 0)

# shape_list =
#  # ["G","DG","BG", "GEP", "BGEP", "cb1R", "cb1L", "cb2"]
