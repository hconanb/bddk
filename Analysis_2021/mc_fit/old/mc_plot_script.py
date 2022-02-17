import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021")
from essentials import *

def plot_mc(ws_name, shape, pullflag = 0):

    fws_base_plot_file = ROOT.TFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_{ws_name}_{shape}.root")
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
    if shape == "GAddBGEP":
        fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{ws_name}_{shape}_fit_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("dg_a") )
        fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{ws_name}_{shape}_fit_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("dg_b") )


    d_chi2 = frame.chiSquare("total_fit", "data")
    print(d_chi2)

    if "z_01" in ws_name:
        le = "D^{-} D^{+} K^{*0}"
        bsp = "B^{0}"
    elif "z_02" in ws_name or "z_03" in ws_name:
        le = "(D^{*-} #rightarrow D- #bar{#pi^{0}}) D^{+} + D^{-} (D^{*+} #rightarrow D+ #pi^{0})"
        bsp = "B^{0}"
    elif "z_04" in ws_name:
        le = "(D^{*-} #rightarrow D- #bar{#pi^{0}}) (D^{*+} #rightarrow D+ #pi^{0})"
        bsp = "B^{0}"
    elif "zz_04" in ws_name:
        le = "(D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) (D^{*+} #rightarrow D^{0} #pi^{+})"
        bsp = "B^{0}"
    elif "zz_07" in ws_name:
        le = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi^{+})"
        bsp = "B^{+}"
    elif "zz_08" in ws_name:
        le = "(#bar{D^{0}} #rightarrow #bar{D^{0}} #bar{#pi^{0}}) (D^{*+} #rightarrow D^{0} #pi^{+})"
        bsp = "B^{+}"
    elif "zz_09" in ws_name:
        le = "(#bar{D^{0}} #D^{0})"
        bsp = "B^{0}"
    elif "zz_10" in ws_name:
        le = "(#bar{D^{0}} #D^{0})"
        bsp = "B^{0}"
    elif "zz_12" in ws_name:
        le = "(#bar{D^{0}} #rightarrow #bar{D^{0}} #bar{#pi^{0}}) (D^{0} #rightarrow D^{0} #pi^{0})"
        bsp = "B^{0}"
    else:
        le = "temp"
        bsp = "B^{0}"

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
        pull.GetXaxis().SetTitle(f"{bsp} [MeV]")

        dtpave_dec = ROOT.TPaveText(0.05, 0.05, 0.7, 0.15, "NB NDC")
        dtpave_dec.SetFillStyle(0)
        dtpave_dec.AddText(f"MC: {le}")
        dtpave_dec.SetBorderSize(0)
        dtpave_dec.Draw()
        frame.addObject(dtpave_dec)


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
