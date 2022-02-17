import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

def plot_all_fit(var, frame, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, ttflag):

    spectrum = frame.GetTitle()
    sspectrum = spectrum.split("_")[0]

    all_data_sets.plotOn(frame, ROOT.RooFit.CutRange("myrange"), ROOT.RooFit.Cut(f"all_cats==all_cats::{spectrum}"), ROOT.RooFit.Name("data"))
    nData = all_data_sets.sumEntries("", "myrange")

    all_fit.plotOn(frame, ROOT.RooFit.Range("myrange"), ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange("myrange"), ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spectrum}_bkg"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg") )
    all_fit.plotOn(frame, ROOT.RooFit.Range("myrange"), ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange("myrange"), ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spectrum}_all_fit"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
    all_fit.plotOn(frame, ROOT.RooFit.Range("myrange"), ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange("myrange"), ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{sspectrum}0_fit"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
    all_fit.plotOn(frame, ROOT.RooFit.Range("myrange"), ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange("myrange"), ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{sspectrum}1_fit"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
    all_fit.plotOn(frame, ROOT.RooFit.Range("myrange"), ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.ProjectionRange("myrange"), ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{sspectrum}2_fit"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))


    d3 = "K^{*0}"
    if sspectrum == "z":
        title = "D^{-} D^{+} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{+}"
    if sspectrum == "zz":
        title = "#bar{D^{0}} D^{0} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{0}"
    if sspectrum == "p":
        title = "#bar{D^{0}} D^{+} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{+}"
    if sspectrum == "m":
        title = "D^{-} D^{0} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{0}"
    if sspectrum == "st":
        title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
        d1 = "D^{-}"
        d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
    if sspectrum == "s":
        title = "D^{-}_{s} D^{+} K^{*0}"
        d1 = "D^{-}_{s}"
        d2 = "D^{+}"
    # if pullflag == 1:
    #     p = residualPlot()
    #     p.pt.cd()
    #     xaxis = frame.GetXaxis()
    #     xaxis.SetTickLength(0)
    #     xaxis.SetNdivisions(0)
    #     xaxis.SetLabelSize(0)
    #     # frame.addObject(txt)
    #     frame.Draw()
    #
    #     legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
    #     legend.SetHeader(title,"C")
    #     legend.AddEntry("pdf","Total Fit","lp")
    #     legend.AddEntry("fit0","Fully Reconstructed Peak","l")
    #     legend.AddEntry("fit1","1 missing particle peak","l")
    #     legend.AddEntry("fit2","2 missing particle peak","l")
    #     legend.AddEntry("bkg","Background","l")
    #
    #     legend.Draw()
    #
    #     p.pb.cd()
    #
    #     hpull = frame.pullHist("data","pdf")
    #     pull = var.frame();
    #     pull.addPlotable(hpull, "P");
    #     pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
    #     pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()*p.blabelratio)
    #     pull.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Y")*(1+p.padratio)/(2*p.padratio))
    #     pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
    #     pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
    #     pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset()/p.blabelratio)
    #     pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength()*p.blabelratio)
    #     pull.GetYaxis().SetTitle("Pull")
    #     pull.GetYaxis().SetNdivisions(202,False)
    #     pull.GetYaxis().SetRangeUser(-3,3)
    #     pull.GetYaxis().CenterTitle()
    #     pull.Draw("AP")
    #
    #     pull.GetXaxis().SetTitle(f"m_{{B}} [MeV]")
    #     now = datetime.datetime.now()
    #     save_pdf(p, "full_fit", f"{spectrum}_{fbkg}_{f0}_{f1}_{f2}_{nflag}_{ttflag}")
    # if pullflag == 0:

    p = ROOT.TCanvas("p1","p1")
    p.cd()

    frame.GetXaxis().SetTitle(f" m({title}) [MeV]")
    frame.Draw()

    legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
    #legend.SetHeader(title,"C")
    legend.AddEntry("pdf","Total Fit","lp")
    legend.AddEntry("fit0","Fully Reconstructed Peak","l")
    legend.AddEntry("fit1","1 missing particle peak","l")
    legend.AddEntry("fit2","2 missing particle peak","l")
    legend.AddEntry("bkg","Background","l")
    # legend.AddEntry("fom",f"fom: {fom}")
    legend.SetTextSize(0.030)
    legend.Draw()
    frame.addObject(legend)

    save_pdf(p, "full_fit", f"{spectrum}_{fbkg}_{f0}_{f1}_{f2}_{nflag}_{ttflag}")
    return frame
def plot_full_histogram(ttflag):

    fws_base_plot_file = ROOT.TFile(f"../workspace_builds/fits/fit_{ttflag}.root")
    fws = fws_base_plot_file.Get(f"fit_{ttflag}")

    b_dtf_m = fws.var("B_DTF_M")
    b_dtf_m.setRange("myrange", bmin, bmax)

    all_cats = fws.cat("all_cats")
    all_data_sets = fws.data("all_data_sets")
    all_fit = fws.pdf("super_fit_Pdf")

    z_frame = b_dtf_m.frame(ROOT.RooFit.Title("z_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    zz_frame = b_dtf_m.frame(ROOT.RooFit.Title("zz_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    p_frame = b_dtf_m.frame(ROOT.RooFit.Title("p_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    m_frame = b_dtf_m.frame(ROOT.RooFit.Title("m_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    st_frame = b_dtf_m.frame(ROOT.RooFit.Title("st_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    # s_frame = b_dtf_m.frame(ROOT.RooFit.Title("s_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(nbins))

    all_cats_Set = ROOT.RooArgSet(all_cats)

    nflag = 0
    z_can = plot_all_fit(b_dtf_m, z_frame, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, ttflag)
    zz_can = plot_all_fit(b_dtf_m, zz_frame, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, ttflag)
    p_can = plot_all_fit(b_dtf_m, p_frame, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, ttflag)
    m_can = plot_all_fit(b_dtf_m, m_frame, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, ttflag)
    st_can = plot_all_fit(b_dtf_m, st_frame, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, ttflag)

plot_full_histogram("ToT")
plot_full_histogram("T")
plot_full_histogram("nTaT")
