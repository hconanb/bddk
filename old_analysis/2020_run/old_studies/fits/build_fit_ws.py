import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

def build_data_fits(dws, bkg_flag, peak_0_flag, peak_1_flag, peak_2_flag, nflag):

    # dws.factory(f"Exponential:s_spectrum_bkg(B_DTF_M, c0_s[0, -2, 2])")
    dws.factory(f"Bernstein:s_spectrum_bkg(B_DTF_M,{{c0_s[1,0,10], c1_s[1,0,10], c2_s[1,0,10], c3_s[1,0,10]}})")
    dws.factory(f"Gaussian::s0_a(B_DTF_M, mean_s0[5366,5350,5375], width_a_s0[8.0,0.01,20])")
    dws.factory(f"Gaussian::s0_b(B_DTF_M, mean_s0, width_b_s0[7.0,0.01,20])")
    dws.factory(f"SUM::s0_fit(s0_a_frac[0.5,0,1]*s0_a, s0_b)")
    # dws.factory(f"BifurGauss::s1_fit(B_DTF_M,mean_s1[5230,5200,5250],width_1_sl[15,0.01,100],width_1_sr[10,0.01,100])")
    # dws.factory(f"BifurGauss::s2_fit(B_DTF_M,mean_s2[5090,5050,5130],width_2_sl[4,0.01,100],width_2_sr[5,0.01,100])")
    dws.factory(f"CBShape::s1_fit(B_DTF_M,mean_s1[5230,5200,5250],width_1_s[10,0,50],alpha_1_s[1,0.1,100], n1_s[3,1,100])")
    dws.factory(f"CBShape::s2_fit(B_DTF_M,mean_s2[5075,5050,5100],width_2_s[10,0,50],alpha_5_s[1,0.1,100], n5_s[3,1,100])")

    speclist = [
    "z",
    "zz",
    "p",
    "m",
    "st",
    ]

    for spec in speclist:

        if bkg_flag == "Exponential":
            dws.factory(f"Exponential:{spec}_spectrum_bkg(B_DTF_M, c0_{spec}[0, -2, 2])")
        if bkg_flag == "Chebychev":
            dws.factory(f"Chebychev:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[0.,-3,3],c1_{spec}[0.,-3,3]}})")
        if bkg_flag == "Bernstein":
            dws.factory(f"Bernstein:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[1,0,10], c1_{spec}[1,0,10], c2_{spec}[1,0,10], c3_{spec}[1,0,10]}})")

        if peak_0_flag == "DG":
            dws.factory(f"Gaussian::{spec}0_a(B_DTF_M, mean_{spec}0[5279,5270,5290], width_a_{spec}0[8.0,0.01,20])")
            dws.factory(f"Gaussian::{spec}0_b(B_DTF_M, mean_{spec}0, width_b_{spec}0[7.0,0.01,20])")
            dws.factory(f"SUM::{spec}0_fit({spec}0_a_frac[0.5,0,1]*{spec}0_a, {spec}0_b)")
        if peak_0_flag == "G":
            dws.factory(f"Gaussian::{spec}0_fit(B_DTF_M, mean_{spec}0[9,-20,20], width_a_{spec}0[400,200,1000])")

        if peak_1_flag == "GEP":
            dws.factory(f"RooGaussExp::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}l[15,0.01,100],alphae_{spec}[10,0.01,100])")
        if peak_1_flag == "BGEP":
            dws.factory(f"RooBifurGaussExp::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}l[15,0.01,100],width_1_{spec}R[15,0.01,100],alphae_{spec}[10,0.01,100],alphae2_{spec}[10,0.01,100])")
        if peak_1_flag == "BG":
            dws.factory(f"BifurGauss::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}l[15,0.01,100],width_1_{spec}r[10,0.01,100])")
        if peak_1_flag == "cb2":
            dws.factory(f"CBShape::{spec}1_a_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}[10,0,50],alpha_1_{spec}[-5,-100,-0.1], n1_{spec}[3,1,100])")
            dws.factory(f"CBShape::{spec}1_b_fit(B_DTF_M,mean_{spec}1,                width_1_{spec},         alpha_2_{spec}[11,0.01,100], n2_{spec}[6,1,100])")
            dws.factory(f"SUM::{spec}1_fit({spec}1_a_frac[0.5,0,1]*{spec}1_a_fit, {spec}1_b_fit)")
        if peak_1_flag == "cb1R":
            dws.factory(f"CBShape::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}[10,0,50],alpha_1_{spec}[-5,-100,-0.1], n1_{spec}[3,1,100])")
        if peak_1_flag == "cb1L":
            dws.factory(f"CBShape::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}[10,0,50],alpha_1_{spec}[1,0.1,100], n1_{spec}[3,1,100])")

        if peak_2_flag == "GEP":
            dws.factory(f"RooGaussExp::{spec}2_fit(B_DTF_M,mean_{spec}2[4999,4950,5020],width_2_{spec}l[15,0.01,100],alpha2e_{spec}[10,0.01,100])")
        if peak_2_flag == "BGEP":
            dws.factory(f"RooBifurGaussExp::{spec}2_fit(B_DTF_M,mean_{spec}2[4999,4950,5020],width_2_{spec}l[15,0.01,100],width_2_{spec}R[15,0.01,100],alpha2e_{spec}[10,0.01,100],alpha2e2_{spec}[10,0.01,100])")
        if peak_2_flag == "BG":
            dws.factory(f"BifurGauss::{spec}2_fit(B_DTF_M,mean_{spec}2[4999,4950,5020],width_2_{spec}l[4,0.01,100],width_2_{spec}r[5,0.01,100])")
        if peak_2_flag == "cb2":
            dws.factory(f"CBShape::{spec}2_a_fit(B_DTF_M,mean_{spec}2[4995,4950,5020],width_2_{spec}[10,0,50],alpha_3_{spec}[-5,-100,-0.1], n3_{spec}[3,1,100])")
            dws.factory(f"CBShape::{spec}2_b_fit(B_DTF_M,mean_{spec}2,                width_2_{spec},         alpha_4_{spec}[11,0.01,100], n4_{spec}[6,1,100])")
            dws.factory(f"SUM::{spec}2_fit({spec}2_a_frac[0.5,0,1]*{spec}2_a_fit, {spec}2_b_fit)")
        if peak_2_flag == "cb1R":
            dws.factory(f"CBShape::{spec}2_fit(B_DTF_M,mean_{spec}2[4995,4950,5020],width_2_{spec}[10,0,50],alpha_5_{spec}[-5,-100,-0.1], n5_{spec}[3,1,100])")
        if peak_2_flag == "cb1L":
            dws.factory(f"CBShape::{spec}2_fit(B_DTF_M,mean_{spec}2[4995,4950,5020],width_2_{spec}[10,0,50],alpha_5_{spec}[1,0.1,100], n5_{spec}[3,1,100])")

    if nflag == "calc":
        dws.factory("SUM::z_spectrum_all_fit(n_1_z*z0_fit,n_23_z*z1_fit,n_4_z*z2_fit,n_z_bkg[100,0,10000]*z_spectrum_bkg)")
        dws.factory("SUM::zz_spectrum_all_fit(n_9_zz*zz0_fit,n_71011_zz*zz1_fit,n_4812_zz*zz2_fit,n_zz_bkg[100,0,100000]*zz_spectrum_bkg)")
        dws.factory("SUM::p_spectrum_all_fit(n_5_p*p0_fit,n_267_p*p1_fit,n_48_p*p2_fit,n_p_bkg[100,0,10000]*p_spectrum_bkg)")
        dws.factory("SUM::m_spectrum_all_fit(n_3_m*m1_fit,n_4_m*m2_fit,n_m_bkg[100,0,10000]*m_spectrum_bkg)")
        dws.factory("SUM::st_spectrum_all_fit(n_7_st*st0_fit,n_48_st*st1_fit,n_st_bkg[100,0,1000]*st_spectrum_bkg)")
        dws.factory("SUM::s_spectrum_all_fit(n_13_s*s0_fit,n_1415_s*s1_fit,n_16_s*s2_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")
    if nflag == "test":
        dws.factory("SUM::z_spectrum_all_fit(n_1_t[100,0,10000]*z0_fit,n_23t[100,0,10000]*z1_fit,n_4at[100,0,10000]*z2_fit,n_z_bkg[100,0,10000]*z_spectrum_bkg)")
        dws.factory("SUM::zz_spectrum_all_fit(n_9t[100,0,10000]*zz0_fit,n_71011t[100,0,10000]*zz1_fit,n_4812t[100,0,10000]*zz2_fit,n_zz_bkg[100,0,100000]*zz_spectrum_bkg)")
        dws.factory("SUM::p_spectrum_all_fit(n_5t[100,0,10000]*p0_fit,n_267t[100,0,10000]*p1_fit,n_48t[100,0,10000]*p2_fit,n_p_bkg[100,0,10000]*p_spectrum_bkg)")
        dws.factory("SUM::m_spectrum_all_fit(n_3t[100,0,10000]*m1_fit,n_4bt[100,0,10000]*m2_fit,n_m_bkg[100,0,10000]*m_spectrum_bkg)")
        dws.factory("SUM::st_spectrum_all_fit(n_7t[100,0,10000]*st0_fit,n_48t[100,0,10000]*st1_fit,n_st_bkg[100,0,1000]*st_spectrum_bkg)")
        dws.factory("SUM::s_spectrum_all_fit(n_13t[100,0,10000]*s0_fit,n_1415t[100,0,10000]*s1_fit,n_16t[100,0,10000]*s2_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")

    dws.Print()
    return dws
def plot_all_fit(var, frame, wsi, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag, pullflag = 0):

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
        legend.SetHeader(title,"C")
        legend.AddEntry("pdf","Total Fit","lp")
        legend.AddEntry("fit0","Fully Reconstructed Peak","l")
        legend.AddEntry("fit1","1 missing particle peak","l")
        legend.AddEntry("fit2","2 missing particle peak","l")
        legend.AddEntry("bkg","Background","l")

        legend.Draw()

        p.pb.cd()

        hpull = frame.pullHist("data","pdf")
        pull = var.frame();
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
        now = datetime.datetime.now()
        if not os.path.exists(f'plots/{now.month}_{now.day}/'):
            os.makedirs(f'plots/{now.month}_{now.day}/')
        p.save(f"plots/{now.month}_{now.day}/{sspectrum}_{fbkg}_{f0}_{f1}_{f2}_{nflag}.pdf")
        p.save(f"plots/{now.month}_{now.day}/{sspectrum}_{fbkg}_{f0}_{f1}_{f2}_{nflag}.png")
    if pullflag == 0:

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

        now = datetime.datetime.now()
        if not os.path.exists(f'plots/{now.month}_{now.day}/fh/'):
            os.makedirs(f'plots/{now.month}_{now.day}/fh/')
        p.SaveAs(f"plots/{now.month}_{now.day}/fh/{spectrum}_{fbkg}_{f0}_{f1}_{f2}_{nflag}_{wsi}.pdf")
        # p.SaveAs(f"plots/{now.month}_{now.day}/fh/{spectrum}_{fbkg}_{f0}_{f1}_{f2}_{nflag}_{wsi}.png")
        return frame
def fom_calc(dws, peak_id, yield_id, type):

    b_dtf_m = dws.var("B_DTF_M")

    sig_pdf = dws.pdf(f"{peak_id}_fit")
    bkg_pdf = dws.pdf(f"{type}_spectrum_bkg")

    n_sig_var = dws.obj(f"{yield_id}")
    n_bkg_var = dws.obj(f"n_{type}_bkg")
    stype = dws.var(f"mean_{peak_id}")
    n_sig = n_sig_var.getVal()
    n_bkg = n_bkg_var.getVal()
    range = f"range_{peak_id}"

    b_dtf_m.setRange(range, stype.getVal()-30, stype.getVal()+30)

    ns = ROOT.RooFit.NormSet(b_dtf_m)
    sig_int = sig_pdf.createIntegral(ROOT.RooArgSet(b_dtf_m), ns, ROOT.RooFit.Range(range))

    ns = ROOT.RooFit.NormSet(b_dtf_m)
    bkg_int = bkg_pdf.createIntegral(ROOT.RooArgSet(b_dtf_m), ns, ROOT.RooFit.Range(range))
    #
    signal = sig_int.getVal()*n_sig
    background = bkg_int.getVal()*n_bkg
    fom = signal/math.sqrt(signal+background)

    print(signal)
    print(background)
    print(fom)
    return(fom)

def fit_ws():

    dws_base_file = ROOT.TFile(f"{data_ws_basepath}/{data_ws_name}.root")
    dws = dws_base_file.Get(data_ws_name)

    if gc_onflag == 1:
        clist = dws.allPdfs()

    all_cats = dws.cat("all_cats")
    b_dtf_m = dws.var("B_DTF_M")
    b_dtf_m.setRange("myrange", bmin, bmax)

    tds = dws.data("all_data_sets")
    fd_cut = 2/10*(i - 1)
    test_fd_cut = f"D1_FDCHI2_ORIVX > {fd_cut} && D2_FDCHI2_ORIVX > {fd_cut}"

    if i != 0:
        all_data_sets = tds.reduce(ROOT.RooArgSet(all_cats,b_dtf_m),f"{dira_cut} && {test_fd_cut}")
    else:
        all_data_sets = tds

    dws = build_data_fits(dws, fbkg, f0, f1, f2, nflag)

    z_model = dws.pdf("z_spectrum_all_fit")
    zz_model = dws.pdf("zz_spectrum_all_fit")
    p_model = dws.pdf("p_spectrum_all_fit")
    m_model = dws.pdf("m_spectrum_all_fit")
    st_model = dws.pdf("st_spectrum_all_fit")
    # s_model = dws.pdf("s_spectrum_all_fit")

    all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)
    all_fit.addPdf(z_model, "z_spectrum")
    all_fit.addPdf(zz_model, "zz_spectrum")
    all_fit.addPdf(p_model, "p_spectrum")
    all_fit.addPdf(m_model, "m_spectrum")
    all_fit.addPdf(st_model, "st_spectrum")
    # all_fit.addPdf(s_model, "s_spectrum")

    if gc_onflag == 0:
        all_fit.fitTo(all_data_sets, ROOT.RooFit.Range("myrange"), ROOT.RooFit.PrintLevel(0))
    if gc_onflag == 1:
        all_fit.fitTo(all_data_sets, ROOT.RooFit.Range("myrange"), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(clist))

    fit_dws = ROOT.RooWorkspace(fit_ws_name)

    getattr(fit_dws,'import')(all_data_sets)
    getattr(fit_dws,'import')(all_fit)

    fit_dws.writeToFile(f"{fit_ws_basepath}/{fit_ws_name}_{i}.root")

for i in range(0,15):
    fit_ws()
