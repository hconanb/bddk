import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
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
def getnset(dws,spec):
    list = dws.allFunctions()
    pset = ROOT.RooArgSet()
    for i in list:
        fname = i.GetName()
        specialcase = spec+"z"
        if spec in fname and specialcase not in fname and "base" not in fname and "norm" not in fname:
            pset.add(i)
    bkg = dws.obj(f"n_{spec}_bkg")
    pset.add(bkg)
    return pset

fom_z_list = []
fom_zz_list = []
fom_p_list = []
fom_m_list = []
fom_st_list = []
fom_v1_list = []
fom_v2_list = []
rows = []
tn=[]


for i in range(0,15):

    rows.append(i)
    dws_base_plot_file = ROOT.TFile(f"{fit_ws_basepath}/{fit_ws_name}_{i}.root")
    dws = dws_base_plot_file.Get(f"{fit_ws_name}")

    # b_dtf_m = dws.var("B_DTF_M")
    file1 = open(f"txts/MyFile_{i}.txt","w")

    file1.write(str(fom_calc(dws, "z0", "n_1_z", "z")) + "\n")
    file1.write(str(fom_calc(dws, "z1", "n_23_z", "z")) + "\n")
    file1.write(str(fom_calc(dws, "z2", "n_4_z", "z")) + "\n")
    file1.write(str(fom_calc(dws, "zz0", "n_9_zz", "zz")) + "\n")
    file1.write(str(fom_calc(dws, "zz1", "n_71011_zz", "zz")) + "\n")
    file1.write(str(fom_calc(dws, "zz2", "n_4812_zz", "zz")) + "\n")

    file1.write(str(fom_calc(dws, "p0", "n_5_p", "p")) + "\n")
    file1.write(str(fom_calc(dws, "p1", "n_267_p", "p")) + "\n")
    file1.write(str(fom_calc(dws, "p2", "n_48_p", "p")) + "\n")

    file1.write(str(fom_calc(dws, "m1", "n_3_m", "m")) + "\n")
    file1.write(str(fom_calc(dws, "m2", "n_4_m", "m")) + "\n")

    file1.write(str(fom_calc(dws, "st0", "n_7_st", "st")) + "\n")
    file1.write(str(fom_calc(dws, "st1", "n_48_st", "st")) + "\n")

    file1.write("")
    file1.close()

    b_dtf_m = dws.var("B_DTF_M")
    b_dtf_m.setRange("myrange", bmin, bmax)

    all_cats = dws.cat("all_cats")
    all_data_sets = dws.data("all_data_sets")
    all_fit = dws.pdf("super_fit_Pdf")

    z_frame = b_dtf_m.frame(ROOT.RooFit.Title("z_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    zz_frame = b_dtf_m.frame(ROOT.RooFit.Title("zz_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    p_frame = b_dtf_m.frame(ROOT.RooFit.Title("p_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    m_frame = b_dtf_m.frame(ROOT.RooFit.Title("m_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
    st_frame = b_dtf_m.frame(ROOT.RooFit.Title("st_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(bbins))
# s_frame = b_dtf_m.frame(ROOT.RooFit.Title("s_spectrum"), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Bins(nbins))

    all_cats_Set = ROOT.RooArgSet(all_cats)

    z_pset = getnset(dws,"z")
    zz_pset = getnset(dws,"zz")
    p_pset = getnset(dws,"p")
    m_pset = getnset(dws,"m")
    st_pset = getnset(dws,"st")
    # s_pset = getnset(dws,"s")

    bset = ROOT.RooArgSet(z_pset, zz_pset)
    b2set = ROOT.RooArgSet(bset, p_pset)
    b3set = ROOT.RooArgSet(b2set, m_pset)
    ballset = ROOT.RooArgSet(b3set, st_pset)

    # z_fom, z_tn = fom_calc(z_pset)
    # zz_fom, zz_tn = fom_calc(zz_pset)
    # p_fom, p_tn = fom_calc(p_pset)
    # m_fom, m_tn = fom_calc(m_pset)
    # st_fom, st_tn = fom_calc(st_pset)
    # fom_v1 = (z_fom + zz_fom + p_fom + m_fom + st_fom)/5
    # fom_v2, v2_tn = fom_calc(ballset)
    # tn.append(z_tn)
    #
    # fom_z_list.append(z_fom)
    # fom_zz_list.append(zz_fom)
    # fom_p_list.append(p_fom)
    # fom_m_list.append(m_fom)
    # fom_st_list.append(st_fom)
    # fom_v1_list.append(fom_v1)
    # fom_v2_list.append(fom_v2)

    z_can = plot_all_fit(b_dtf_m, z_frame, i,  all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag)
    zz_can = plot_all_fit(b_dtf_m, zz_frame, i, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag)
    p_can = plot_all_fit(b_dtf_m, p_frame,  i,  all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag)
    m_can = plot_all_fit(b_dtf_m, m_frame,  i,  all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag)
    st_can = plot_all_fit(b_dtf_m, st_frame, i, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag)

    bigcan = ROOT.TCanvas("cc1","cc1")
    bigcan.Divide(2,3)
    bigcan.cd(1)
    z_can.Draw()
    bigcan.cd(2)
    zz_can.Draw()
    bigcan.cd(3)
    p_can.Draw()
    bigcan.cd(4)
    m_can.Draw()
    bigcan.cd(5)
    st_can.Draw()
    bigcan.cd(6)

    # fompt = ROOT.TPaveText(.05,.1,.95,.8)
    # fompt.SetFillStyle(0)
    # fompt.AddText(f"fom_v1: {fom_v1}")
    # fompt.AddText(f"fom_v2: {fom_v2}")
    # fompt.Draw()
    #
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}/'):
        os.makedirs(f'plots/{now.month}_{now.day}/')
    bigcan.SaveAs(f"plots/{now.month}_{now.day}/all_{i}.pdf")


# import pandas as pd
#
# #, fom_p_list, fom_m_list, fom_st_list, fom_v1_list, fom_v2_list
# #  'p', 'm', 'st', 'v1', 'v2'
#
# colums = [[fom_z_list, fom_zz_list]]
#
# df = pd.DataFrame(list(zip(*[fom_z_list, fom_zz_list, fom_p_list, fom_m_list, fom_st_list, fom_v1_list, fom_v2_list]))).add_prefix('Col')
#
# df.to_excel("out.xlsx")
#
# print(tn)
#
# plot_all_fit(b_dtf_m, s_frame, s_pset, all_data_sets, all_fit, all_cats, all_cats_Set, fbkg, f0, f1, f2, nflag)
