import sys
from createXFD import *

sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/")
from essentials import *


def plot_d_fit(spec, dvar, dstring, fancy_dstring, frame, dstrat):

    uid = f"{spec}_{dstring}"
    p = residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)

    leg = ROOT.TLegend(0.75,0.5,0.95,0.87)
    # leg.SetFillColor(kWhite)
    # leg.SetLineColor(kWhite)
    leg.AddEntry(f"{uid}_data","Data", "P");
    leg.AddEntry(frame.findObject(f"{uid}_fit"),"Signal + background","L");
    leg.AddEntry(frame.findObject(f"{uid}_bkg"),"Background only", "L");
    leg.AddEntry(frame.findObject(f"{uid}_signal"),"Signal only", "L");
    if dstrat == "e_dg":
        leg.AddEntry(frame.findObject(f"{uid}_siga"),"Signal A", "L");
        leg.AddEntry(frame.findObject(f"{uid}_sigb"),"Signal B", "L");
    frame.Draw()
    leg.Draw()
    p.pb.cd()

    hpull = frame.pullHist(f"{uid}_data", f"{uid}_fit")
    pull = dvar.frame(ROOT.RooFit.Range(f"{uid}_range"))
    pull.addPlotable(hpull, "P")
    pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
    pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize() * p.blabelratio)
    pull.GetYaxis().SetLabelSize(
        ROOT.gStyle.GetLabelSize("Y") * (1 + p.padratio) / (2 * p.padratio)
    )
    pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
    pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
    pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset() / p.blabelratio)
    pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength() * p.blabelratio)
    pull.GetYaxis().SetTitle("Pull")
    pull.GetYaxis().SetNdivisions(202, False)
    pull.GetYaxis().SetRangeUser(-3, 3)
    pull.GetYaxis().CenterTitle()
    pull.Draw("AP")

    pull.GetXaxis().SetTitle(f"{fancy_dstring} Mass [MeV]")
    save_png(p, "d_window_plots", f"{uid}_{dstrat}", rpflag=1)


def build_d_fit(spec, dvar, dstring, mydmass, data_set, dstrat):
    uid = f"{spec}_{dstring}"
    dws.factory(f"StringVar:{uid}_strat({dstrat})")
    if dstrat == "e_g":
        dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -2, 3])")
        dws.factory(
            f"Gaussian::{uid}_signal({dstring}_M, mean_{uid}[{mydmass},{mydmass-10},{mydmass+10}], width_{uid}[10.0,0.01,20])"
        )
        dws.factory(
            f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal,n_{uid}_bkg[100000,0,100000000]*{uid}_bkg)"
        )
        bkgfit = "Exponential"
        sigfit = "Gaussian"
    if dstrat == "e_bg":
        dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -2, 3])")
        dws.factory(
            f"BifurGauss::{uid}_signal({dstring}_M, mean_{uid}[{mydmass},{mydmass}-15,{mydmass}+15], width_{uid}_a[4.0,0.01,20], width_{uid}_b[10.0,0.01,20])"
        )
        dws.factory(
            f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal,n_{uid}_bkg[200000,0,100000000]*{uid}_bkg)"
        )
        bkgfit = "Exponential"
        sigfit = "Bifurcated Gaussian"
    if dstrat == "e_dg":
        dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -2, 3])")
        dws.factory(
            f"Gaussian::{uid}0_a({dstring}_M, mean_{uid}[{mydmass},{mydmass-10},{mydmass+10}], width_a_{uid}0[8.0,0.01,20])"
        )
        dws.factory(
            f"Gaussian::{uid}0_b({dstring}_M, mean_{uid}, width_b_{uid}0[7.0,0.01,20])"
        )
        dws.factory(f"SUM::{uid}_signal({uid}0_a_frac[0.5,0,1]*{uid}0_a, {uid}0_b)")
        dws.factory(
            f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal,n_{uid}_bkg[200000,0,100000000]*{uid}_bkg)"
        )
        bkgfit = "Exponential"
        sigfit = "Double Gaussian"
    if dstrat == "e_dg_st":
        dws.factory(
            f"Gaussian::{uid}0_a({dstring}_M, mean_{uid}[{mydmass},{mydmass-10},{mydmass+10}], width_a_{uid}0[8.0,0.01,20])"
        )
        dws.factory(
            f"Gaussian::{uid}0_b({dstring}_M, mean_{uid}, width_b_{uid}0[7.0,0.01,20])"
        )
        dws.factory(f"SUM::{uid}_signal({uid}0_a_frac[0.5,0,1]*{uid}0_a, {uid}0_b)")
        dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal)")
        bkgfit = "None"
        sigfit = "Double Gaussian"
    if dstrat == "e_cb_st":
        dws.factory(
            f"CBShape::{uid}_signal({dstring}_M, mean_{uid}[{mydmass},{mydmass-10},{mydmass+10}], width_1_{spec}[6.0,0.01,10], alpha_1_{uid}[1,0.1,10], n1_{uid}[1,0,10])"
        )
        dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal)")
        bkgfit = "None"
        sigfit = "Crystal Ball Right Sided"
    if dstrat == "e_bg_st":
        dws.factory(
            f"BifurGauss::{uid}_signal({dstring}_M, mean_{uid}[{mydmass},{mydmass}-15,{mydmass}+15], width_{uid}_a[4.0,0.01,20], width_{uid}_b[10.0,0.01,20])"
        )
        dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal)")
        bkgfit = "None"
        sigfit = "Bifurcated Gaussian"
    if dstrat == "e_ge_st":
        dws.factory(
            f"RooGaussExp::{uid}_signal({dstring}_M,mean_{uid}[{mydmass},{mydmass}-3,{mydmass}+3],width_{uid}[15,0.01,10],alpha2e_{uid}[10,0.01,10])"
        )
        dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal)")
        bkgfit = "None"
        sigfit = "Gaussian Exponential"
    if dstrat == "e_bge_st":
        dws.factory(
            f"RooBifurGaussExp::{uid}_signal({dstring}_M,mean_{uid}[{mydmass},{mydmass}-3,{mydmass}+3], width_{uid}[15,0.01,10], width_2_{uid}[10,0.01,100], alpha2e_{uid}[10,0.01,10], alpha2e2_{uid}[5,0.01,10])"
        )
        dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal)")
        bkgfit = "None"
        sigfit = "Gaussian Exponential"
    print(f"built {uid}")

    d_model = dws.pdf(f"{uid}_fit")
    d_model.fitTo(data_set, ROOT.RooFit.PrintLevel(0))
    d_frame = dvar.frame()

    data_set.plotOn(d_frame, ROOT.RooFit.Name(f"{uid}_data"))

    d_model.plotOn(
        d_frame,
        ROOT.RooFit.Components(f"{uid}_fit"),
        ROOT.RooFit.LineStyle(ROOT.kDashed),
        ROOT.RooFit.LineColor(ROOT.kGreen),
        ROOT.RooFit.Name(f"{uid}_fit"),
    )
    d_model.plotOn(
        d_frame,
        ROOT.RooFit.Components(f"{uid}_signal"),
        ROOT.RooFit.LineStyle(ROOT.kDashed),
        ROOT.RooFit.LineColor(ROOT.kBlue),
        ROOT.RooFit.Name(f"{uid}_signal"),
    )
    if dstrat == "e_dg":
        d_model.plotOn(
            d_frame,
            ROOT.RooFit.Components(f"{uid}0_a"),
            ROOT.RooFit.LineStyle(ROOT.kSolid),
            ROOT.RooFit.LineColor(ROOT.kOrange),
            ROOT.RooFit.Name(f"{uid}_siga"),
        )
        d_model.plotOn(
            d_frame,
            ROOT.RooFit.Components(f"{uid}0_b"),
            ROOT.RooFit.LineStyle(ROOT.kSolid),
            ROOT.RooFit.LineColor(ROOT.kViolet),
            ROOT.RooFit.Name(f"{uid}_sigb"),
        )
    d_model.plotOn(
        d_frame,
        ROOT.RooFit.Components(f"{uid}_bkg"),
        ROOT.RooFit.LineStyle(ROOT.kDashed),
        ROOT.RooFit.LineColor(ROOT.kRed),
        ROOT.RooFit.Name(f"{uid}_bkg"),
    )

    d_chi2 = d_frame.chiSquare(f"{uid}_fit", f"{uid}_data")
    d_mean = dws.obj(f"mean_{uid}")
    d_nyield = dws.var(f"n_{uid}_signal")
    d_sigfrac = dws.var(f"{uid}0_a_frac")
    d_stda = dws.var(f"width_a_{uid}0")
    d_stdb = dws.var(f"width_b_{uid}0")

    # print(dstrat)
    if "_st" not in dstrat:
        d_nbkg = dws.var(f"n_{uid}_bkg")

    dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.AddText(f"mean: {round(d_mean.getValV(),3)}")
    dtpave.AddText(f"nsig_a: {round(d_nyield.getValV()*d_sigfrac.getValV(),3)}")
    dtpave.AddText(f"nsig_b: {round(d_nyield.getValV()*(1-d_sigfrac.getValV()),3)}")
    dtpave.AddText(f"width_a: {round(d_stda.getValV(),3)}")
    dtpave.AddText(f"width_b: {round(d_stdb.getValV(),3)}")
    if round(d_nyield.getValV()*d_sigfrac.getValV(),3) > round(d_nyield.getValV()*(1-d_sigfrac.getValV()),3):
        dtpave.AddText(f"D Window opt 1: {round(d_stda.getValV(),3)*2}")
    else:
        dtpave.AddText(f"D Window opt 1: {round(d_stdb.getValV(),3)*2}")
    dtpave.AddText(f"D Window opt 2: {round(d_stda.getValV()*d_sigfrac.getValV() + d_stdb.getValV()*(1-d_sigfrac.getValV()),3)*2}")

    # dtpave.AddText(f"nsig: {round(d_nyield.getValV(),3)}")
    # if "_st" not in dstrat:
    #     print(dstrat)
    #     dtpave.AddText(f"nbkg: {round(d_nbkg.getValV(),3)}")
    # dtpave.AddText(f"Signal Fit: {sigfit}")
    # dtpave.AddText(f"Background Fit: {bkgfit}")


    d_frame.addObject(dtpave)
    return d_frame


speclist = ["st"]
#, "s", "norm7", "norm8"]
#["z", "zz", "p", "m",]


for spec in speclist:

    dname = f"d_{spec}_mass_fits"
    dws = ROOT.RooWorkspace(dname)
    if spec == "z":
        myd1mass = dpmass
        myd2mass = dpmass
        d1id = "D^{-}"
        d2id = "D^{+}"
        d1strat = "e_dg"
        d2strat = "e_dg"
    if spec == "zz":
        myd1mass = d0mass
        myd2mass = d0mass
        d1id = "#bar{D^{0}}"
        d2id = "D^{0}"
        d1strat = "e_dg"
        d2strat = "e_dg"
    if spec == "p":
        myd1mass = d0mass
        myd2mass = dpmass
        d1id = "#bar{D^{0}}"
        d2id = "D^{+}"
        d1strat = "e_dg"
        d2strat = "e_dg"
    if spec == "m":
        myd1mass = dpmass
        myd2mass = d0mass
        d1id = "D^{-}"
        d2id = "D^{0}"
        d1strat = "e_dg"
        d2strat = "e_dg"
    if spec == "st":
        myd1mass = d0mass
        myd2mass = d0mass
        myd2stmass = dstpmass
        myd2stmmass = 145.5
        d1id = "#bar{D^{0}}"
        d2id = "D^{0}"
        d2stid = "(D^{*+} #rightarrow D^{0} #pi+)"
        d2stmid = "D^{*+} - D^{0}"
        d1strat = "e_dg_st"
        d2strat = "e_dg_st"
        d2stmstrat = "e_dg_st"
    if spec == "s":
        myd1mass = dsmass
        myd2mass = dpmass
        d1id = "D^{-}_{s}"
        d2id = "D^{+}"
        d1strat = "e_dg"
        d2strat = "e_dg"
    if spec == "norm7":
        myd1mass = d0mass
        myd2mass = d0mass
        d1id = "#bar{D^{0}}"
        d2id = "D^{0} #rightarrow k#pi#pi#pi"
        d1strat = "e_dg"
        d2strat = "e_dg"
    if spec == "norm8":
        myd1mass = dpmass
        myd2mass = d0mass
        d1id = "D^{-}"
        d2id = "D^{0} #rightarrow K#pi#pi#pi"
        d1strat = "e_dg"
        d2strat = "e_dg"
    dfithalf = 40
    dws.factory(f"d1_strat[{d1strat}]")
    dws.factory(f"d2_strat[{d2strat}]")
    dws.factory(f"D1_M[{myd1mass - dfithalf},{myd1mass + dfithalf}]")
    dws.factory(f"D2_M[{myd2mass - dfithalf},{myd2mass + dfithalf}]")
    if spec == "st":
        dws.factory(f"D2st_M[{myd2stmass - dfithalf},{myd2stmass + dfithalf}]")
        dws.factory(f"DstmD_M[{142},{150}]")

    d1_mass = dws.var("D1_M")
    d2_mass = dws.var("D2_M")
    d2st_mass = dws.var("D2st_M")
    d2stm_mass = dws.var("DstmD_M")

    if spec != "st":
        data_args = ROOT.RooArgSet(d1_mass, d2_mass)
        data_file = ROOT.TFile(data_basepath + f"{spec}_spectrum.root")
    if spec == "st":
        CreateXFD(
            f"/mnt/c/Users/Harris/Desktop/rootfiles/data_run2/{spec}_spectrum_prexfd.root",
            "DecayTreeTuple",
            f"/mnt/c/Users/Harris/Desktop/rootfiles/data_run2/{spec}_spectrum.root",
            "DecayTreeTuple",
        )
        print(
            f"created new st tuple: /mnt/c/Users/Harris/Desktop/rootfiles/data_run2/{spec}_spectrum.root"
        )
        data_args = ROOT.RooArgSet(d1_mass, d2_mass, d2st_mass, d2stm_mass)
        data_file = ROOT.TFile(data_basepath + f"{spec}_spectrum.root")

    tree = data_file.Get("DecayTreeTuple")

    data_set = ROOT.RooDataSet(f"{spec}_data", f"{spec}_data", tree, data_args)

    d1f = build_d_fit(spec, d1_mass, "D1", myd1mass, data_set, d1strat)
    plot_d_fit(spec, d1_mass, "D1", d1id, d1f, d1strat)
    d2f = build_d_fit(spec, d2_mass, "D2", myd2mass, data_set, d2strat)
    plot_d_fit(spec, d2_mass, "D2", d2id, d2f, d2strat)
    if spec == "st":
        d2stmf = build_d_fit(
            spec, d2stm_mass, "DstmD", myd2stmmass, data_set, d2stmstrat
        )
        plot_d_fit(spec, d2stm_mass, "DstmD", d2stmid, d2stmf, d2stmstrat)

    dws.writeToFile(f"d_window_root_files/{dname}.root")
    dws.Print()
