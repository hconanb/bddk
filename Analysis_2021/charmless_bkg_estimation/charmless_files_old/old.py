import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
import os
import ROOT
import pandas as pd
import uproot
import glob
import datetime
import uproot

class spectrum_class:

  def __init__(self, name, d1_string, d2_string, d1_mass, d2_mass, b_string):
    self.spec = name

    self.d1_string = d1_string
    self.d2_string = d2_string

    self.d1_mass = d1_mass
    self.d2_mass = d2_mass

    if self.spec == "st":
        self.d3_string = "(D^{*+} #rightarrow D^{0} #pi+)"
        self.d3_mass = 150

    self.b_string = b_string

z_c = spectrum_class("Z_m_p", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
# zz_c = spectrum_class("Z_z_z", "#barD^{0}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D^{0} K*0")
# p_c = spectrum_class("P_z_p", "#bar{D^{0}}", "D^{+}", d0mass, dpmass, "B #rightarrow #bar{D^{0}} D+ K*0")
# m_c = spectrum_class("M_m_z", "D^{-}", "D^{0}", dpmass, d0mass,  "B #rightarrow D- D^{0} K*0")
z_c_01 = spectrum_class("01_Z_m_p_11198006", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
# st_c = spectrum_class("P_z_pst","#bar{D^{0}}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D*+ K*0")
# # s_c = spectrum_class("s","D^{-}_{s}", "D^{+}" ,dsmass, dpmass, "e_g")
# # n7_c = spectrum_class("norm7","#bar{D^{0}}", "D^{0} #rightarrow k#pi#pi#pi", d0mass, d0mass, "e_dg_b")
# n8_c = spectrum_class("norm8","D^{-}", "D^{0} #rightarrow k#pi#pi#pi", dpmass, d0mass, "B #rightarrow D-D+K*0")


def fit_est(sc, strat = "2peak"):

    ws_base_file = ROOT.TFile(f"charmless_files/{sc.spec}_e.root")
    tree = ws_base_file.Get(f"DecayTreeTuple")

    ws = ROOT.RooWorkspace(f"fit_{sc.spec}")
    if strat == "3peak":
        ws.factory("B_M[4800, 5600]")
    if strat == "2peak":
        ws.factory("B_M[5000, 5600]")
    # ws.factory(f"Exponential:{sc.spec}_spectrum_bkg(B_M, c0_{sc.spec}[0, -2, 2])")
    ws.factory(f"Chebychev:{sc.spec}_spectrum_bkg(B_M,{{c0_{sc.spec}[0.,-3,3],c1_{sc.spec}[0.,-3,3], c2_{sc.spec}[0.,-3,3], c3_{sc.spec}[0.,-3,3]}})")
    ws.factory(f"Gaussian::{sc.spec}_fit_1(B_M, mean_1[5280,5255,5295], width_1[5.0,0.1,50.0])")
    ws.factory(f"Gaussian::{sc.spec}_fit_2(B_M, mean_2[5100,5075,5125], width_2[5.0,0.1,50.0])")
    ws.factory(f"Gaussian::{sc.spec}_fit_3(B_M, mean_3[4900,4850,4950], width_3[5.0,0.1,50.0])")

    if strat == "3peak":
        ws.factory(f"SUM::{sc.spec}_all_fit(n_1[100,0,10000]*{sc.spec}_fit_1,n_2[100,0,10000]*{sc.spec}_fit_2,n_3[100,0,10000]*{sc.spec}_fit_3,n_bkg[100,0,10000]*{sc.spec}_spectrum_bkg)")
    if strat == "2peak":
        ws.factory(f"SUM::{sc.spec}_all_fit(n_1[100,0,10000]*{sc.spec}_fit_1,n_2[100,0,10000]*{sc.spec}_fit_2,n_bkg[100,0,10000]*{sc.spec}_spectrum_bkg)")

    b_m = ws.var("B_M")
    data_args = ROOT.RooArgSet(b_m)
    fit = ws.pdf(f"{sc.spec}_all_fit")
    data_set = ROOT.RooDataSet(f"{sc.spec}_data", f"{sc.spec}_data", tree, data_args)
    all_fit = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    ws.Import(all_fit)
    fitresult = ws.obj(f"fitresult_{sc.spec}_all_fit_{sc.spec}_data")

    frame = b_m.frame(ROOT.RooFit.Title(f"{sc.spec}"))
    data_set.plotOn(frame, ROOT.RooFit.Name("data"))
    fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))
    fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_all_fit"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
    fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_fit_1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("f1"))
    fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_fit_2"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet), ROOT.RooFit.Name("f2"))
    if strat == "3peak":
        fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_fit_3"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kMagenta), ROOT.RooFit.Name("f3"))

    p = ROOT.TCanvas("p1","p1")
    p.cd()

    d_nyield_1 = ws.obj(f"n_1")
    d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)

    d_nyield_2 = ws.obj(f"n_2")
    d_nyield_2_err = d_nyield_2.getPropagatedError(fitresult)
    if strat == "3peak":
        d_nyield_3 = ws.obj(f"n_3")
        d_nyield_3_err = d_nyield_3.getPropagatedError(fitresult)

    d_nyield_bkg = ws.obj(f"n_bkg")
    d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

    # dtpave = ROOT.TPaveText(0.20, 0.65, 0.50, 0.95, "NB NDC")
    dtpave = ROOT.TPaveText(0.70, 0.65, 0.90, 0.85, "NB NDC")

    dtpave.SetFillStyle(0)
    dtpave.AddText(f"Yield Right: {round(d_nyield_1.getValV(),3)} #pm {round(d_nyield_1_err,3)}")
    dtpave.AddText(f"Yield Middle: {round(d_nyield_2.getValV(),3)} #pm {round(d_nyield_2_err,3)}")
    if strat == "3peak":
        dtpave.AddText(f"Yield Left: {round(d_nyield_3.getValV(),3)} #pm {round(d_nyield_3_err,3)}")
    dtpave.AddText(f"Yield BKG: {round(d_nyield_bkg.getValV(),3)} #pm {round(d_nyield_bkg_err,3)}")
    dtpave.Draw()
    frame.addObject(dtpave)

    frame.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K*0 [MeV]")
    frame.Draw()

    save_png(p, f"sb_fit_tests", f"{sc.spec}_sbfit_{strat}", rpflag = 0)

def build_d_window_ws(sc, type):

    files = glob.glob(f"charmless_files/{sc.spec}_2*.root")

    tc_sig = ROOT.TChain(f"DecayTreeTuple_SIG")
    tc_sb = ROOT.TChain(f"DecayTreeTuple_SB")

    for file_name in files:
        tc_sig.Add(file_name)
        tc_sb.Add(file_name)


    rdf_sig_base = RDF(tc_sig)
    rdf_sb_base = RDF(tc_sb)

    ROOT.gStyle.SetOptStat("ne")
    ROOT.gStyle.SetOptTitle(ROOT.kTRUE);
    ROOT.gStyle.SetPalette(ROOT.kBird)

    d1_mstart, d1_std, d2_mstart, d2_std = get_dwindow_values(spec, d1_flag, d2_flag, dst_flag, rflag = "print")

    d1_sig_max = 2*d1_std
    d2_sig_max = 2*d2_std

    d1_sb_min = 3*d1_std
    d2_sb_min = 3*d2_std

    d1_sb_max = 5*d1_std
    d2_sb_max = 5*d2_std

    d1_sb_min_line = f"(abs({d1_mstart} - D1_M) > {d1_sb_min})"
    d2_sb_min_line = f"(abs({d2_mstart} - D2_M) > {d2_sb_min})"

    d1_sb_max_line = f"(abs({d1_mstart} - D1_M) < {d1_sb_max})"
    d2_sb_max_line = f"(abs({d2_mstart} - D2_M) < {d2_sb_max})"

    d1_sb_line = f"({d1_sb_min_line} && {d1_sb_max_line})"
    d2_sb_line = f"({d2_sb_min_line} && {d2_sb_max_line})"

    d1_sig_line = f"(abs({d1_mstart} - D1_M) < {d1_sig_max})"
    d2_sig_line = f"(abs({d2_mstart} - D2_M) < {d2_sig_max})"


    for name, rdf in zip(["sig", "sb"], [rdf_sig_base, rdf_sb_base]):

        d1xbins = 100
        d2ybins = 100
        d1xmin = d1_mstart-d1_sb_max
        d1xmax = d1_mstart+d1_sb_max
        d2ymin = d2_mstart-d2_sb_max
        d2ymax = d2_mstart+d2_sb_max

        xbinw = (d1xmax - d1xmin) / d1xbins
        ybinw = (d1xmax - d1xmin) / d1xbins

        hist_all_sb = rdf.Histo2D((f"d1d2_{name}", f"d1d2_{name}", d1xbins, d1xmin, d1xmax, d2ybins, d2ymin, d2ymax), "D1_M", "D2_M")
        hist_bm = rdf.Histo1D((f"bm_{name}", f"bm_{name}", 100, 4800, 5600), 'B_M')
        hist_bdtfm = rdf.Histo1D((f"bdtfm_{name}", f"bdtfm_{name}", 100, 4800, 5600), 'B_DTF_M')

        c1 = ROOT.TCanvas("c1","c1")
        c1.SetRightMargin(0.2)
        # c1.SetLeftMargin(0.9)
        # c1.SetTopMargin(0.1)
        # c1.SetBottomMargin(0.9)

        hist_all_sb.GetXaxis().SetTitle(f"{sc.d1_string} [MeV]")
        hist_all_sb.GetYaxis().SetTitle(f"{sc.d2_string} [MeV]")
        hist_all_sb.GetZaxis().SetTitle(f"Candidates / {xbinw*ybinw:.2f} [MeV^{{2}}]")
        hist_all_sb.GetYaxis().SetTitleOffset(1)
        hist_all_sb.GetZaxis().SetTitleOffset(0.75)
        hist_all_sb.Draw("COLZ")
        save_png(c1, f"charmless_plots_2d/{sc.spec}", f"{sc.spec}_hist_all_{name}", rpflag = 0)

        c2 = ROOT.TCanvas("c2","c2")
        hist_bdtfm.SetTitle(f"{sc.b_string}")
        hist_bdtfm.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} DTF [MeV]")
        hist_bdtfm.Draw()
        save_png(c2, f"charmless_plots_dtf/{sc.spec}", f"{sc.spec}_hist_bm_DTF_{name}", rpflag = 0)

        c3 = ROOT.TCanvas("c3","c3")
        hist_bm.SetTitle(f"{sc.b_string}")
        hist_bm.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} [MeV]")
        hist_bm.Draw()
        save_png(c3, f"charmless_plots_m/{sc.spec}", f"{sc.spec}_hist_bm_{name}", rpflag = 0)

build_d_window_ws(z_c_01, type = "mc")

# files = glob.glob(f"charmless_files/{sc.spec}_sb2.root")
# getBestCand(files[0], sc.spec, "Data")
# fit_est(sc, strat = "2peak")
# fit_est(sc, strat = "3peak")
