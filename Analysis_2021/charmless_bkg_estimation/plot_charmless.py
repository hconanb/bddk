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
import math
#####################################
PE_mcdtt = "(pow(D1_TRUEP_E + D2_TRUEP_E + KST_TRUEP_E,2))"
PX_mcdtt = "(pow(D1_TRUEP_X + D2_TRUEP_X + KST_TRUEP_X,2))"
PY_mcdtt = "(pow(D1_TRUEP_Y + D2_TRUEP_Y + KST_TRUEP_Y,2))"
PZ_mcdtt = "(pow(D1_TRUEP_Z + D2_TRUEP_Z + KST_TRUEP_Z,2))"
M2_mcdtt = f"{PE_mcdtt} - ({PX_mcdtt} + {PY_mcdtt} + {PZ_mcdtt})"
True_BM = f"sqrt({M2_mcdtt})"

PE_mcdtt_d1_t = "(pow(D1H1_TRUEP_E + D1H2_TRUEP_E + D1H3_TRUEP_E,2))"
PX_mcdtt_d1_t = "(pow(D1H1_TRUEP_X + D1H2_TRUEP_X + D1H3_TRUEP_X,2))"
PY_mcdtt_d1_t = "(pow(D1H1_TRUEP_Y + D1H2_TRUEP_Y + D1H3_TRUEP_Y,2))"
PZ_mcdtt_d1_t = "(pow(D1H1_TRUEP_Z + D1H2_TRUEP_Z + D1H3_TRUEP_Z,2))"
M2_mcdtt_d1_t = f"{PE_mcdtt_d1_t} - ({PX_mcdtt_d1_t} + {PY_mcdtt_d1_t} + {PZ_mcdtt_d1_t})"
True_D1_M_t = f"sqrt({M2_mcdtt_d1_t})"

PE_mcdtt_d1 = "(pow(D1_TRUEP_E,2))"
PX_mcdtt_d1 = "(pow(D1_TRUEP_X,2))"
PY_mcdtt_d1 = "(pow(D1_TRUEP_Y,2))"
PZ_mcdtt_d1 = "(pow(D1_TRUEP_Z,2))"
M2_mcdtt_d1 = f"{PE_mcdtt_d1} - ({PX_mcdtt_d1} + {PY_mcdtt_d1} + {PZ_mcdtt_d1})"
True_D1_M = f"sqrt({M2_mcdtt_d1})"
#####################################
#####################################
PE_mcdtt_d2_t = "(pow(D2H1_TRUEP_E + D2H2_TRUEP_E + D2H3_TRUEP_E,2))"
PX_mcdtt_d2_t = "(pow(D2H1_TRUEP_X + D2H2_TRUEP_X + D2H3_TRUEP_X,2))"
PY_mcdtt_d2_t = "(pow(D2H1_TRUEP_Y + D2H2_TRUEP_Y + D2H3_TRUEP_Y,2))"
PZ_mcdtt_d2_t = "(pow(D2H1_TRUEP_Z + D2H2_TRUEP_Z + D2H3_TRUEP_Z,2))"
M2_mcdtt_d2_t = f"{PE_mcdtt_d2_t} - ({PX_mcdtt_d2_t} + {PY_mcdtt_d2_t} + {PZ_mcdtt_d2_t})"
True_D2_M_t = f"sqrt({M2_mcdtt_d2_t})"

PE_mcdtt_d2 = "(pow(D2_TRUEP_E,2))"
PX_mcdtt_d2 = "(pow(D2_TRUEP_X,2))"
PY_mcdtt_d2 = "(pow(D2_TRUEP_Y,2))"
PZ_mcdtt_d2 = "(pow(D2_TRUEP_Z,2))"
M2_mcdtt_d2 = f"{PE_mcdtt_d2} - ({PX_mcdtt_d2} + {PY_mcdtt_d2} + {PZ_mcdtt_d2})"
True_D2_M = f"sqrt({M2_mcdtt_d2})"
#############################################################
#############################################################
PE_mcdtt_kst_t = "(pow(KSTH1_TRUEP_E + KSTH2_TRUEP_E,2))"
PX_mcdtt_kst_t = "(pow(KSTH1_TRUEP_X + KSTH2_TRUEP_X,2))"
PY_mcdtt_kst_t = "(pow(KSTH1_TRUEP_Y + KSTH2_TRUEP_Y,2))"
PZ_mcdtt_kst_t = "(pow(KSTH1_TRUEP_Z + KSTH2_TRUEP_Z,2))"
M2_mcdtt_kst_t = f"{PE_mcdtt_kst_t} - ({PX_mcdtt_kst_t} + {PY_mcdtt_kst_t} + {PZ_mcdtt_kst_t})"
True_KST_M_t = f"sqrt({M2_mcdtt_kst_t})"
#
PE_mcdtt_kst = "(pow(KST_TRUEP_E,2))"
PX_mcdtt_kst = "(pow(KST_TRUEP_X,2))"
PY_mcdtt_kst = "(pow(KST_TRUEP_Y,2))"
PZ_mcdtt_kst = "(pow(KST_TRUEP_Z,2))"
M2_mcdtt_kst = f"{PE_mcdtt_kst} - ({PX_mcdtt_kst} + {PY_mcdtt_kst} + {PZ_mcdtt_kst})"
True_KST_M = f"sqrt({M2_mcdtt_kst})"

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

B0_ID = 511
Bp_ID = 521
Bs_ID = 531

Dp_ID = 411
D0_ID = 421
Ds_ID = 431

Dpst_ID = 413
D0st_ID = 423
Dsst_ID = 433

k_ID = 321
pi_ID = 211
kst0_ID = 313

z_c = spectrum_class("Z_m_p", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
# zz_c = spectrum_class("Z_z_z", "#barD^{0}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D^{0} K*0")
# p_c = spectrum_class("P_z_p", "#bar{D^{0}}", "D^{+}", d0mass, dpmass, "B #rightarrow #bar{D^{0}} D+ K*0")
# m_c = spectrum_class("M_m_z", "D^{-}", "D^{0}", dpmass, d0mass,  "B #rightarrow D- D^{0} K*0")
z_c_01 = spectrum_class("01_Z_m_p_11198006", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
z_c_02 = spectrum_class("02_Z_m_p_11198400", "D^{*-}", "D^{+}", dpmass, dpmass, "B #rightarrow D*-D+K*0")
z_c_04 = spectrum_class("04_Z_m_p_11198401", "D^{*-}", "D^{*+}", dpmass, dpmass, "B #rightarrow D*-D*+K*0")
z_c_09 = spectrum_class("09_Z_z_z_11196", "D^{*-}", "D^{*+}", dpmass, dpmass, "B #rightarrow D*-D*+K*0")

# st_c = spectrum_class("P_z_pst","#bar{D^{0}}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D*+ K*0")
# # s_c = spectrum_class("s","D^{-}_{s}", "D^{+}" ,dsmass, dpmass, "e_g")
# # n7_c = spectrum_class("norm7","#bar{D^{0}}", "D^{0} #rightarrow k#pi#pi#pi", d0mass, d0mass, "e_dg_b")
# n8_c = spectrum_class("norm8","D^{-}", "D^{0} #rightarrow k#pi#pi#pi", dpmass, d0mass, "B #rightarrow D-D+K*0")
def build_truth_base(b_spec, c1_spec, c2_spec, norm_flag = 0, st_flag = 0):
    truth_MC = f"abs(B_TRUEID) == {b_spec}"
    truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}"
    truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID}"
    if c1_spec == Dp_ID:
        truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID}"
    truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}"
    truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID}"
    if c2_spec == Dp_ID:
        truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID}"
    if norm_flag == 1:
        truth_MC = f"{truth_MC} && abs(D2H4_TRUEID) == {pi_ID}"
        truth_MC = f"{truth_MC} && abs(K_TRUEID) == {k_ID}"
    if norm_flag == 0:
        truth_MC = f"{truth_MC} && abs(KSTH1_TRUEID) == {k_ID}"
        truth_MC = f"{truth_MC} && abs(KSTH2_TRUEID) == {pi_ID}"
    return truth_MC

def build_truth_mom(b_spec, c1_spec, c2_spec, norm_flag = 0, st_flag = 0):
    truth_MC = f"abs(D1_TRUEID) == {c1_spec}"
    truth_MC = f"{truth_MC} && abs(D2_TRUEID) == {c2_spec}"
    if norm_flag == 0:
        truth_MC = f"{truth_MC} && abs(KST_TRUEID) == {kst0_ID}"
    return truth_MC

def fit_est(sc, strat = "2peak"):

    files = glob.glob(f"charmless_files/{sc.spec}_2*.root")

    tc_sig = ROOT.TChain(f"DecayTreeTuple_SIG")
    tc_sb = ROOT.TChain(f"DecayTreeTuple_SB")

    for file_name in files:
        tc_sig.Add(file_name)
        tc_sb.Add(file_name)

    ws = ROOT.RooWorkspace(f"fit_{sc.spec}")

    if strat == "3peak":
        ws.factory("B_M[4800, 5600]")
    if strat == "2peak":
        ws.factory("B_M[5000, 5400]")


    ws.factory(f"Gaussian::{sc.spec}_mpeak_fit_1(B_M, mean_mpeak[5100,5075,5125], width_mpeak[5.0,0.1,100.0])")

    ws.factory(f"Exponential:{sc.spec}_sig_spectrum_bkg(B_M, c0_{sc.spec}_sig[0, -2, 2])")
    # ws.factory(f"Chebychev:{sc.spec}_sig_spectrum_bkg(B_M,{{c0_{sc.spec}_sig[0.,-3,3],c1_{sc.spec}_sig[0.,-3,3], c2_{sc.spec}_sig[0.,-3,3], c3_{sc.spec}_sig[0.,-3,3]}})")
    ws.factory(f"Gaussian::{sc.spec}_sig_fit_1(B_M, mean_1[5280,5255,5295], width_1[5.0,0.1,50.0])")
    ws.factory(f"Gaussian::{sc.spec}_sig_fit_2(B_M, mean_2[5100,5075,5150], width_2[5.0,0.1,50.0])")
    ws.factory(f"SUM::{sc.spec}_sig_all_fit(n_1_sig[100,0,10000]*{sc.spec}_sig_fit_1,n_2_sig[100,0,10000]*{sc.spec}_sig_fit_2,n_bkg_sig[100,0,10000]*{sc.spec}_sig_spectrum_bkg)")

    ws.factory(f"Exponential:{sc.spec}_sb_spectrum_bkg(B_M, c0_{sc.spec}_sb[0, -2, 2])")

    # ws.factory(f"Chebychev:{sc.spec}_sb_spectrum_bkg(B_M,{{c0_{sc.spec}_sb[0.,-3,3],c1_{sc.spec}_sb[0.,-3,3], c2_{sc.spec}_sb[0.,-3,3], c3_{sc.spec}_sb[0.,-3,3]}})")
    ws.factory(f"Gaussian::{sc.spec}_sb_fit_1(B_M, mean_1, width_1)")
    ws.factory(f"Gaussian::{sc.spec}_sb_fit_2(B_M, mean_2, width_2)")
    ws.factory(f"SUM::{sc.spec}_sb_all_fit(n_mpeak[100,0,10000]*{sc.spec}_mpeak_fit_1, n_1_sb[100,0,10000]*{sc.spec}_sb_fit_1,n_2_sb[100,0,10000]*{sc.spec}_sb_fit_2,n_bkg_sb[100,0,10000]*{sc.spec}_sb_spectrum_bkg)")

    # ws.factory(f"Gaussian::{sc.spec}_fit_3(B_M, mean_3[4900,4850,4950], width_3[5.0,0.1,50.0])")

    # if strat == "3peak":
    #     ws.factory(f"SUM::{sc.spec}_all_fit(n_1[100,0,10000]*{sc.spec}_fit_1,n_2[100,0,10000]*{sc.spec}_fit_2,n_3[100,0,10000]*{sc.spec}_fit_3,n_bkg[100,0,10000]*{sc.spec}_spectrum_bkg)")
    # if strat == "2peak":
        # ws.factory(f"SUM::{sc.spec}_all_fit(n_1[100,0,10000]*{sc.spec}_fit_1,n_2[100,0,10000]*{sc.spec}_fit_2,n_bkg[100,0,10000]*{sc.spec}_spectrum_bkg)")

    b_m = ws.var("B_M")
    data_args = ROOT.RooArgSet(b_m)

    ##########
    # files_mpeak = glob.glob(f"charmless_files/02_Z_m_p*.root")
    # tc_sb_mpeak = ROOT.TChain(f"DecayTreeTuple_SB")
    #
    # for file_name_mp in files_mpeak:
    #     tc_sb_mpeak.Add(file_name_mp)
    #     print(file_name_mp)

    # rdf_mpeak = RDF(tc_sb_mpeak)
    # b_m.setBins(10)
    # rooDataSet = rdf_mpeak.Book(ROOT.std.move(ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(b_m))), ("B_M"))
    #

    # Method 2:
    # We first declare the RooDataHistMaker
    # rdhMaker = ROOT.RooDataSetHelper("dataset", "Title of dataset", ROOT.RooArgSet(b_m))
    # Then, we move it into the RDataFrame action:
    # data_set_mpeak = rdf_mpeak.Book(ROOT.std.move(rdhMaker), ("B_M"))

    # hist_bm = rdf_mpeak.Histo1D((f"bm", f"bm", 100, 4800, 5600), 'B_M')
    # ws.Import(hist_bm)
    # hist = ws.obj("bm")
    # # hist_bm_copy = hist_bm.GetPtr()
    # # print(hist_bm_copy)
    # data_set_mpeak = ROOT.RooDataHist("data","dataset with x", data_args, ROOT.RooFit.Import("A",hist))

    # mpeak_fit = ws.pdf(f"{sc.spec}_mpeak_fit_1")
    # data_set_mpeak = ROOT.RooDataSet(f"{sc.spec}_data_mpeak", f"{sc.spec}_data_mpeak", tc_sb_mpeak, data_args)
    # all_fit_mpeak = mpeak_fit.fitTo(data_set_mpeak, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    #
    # frame = b_m.frame(ROOT.RooFit.Title(f"{sc.spec}"))
    # data_set_mpeak.plotOn(frame, ROOT.RooFit.Name("data"))
    # # mpeak_fit.plotOn(frame,  ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))
    #
    # p = ROOT.TCanvas("p1","p1")
    # p.cd()
    #
    # frame.GetXaxis().SetTitle(f"mpeak {sc.d1_string}{sc.d2_string}K*0 [MeV]")
    # frame.Draw()
    #
    # save_png(p, f"sb_fit_tests", f"{sc.spec}_mpeak_fit_{strat}", rpflag = 0)

    ##################

    sig_fit = ws.pdf(f"{sc.spec}_sig_all_fit")
    data_set_sig = ROOT.RooDataSet(f"{sc.spec}_data_sig", f"{sc.spec}_data_sig", tc_sig, data_args)
    all_fit_sig = sig_fit.fitTo(data_set_sig, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    ws.Import(all_fit_sig)
    fitresult_sig = ws.obj(f"fitresult_{sc.spec}_sig_all_fit_{sc.spec}_data_sig")

    vars = ws.allVars()
    for i in vars:
        if i.GetName() in [f"mean_1",f"mean_2","width_1","width_2"]:
            temp = ws.var(i.GetName())
            temp.setConstant(True)
            print(f"{temp} is constant")

    data_set_sb = ROOT.RooDataSet(f"{sc.spec}_data_sb", f"{sc.spec}_data_sb", tc_sb, data_args)
    sb_fit = ws.pdf(f"{sc.spec}_sb_all_fit")
    all_fit_sb = sb_fit.fitTo(data_set_sb, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    ws.Import(all_fit_sb)
    fitresult_sb = ws.obj(f"fitresult_{sc.spec}_sb_all_fit_{sc.spec}_data_sb")

    for data_set, fit, fitresult, name in zip([data_set_sig, data_set_sb], [sig_fit, sb_fit], [fitresult_sig, fitresult_sb], ["sig","sb"]):

        frame = b_m.frame(ROOT.RooFit.Title(f"{sc.spec}"))
        data_set.plotOn(frame, ROOT.RooFit.Name("data"))
        fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_{name}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))
        fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_{name}_all_fit"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
        fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_{name}_fit_1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("f1"))
        fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_{name}_fit_2"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet), ROOT.RooFit.Name("f2"))
        if name == "sb":
            fit.plotOn(frame, ROOT.RooFit.Components(f"{sc.spec}_mpeak_fit_1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kMagenta), ROOT.RooFit.Name("f3"))

        p = ROOT.TCanvas("p1","p1")
        p.cd()

        d_nyield_1 = ws.obj(f"n_1_{name}")
        d_nyield_1_err = d_nyield_1.getPropagatedError(fitresult)

        d_nyield_2 = ws.obj(f"n_2_{name}")
        d_nyield_2_err = d_nyield_2.getPropagatedError(fitresult)
        # if strat == "3peak":
        #     d_nyield_3 = ws.obj(f"n_3_{name}")
        #     d_nyield_3_err = d_nyield_3.getPropagatedError(fitresult)

        d_nyield_bkg = ws.obj(f"n_bkg_{name}")
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

        save_png(p, f"sb_fit_tests", f"{sc.spec}_{name}_fit_{strat}", rpflag = 0)

def build_d_window_ws(sc, type):
    if type == "data":
        files = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/*/final_sample/{sc.spec}.root")
    if type == "mc":
        files = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/*/final_sample/{sc.spec}.root")

    tc_sig = ROOT.TChain(f"DecayTreeTuple_SIG")
    tc_sb = ROOT.TChain(f"DecayTreeTuple_SB")

    for file_name in files:
        tc_sig.Add(file_name)
        tc_sb.Add(file_name)

    rdf_sig_base = RDF(tc_sig)
    rdf_sb_base = RDF(tc_sb)

    if type == "mc":
        truth_base_string = build_truth_base(B0_ID, Dp_ID, Dp_ID)
        print(truth_base_string)
        rdf_sig_base = rdf_sig_base.Filter(truth_base_string, "truth_base_1")
        rdf_sb_base = rdf_sb_base.Filter(truth_base_string, "truth_base_2")

    ROOT.gStyle.SetOptStat("ne")
    ROOT.gStyle.SetOptTitle(ROOT.kTRUE);
    ROOT.gStyle.SetPalette(ROOT.kBird)

    d1_flag = "mp"
    d2_flag = "mp"
    dst_flag = False

    d1_mstart, d1_std, d2_mstart, d2_std = get_dwindow_values(sc.spec, d1_flag, d2_flag, dst_flag, rflag = "print")

    d1_sig_max = 2*d1_std
    d2_sig_max = 2*d2_std

    d1_sb_min = 3*d1_std
    d2_sb_min = 3*d2_std

    d1_sb_max = 5*d1_std
    d2_sb_max = 5*d2_std

    rdf_sig_base.Define("D1_M_Plot",f"D1_M - {d1_mstart}") \
                .Define("D2_M_Plot",f"D2_M - {d2_mstart}")

    rdf_sb_base.Define("D1_M_Plot",f"D1_M - {d1_mstart}") \
               .Define("D2_M_Plot",f"D2_M - {d2_mstart}")

    for name, rdf in zip(["sig", "sb"], [rdf_sig_base, rdf_sb_base]):
        # for fd_cut in [0, 0.5, 1, 1.5, 2]:

        d1xbins = 100
        d2ybins = 100
        d1xmin = d1_mstart-d1_sb_max
        d1xmax = d1_mstart+d1_sb_max
        d2ymin = d2_mstart-d2_sb_max
        d2ymax = d2_mstart+d2_sb_max

        xbinw = (d1xmax - d1xmin) / d1xbins
        ybinw = (d1xmax - d1xmin) / d1xbins

        bins_cust = 6
        d1_binw_array = [-5*d1_std, -3*d1_std, 2*d1_std, 0, 2*d1_std, 3*d1_std, 5*d1_std,]

        # rdf = rdf.Filter(f"D1_FDCHI2_ORIVX > {fd_cut} && D2_FDCHI2_ORIVX > {fd_cut}")
        hist_all_sb = rdf.Histo2D((f"d1d2_{name}", f"d1d2_{name}", d1xbins, d1xmin, d1xmax, d2ybins, d2ymin, d2ymax), "D1_M", "D2_M")
        hist_bm = rdf.Histo1D((f"bm_{name}", f"bm_{name}", 100, 5000, 5600), 'B_M')
            # hist_bdtfm = rdf.Histo1D((f"bdtfm_{name}", f"bdtfm_{name}", 100, 4800, 5600), 'B_DTF_M')


            # if type == "mc":
            #     rdf_0 = rdf.Filter("B_BKGCAT == 0")
            #     rdf_50 = rdf.Filter("B_BKGCAT == 50")
            #     rdf_60 = rdf.Filter("B_BKGCAT == 60")
            #
            #     for bkg_name, rdf_bkg in zip(["sig", "lowmassbkg", "ghost"], [rdf_0, rdf_50, rdf_60]):
            #
            #         hist_all_bkgcat = rdf_bkg.Histo2D((f"d1d2_{bkg_name}", f"d1d2_{bkg_name}", d1xbins, d1xmin, d1xmax, d2ybins, d2ymin, d2ymax), "D1_M", "D2_M")
            #         hist_bm_bkgcat =  rdf_bkg.Histo1D((f"bm_{bkg_name}", f"bm_{bkg_name}", 100, 4800, 5600), 'B_M')
            #         hist_bdtfm_bkgcat =  rdf_bkg.Histo1D((f"bdtfm_{bkg_name}", f"bdtfm_{bkg_name}", 100, 4800, 5600), 'B_DTF_M')
            #
            #         c1 = ROOT.TCanvas("c1","c1")
            #         c1.SetRightMargin(0.2)
            #         # c1.SetLeftMargin(0.9)
            #         # c1.SetTopMargin(0.1)
            #         # c1.SetBottomMargin(0.9)
            #
            #         xax = hist_all_bkgcat.GetXaxis()
            #         yax = hist_all_bkgcat.GetYaxis()
            #         xax.SetTitle(f"{sc.d1_string} [MeV]")
            #         yax.SetTitle(f"{sc.d2_string} [MeV]")
            #
            #         xax.SetNdivisions(-10)
            #         yax.SetNdivisions(-10)
            #         # label_list = ["-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma"]
            #
            #         for i in range(1, 12):
            #             xax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")
            #             yax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")
            #
            #
            #         hist_all_bkgcat.GetZaxis().SetTitle(f"Candidates / {xbinw*ybinw:.2f} [MeV^{{2}}]")
            #         hist_all_bkgcat.GetYaxis().SetTitleOffset(1)
            #         hist_all_bkgcat.GetZaxis().SetTitleOffset(0.75)
            #         hist_all_bkgcat.Draw("COLZ")
            #
            #         save_png(c1, f"charmless_plots_bkg/{sc.spec}", f"{sc.spec}_hist_all_{name}_{bkg_name}", rpflag = 0)
            #
            #         c2 = ROOT.TCanvas("c2","c2")
            #         hist_bdtfm_bkgcat.SetTitle(f"{sc.b_string}")
            #         hist_bdtfm_bkgcat.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} DTF [MeV]")
            #         hist_bdtfm_bkgcat.Draw()
            #         save_png(c2, f"charmless_plots_bkg_dtf/{sc.spec}", f"{sc.spec}_hist_bm_DTF_{name}_{bkg_name}", rpflag = 0)
            #
            #         c3 = ROOT.TCanvas("c3","c3")
            #         hist_bm_bkgcat.SetTitle(f"{sc.b_string}")
            #         hist_bm_bkgcat.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} [MeV]")
            #         hist_bm_bkgcat.Draw()
            #         save_png(c3, f"charmless_plots_bkg_m/{sc.spec}", f"{sc.spec}_hist_bm_{name}_{bkg_name}", rpflag = 0)
            #
            #
            #
        c1 = ROOT.TCanvas("c1","c1")
        c1.SetRightMargin(0.2)
        # c1.SetLeftMargin(0.9)
        # c1.SetTopMargin(0.1)
        # c1.SetBottomMargin(0.9)

        xax = hist_all_sb.GetXaxis()
        yax = hist_all_sb.GetYaxis()
        xax.SetTitle(f"{sc.d1_string} [MeV]")
        yax.SetTitle(f"{sc.d2_string} [MeV]")

        xax.SetNdivisions(-10)
        yax.SetNdivisions(-10)
        # label_list = ["-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma","-5#sigma"]

        for i in range(1, 12):
            xax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")
            yax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")


        hist_all_sb.GetZaxis().SetTitle(f"Candidates / {xbinw*ybinw:.2f} [MeV^{{2}}]")
        hist_all_sb.GetYaxis().SetTitleOffset(1)
        hist_all_sb.GetZaxis().SetTitleOffset(0.75)
        hist_all_sb.Draw("COLZ")
        save_png(c1, f"charmless_plots/{sc.spec}", f"{sc.spec}_{name}_ddmass", rpflag = 0)
            #
            # c2 = ROOT.TCanvas("c2","c2")
            # hist_bdtfm.SetTitle(f"{sc.b_string}")
            # hist_bdtfm.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} DTF [MeV]")
            # hist_bdtfm.Draw()
            # save_png(c2, f"charmless_plots_dtf/{sc.spec}", f"{sc.spec}_hist_bm_DTF_{name}", rpflag = 0)

        c3 = ROOT.TCanvas("c3","c3")
        hist_bm.SetTitle(f"{sc.b_string}")
        hist_bm.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} [MeV]")
        hist_bm.Draw()
        save_png(c3, f"charmless_plots/{sc.spec}", f"{sc.spec}_{name}_bmass", rpflag = 0)

def build_d_window_wmultiplecans_ws(sc, type):

    files = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/*/{sc.spec}.root")
    # files_sb = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/sb_d/{sc.spec}.root")
    # files_r23 = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/r23_d/{sc.spec}.root")

    tc_base = ROOT.TChain(f"DecayTreeTuple")

    for file in files:
        print(file)
        tc_base.Add(file)

    rdf_base = RDF(tc_base)

    ROOT.gStyle.SetOptStat("ne")
    ROOT.gStyle.SetOptTitle(ROOT.kTRUE);
    ROOT.gStyle.SetPalette(ROOT.kBird)

    d1_flag = "mp"
    d2_flag = "mp"
    dst_flag = False

    d1_mstart, d1_std, d2_mstart, d2_std = get_dwindow_values(sc.spec, d1_flag, d2_flag, dst_flag, rflag = "print")

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

    sig_sig_line = f"({d1_sig_line} && {d2_sig_line})"
    sig_sb_line = f"({d1_sig_line} && {d2_sb_line})"
    sb_sig_line = f"({d1_sb_line} && {d2_sig_line})"
    sb2_line = f"({sig_sb_line} || {sb_sig_line})"

    rdf_base = rdf_base.Define("D1_M_Plot",f"D1_M - {d1_mstart}") \
                       .Define("D2_M_Plot",f"D2_M - {d2_mstart}")

    rdf_base = rdf_base.Define("True_P1", f"sqrt(pow(D1_TRUEP_X,2) + pow(D1_TRUEP_Y,2) + pow(D1_TRUEP_Z,2))") \
                       .Define("P1_Res", f"(True_P1 - D1_P)/True_P1") \
                       .Define("True_P2", f"sqrt(pow(D2_TRUEP_X,2) + pow(D2_TRUEP_Y,2) + pow(D2_TRUEP_Z,2))") \
                       .Define("P2_Res", f"(True_P2 - D2_P)/True_P2") \
                       .Define("True_D1_M", True_D1_M) \
                       .Define("True_D2_M", True_D2_M) \
                       .Define("True_KST_M", True_KST_M) \
                       .Define("True_D1_M_t", True_D1_M_t) \
                       .Define("True_D2_M_t", True_D2_M_t) \
                       .Define("True_KST_M_t", True_KST_M_t) \
                       .Define("D1_diff", f"{True_D1_M} - {True_D1_M_t}") \
                       .Define("D2_diff", f"{True_D2_M} - {True_D2_M_t}") \
                       .Define("KST_diff", f"{True_KST_M} - {True_KST_M_t}")

    truth_base_string = build_truth_base(B0_ID, Dp_ID, Dp_ID)
    truth_mom_string = build_truth_mom(B0_ID, Dp_ID, Dp_ID)

    print(truth_base_string)
    print(truth_mom_string)

    rdf_truthbase = rdf_base.Filter(truth_base_string, "truth_base")
    rdf_mom = rdf_base.Filter(truth_mom_string, "truth_mom")
    rdf_fulltmc = rdf_base.Filter(f"({truth_mom_string}) && ({truth_base_string})", "truth_all")

#.Filter(sb2_line)
#.Filter(sb2_line)

    rdf_sb_base = rdf_base
    rdf_sb_truth = rdf_truthbase
    rdf_sb_truthall = rdf_fulltmc.Filter(sb2_line)

    clist = rdf_sb_truthall.GetColumnNames()
    # rdfsnaptt = rdf_sb_truthall.Snapshot(f"DecayTreeTuple", f"sb_alltmc_{sc.spec}.root", clist)

    d1xbins = 100
    d2ybins = 100
    d1xmin = d1_mstart-d1_sb_max
    d1xmax = d1_mstart+d1_sb_max
    d2ymin = d2_mstart-d2_sb_max
    d2ymax = d2_mstart+d2_sb_max

    xbinw = (d1xmax - d1xmin) / d1xbins
    ybinw = (d1xmax - d1xmin) / d1xbins

    from array import array

    bins_cust = 5
    d1_binw_array = array('d', [-5*d1_std, -3*d1_std, -2*d1_std, 2*d1_std, 3*d1_std, 5*d1_std])
    d2_binw_array = array('d', [-5*d2_std, -3*d2_std, -2*d2_std, 2*d2_std, 3*d2_std, 5*d2_std])
#"base", "truthbase", "mom", "fulltmc",
#rdf_base, rdf_truthbase, rdf_mom, rdf_fulltmc,
    for truth_name, rdf_bkg in zip(["tmc2_all", "tm_sb"], [rdf_fulltmc, rdf_sb_truthall]):

        # hist_all_truth = rdf_bkg.Histo2D((f"d1d2_{truth_name}", f"d1d2_{truth_name}", bins_cust, d1_binw_array, bins_cust, d2_binw_array), "D1_M_Plot", "D2_M_Plot")
        # hist_all_truth = rdf_bkg.Histo2D((f"d1d2_{truth_name}", f"d1d2_{truth_name}", d1xbins, d1xmin, d1xmax, d2ybins, d2ymin, d2ymax), "D1_M", "D2_M")
        hist_bm_bkgcat =  rdf_bkg.Histo1D((f"bm_{truth_name}", f"bm_{truth_name}", 100, 5000, 5600), 'B_M')
        # hist_bdtfm_bkgcat =  rdf_bkg.Histo1D((f"bdtfm_{truth_name}", f"bdtfm_{truth_name}", 100, 4800, 5600), 'B_DTF_M')

        # hist_p1 =  rdf_bkg.Histo1D((f"{truth_name}_p1", f"bm_{truth_name}_p1", 100, -1, 1), 'P1_Res')
        # hist_p2 =  rdf_bkg.Histo1D((f"{truth_name}_p2", f"bm_{truth_name}_p2", 100, -1, 1), 'P2_Res')
        #
        # hist_d1 =  rdf_bkg.Histo1D((f"{truth_name}_D1d", f"bm_{truth_name}_D1d", 100, -0.5, 100), 'D1_diff')
        # hist_d2 =  rdf_bkg.Histo1D((f"{truth_name}_D2d", f"bm_{truth_name}_D2d", 100, -0.5, 100), 'D2_diff')
        # hist_kst =  rdf_bkg.Histo1D((f"{truth_name}_KSTd", f"bm_{truth_name}_KSTd", 100, -0.5, 100), 'KST_diff')

#         c1 = ROOT.TCanvas("c1","c1")
#         c1.SetRightMargin(0.2)
#         # c1.SetLeftMargin(0.9)
#         # c1.SetTopMargin(0.1)
#         # c1.SetBottomMargin(0.9)
#
#         xax = hist_all_truth.GetXaxis()
#         yax = hist_all_truth.GetYaxis()
#         xax.SetTitle(f"{sc.d1_string} [MeV]")
#         yax.SetTitle(f"{sc.d2_string} [MeV]")
#
#         xax.SetNdivisions(-10)
#         yax.SetNdivisions(-10)
#
#         for i in range(1, 12):
#             xax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")
#             yax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")
#
#         hist_all_truth.GetZaxis().SetTitle(f"Candidates [MeV^{{2}}]")
#         hist_all_truth.GetYaxis().SetTitleOffset(1)
#         hist_all_truth.GetZaxis().SetTitleOffset(0.75)
#         hist_all_truth.Draw("COLZ")
#
#         save_png(c1, f"charmless_plots_truth", f"{sc.spec}_{truth_name}", rpflag = 0)
#
#         c2 = ROOT.TCanvas("c2","c2")
#         hist_bdtfm_bkgcat.SetTitle(f"{sc.b_string}")
#         hist_bdtfm_bkgcat.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} DTF [MeV]")
#         hist_bdtfm_bkgcat.Draw()
#         save_png(c2, f"charmless_plots_bkg_dtf/{sc.spec}", f"{sc.spec}_hist_bm_DTF_{truth_name}", rpflag = 0)
#
        c3 = ROOT.TCanvas("c3","c3")
        hist_bm_bkgcat.SetTitle(f"{sc.b_string}")
        hist_bm_bkgcat.GetXaxis().SetTitle(f"{sc.d1_string}{sc.d2_string}K^{{*0}} [MeV]")
        hist_bm_bkgcat.Draw()
        save_png(c3, f"charmless_plots_bkg_m/{sc.spec}", f"{sc.spec}_hist_bm_{truth_name}", rpflag = 0)

#         c4 = ROOT.TCanvas("c4","c4")
#         hist_p1.GetXaxis().SetTitle(f"{sc.d1_string} Pres")
#         hist_p1.Draw()
#         save_png(c4, f"charmless_plots_pres/{sc.spec}", f"{sc.spec}_p1res_{truth_name}", rpflag = 0)
#
#         c5 = ROOT.TCanvas("c5","c5")
#         hist_p2.GetXaxis().SetTitle(f"{sc.d1_string} Pres")
#         hist_p2.Draw()
#         save_png(c5, f"charmless_plots_pres/{sc.spec}", f"{sc.spec}_p2res_{truth_name}", rpflag = 0)

        # c6 = ROOT.TCanvas("c6","c6")
        # c6.SetLogy()
        # hist_d1.GetXaxis().SetTitle(f"{sc.d1_string} D1d")
        # hist_d1.Draw()
        # save_png(c6, f"charmless_plots_mtest/{sc.spec}", f"{sc.spec}_d1d_{truth_name}", rpflag = 0)
        #
        # c7 = ROOT.TCanvas("c7","c7")
        # c7.SetLogy()
        # hist_d2.GetXaxis().SetTitle(f"{sc.d2_string} D2d")
        # hist_d2.Draw()
        # save_png(c7, f"charmless_plots_mtest/{sc.spec}", f"{sc.spec}_d2d_{truth_name}", rpflag = 0)
        #
        # c8 = ROOT.TCanvas("c8","c8")
        # c8.SetLogy()
        # hist_kst.GetXaxis().SetTitle(f"kstd")
        # hist_kst.Draw()
        # save_png(c8, f"charmless_plots_mtest/{sc.spec}", f"{sc.spec}_kstd_{truth_name}", rpflag = 0)

# build_d_window_ws(z_c, type = "data")
build_d_window_ws(z_c_01, type = "mc")
build_d_window_ws(z_c_02, type = "mc")
# build_d_window_wmultiplecans_ws(z_c_01, type = "mc")
# build_d_window_wmultiplecans_ws(z_c_02, type = "mc")
# build_d_window_wmultiplecans_ws(zz_c_09, type = "mc")

# build_d_window_ws(z_c_04, type = "mc")
# build_d_window_ws(z_c, type = "data")
#
# fit_est(z_c, strat = "2peak")
# fit_est(sc, strat = "3peak")
