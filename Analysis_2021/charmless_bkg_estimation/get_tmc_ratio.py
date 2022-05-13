import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
import rootutils
# from Analysis_2021.essentials import get_shapes_bkg, get_free_shapes, save_png, MakeSWeights
import glob as glob
# from Analysis_2021 import rootutils
from uncertainties import ufloat
from uncertainties import covariance_matrix, correlation_matrix
from uncertainties import unumpy, ufloat_fromstr
import math as math
import pandas
import numpy as np
RDF = ROOT.ROOT.RDataFrame

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

def get_tmc_ratio(path, mc_spec, data_ws, plot_flag = False, dict_flag = False, snap_sideband_flag = False, fit_sideband_flag = True):

    sig_files = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/sig_d/{path}.root")
    sb_files = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/sb_d/{path}.root")

    tc_sig = ROOT.TChain(f"DecayTreeTuple")
    tc_sb = ROOT.TChain(f"DecayTreeTuple")

    for file in sig_files:
        tc_sig.Add(file)

    for file in sb_files:
        tc_sb.Add(file)

    rdf_sig = RDF(tc_sig)
    rdf_sb = RDF(tc_sb)

    truth_base_string = build_truth_base(B0_ID, Dp_ID, Dp_ID)
    # truth_mom_string = build_truth_mom(B0_ID, Dp_ID, Dp_ID)

    # print(truth_base_string)
    # print(truth_mom_string)
    #({truth_mom_string}) &&
    # ({truth_mom_string}) &&

    rdf_tmc_sig = rdf_sig.Filter(f"({truth_base_string})", "truth_all")
    rdf_tmc_sb = rdf_sb.Filter(f"({truth_base_string})", "truth_all")

    tmc_sig_count = rdf_tmc_sig.Count().GetValue()
    tmc_sb_count = rdf_tmc_sb.Count().GetValue()
    mean_start = rdf_tmc_sb.Mean("B_M").GetValue()

    if snap_sideband_flag:
        clist = rdf_tmc_sb.GetColumnNames()
        rdf_tmc_sb_snap = rdf_tmc_sb.Snapshot(f"DecayTreeTuple", f"sb_tmc1_files/{mc_spec}.root", clist)
        print(f"finished snapshot for tmc1 mc")

    if fit_sideband_flag:

        sbws = ROOT.RooWorkspace("sbws")
        bmin = 5000
        bmax = 5400
        sbws.factory(f"B_M[{bmin},{bmax}]")
        b_m = sbws.var("B_M")
        data_args = ROOT.RooArgSet(b_m)
        sbfile = ROOT.TFile(f"sb_tmc1_files/{mc_spec}.root")
        sbtree = sbfile.Get("DecayTreeTuple")
        lmbdata = ROOT.RooDataSet(f"{mc_spec}_lmbsbmc", f"{mc_spec}_lmbsbmc", sbtree, data_args)

        sbws.Import(lmbdata)
        sbws.factory(f"RooKeysPdf::{mc_spec}_lmb_fit(B_M, {mc_spec}_lmbsbmc)")
        lmbfit = sbws.pdf(f"{mc_spec}_lmb_fit")
        data_ws.Import(lmbfit)

    if plot_flag:

        hist_bm_sig = rdf_tmc_sig.Histo1D((f"bm_sig", f"bm_sig", 100, 4800, 5600), 'B_M')
        hist_bm_sb = rdf_tmc_sb.Histo1D((f"bm_sb", f"bm_sb", 100, 4800, 5600), 'B_M')

        c1 = ROOT.TCanvas("c1","c1")
        hist_bm_sig.Draw()
        save_png(c1, f"test/{path}", f"{path}_sig", rpflag = 0)

        c2 = ROOT.TCanvas("c2","c2")
        hist_bm_sb.Draw()
        save_png(c2, f"test/{path}", f"{path}_sb", rpflag = 0)

        frame = b_m.frame(ROOT.RooFit.Title(f"{mc_spec}_mcsb_spectrum"), ROOT.RooFit.Bins(40))

        lmbdata.plotOn(frame, ROOT.RooFit.Name("data"))
        lmbfit.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))

        p = rootutils.residualPlot()
        p.pt.cd()
        xaxis = frame.GetXaxis()
        xaxis.SetTickLength(0)
        xaxis.SetNdivisions(0)
        xaxis.SetLabelSize(0)

        frame.Draw()
        p.pb.cd()

        hpull = frame.pullHist(f"data", f"pdf")
        pull = b_m.frame()
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

        # pull.GetXaxis().SetTitle(f" m({title} {varoi}) [MeV]")
        # frame.GetXaxis().SetTitle(f" m({title}) [MeV]")

        save_png(p, f"sb_keys_fit_tests", f"{mc_spec}_sb", rpflag = 1)
    if dict_flag:
        dict = {
            "Scheme ID" : path,
            "TMC in Sig": tmc_sig_count,
            "TMC in sb": tmc_sb_count,
        }

        folder_path = '/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/charmless_bkg_estimation/charmless_ratio'
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)

        df = pandas.DataFrame(dict , index=[0])
        df.to_csv(f"{folder_path}/{path}_charm_txt")

    return (tmc_sb_count/tmc_sig_count, mean_start)

def get_shape(dws, data_spec, mc_spec, run_name, varoi):

    nn_ws_base = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/{data_spec}_{run_name}_{varoi}.root")
    nn_ws = nn_ws_base.Get(run_name)

    if mc_spec == "Z_m_p_01":
        file = "01_Z_m_p_11198006"
    if mc_spec == "Z_m_p_0203":
        file = "02_Z_m_p_11198400"
    # Get Signal Yield
    y = nn_ws.var(f"{mc_spec}_yield")
    # Get TMC sb ratio
    tmc_ratio, mean_lmb = get_tmc_ratio(file, mc_spec, dws)
    # Get mean for charmless
    mean_charmless = nn_ws.var(f"mean_{mc_spec}")
    dws.Import(mean_charmless)
    lmb_yield = tmc_ratio*y.getValV()
    dws.factory(f"{mc_spec}_lmb_yield[{lmb_yield}]")
    dws.factory(f"Gaussian::{mc_spec}_charmless_fit(B_M, mean_{mc_spec}, width_{mc_spec}_charmless[5.0,0.1,15.0])")


def build_charmless_fit(run_name, data_spec):

    dws = ROOT.RooWorkspace(run_name)
    bmin = 5000
    bmax = 5400
    dws.factory(f"B_M[{bmin},{bmax}]")

    get_shape(dws, data_spec, "Z_m_p_01", run_name, "B_M")
    get_shape(dws, data_spec, "Z_m_p_0203", run_name, "B_M")
    get_shapes_bkg(data_spec, "Exponential", dws, varoi = "B_M")
    # get_shapes_bkg(data_spec, "Bernstein", dws, varoi = "B_M")

    dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_charmless_yield[500,0,10000]*Z_m_p_01_charmless_fit, \
                                             Z_m_p_01_lmb_yield*Z_m_p_01_lmb_fit, \
                                             Z_m_p_0203_charmless_yield[500,0,10000]*Z_m_p_0203_charmless_fit, \
                                             Z_m_p_0203_lmb_yield*Z_m_p_0203_lmb_fit, \
                                             Z_m_p_bkg_yield[100,0,100000]*Z_m_p_spectrum_bkg)")

    dws.Print()
    b_m = dws.var("B_M")
    data_args = ROOT.RooArgSet(b_m)

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{data_spec}.root")
    tchain = ROOT.TChain("DecayTreeTuple_SB")
    for file_name in file_list:
        tchain.Add(file_name)
    #
    data = ROOT.RooDataSet(f"{data_spec}_final_data", f"{data_spec}_final_data", tchain, data_args)
    model = dws.pdf(f"{data_spec}_spectrum_all_fit")
    fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    dws.Import(fit)
    fitresult = dws.obj(f"fitresult_{data_spec}_spectrum_all_fit_{data_spec}_final_data")

    frame = b_m.frame(ROOT.RooFit.Title(f"{data_spec}_spectrum"), ROOT.RooFit.Bins(40))

    ylist = ["Z_m_p_01_charmless_fit","Z_m_p_01_lmb_fit","Z_m_p_0203_charmless_fit","Z_m_p_0203_lmb_fit"]

    data.plotOn(frame, ROOT.RooFit.Name("data"))
    model.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
    model.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))

    color_list = [ROOT.kBlue, ROOT.kViolet, ROOT.kAzure, ROOT.kMagenta, ROOT.kCyan]

    for y, color in zip(ylist, color_list):
        model.plotOn(frame, ROOT.RooFit.Components(y), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(color), ROOT.RooFit.Name(f"fit_{y}"))


    legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.905)
    name_list =  ["01_charmless_fit","01_lmb_fit","0203_charmless_fit","0203_lmb_fit"]

    legend.SetFillStyle(1001)

    d_nyield_bkg = dws.obj(f"{data_spec}_bkg_yield")
    d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

    legend.AddEntry(frame.findObject("data"),"Run 2 Data","ep")
    legend.AddEntry(frame.findObject("pdf"),"Total Fit PDF","l")
    legend.AddEntry(frame.findObject("bkg"),f"Background PDF : {round(d_nyield_bkg.getValV(),1)} #pm {round(d_nyield_bkg_err,1)}","l")
    for y, name in zip(ylist, name_list):
        yname = y.split("_fit")[0]
        print(yname)
        d_nyield = dws.obj(f"{yname}_yield")
        d_nyield_err = d_nyield.getPropagatedError(fitresult)
        legend.AddEntry(frame.findObject(f"fit_{y}"), f"{name}: {round(d_nyield.getValV(),1)} #pm {round(d_nyield_err,1)}", "l")
    legend.SetTextSize(0.020)

    p = rootutils.residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)

    frame.Draw()
    legend.Draw()
    p.pb.cd()

    hpull = frame.pullHist(f"data", f"pdf")
    pull = b_m.frame()
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

    # pull.GetXaxis().SetTitle(f" m({title} {varoi}) [MeV]")
    # frame.GetXaxis().SetTitle(f" m({title}) [MeV]")

    save_pdf(p, f"fit_tests", f"{data_spec}_{run_name}_nl", rpflag = 1)

build_charmless_fit("charmlesstest", "Z_m_p")
