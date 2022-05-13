import sys
import os

basedir = os.getcwd().split('sys_studies')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

from hep_ml import reweight
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from hep_ml.metrics_utils import ks_2samp_weighted
import matplotlib.pyplot as plt
import warnings

def build_sb_plot(data_tchain, data_spec):

    rdf_base = RDF(data_tchain)

    d1_mstart, d2_mstart = get_dwindow_values(data_spec, "charmless_start")
    d1_sb_min, d2_sb_min, d1_sb_max, d2_sb_max = get_dwindow_values(data_spec, "charmless_2d")

    bm_max = 5600
    bm_min = 5000

    if "norm" in data_spec:
        bm_max = 5330
        bm_min = 5230

    rdf_base = rdf_base.Define("D1_M_Plot",f"D1_M - {d1_mstart}") \
                       .Define("D2_M_Plot",f"D2_M - {d2_mstart}") \
                       .Filter(f"B_M < {bm_max} && B_M > {bm_min}")

    d1xbins = 75
    d2ybins = 75

    xbinw = (d1_sb_max - d1_sb_min) / d1xbins
    ybinw = (d2_sb_max - d2_sb_min) / d2ybins

    hist_bm_sb = rdf_base.Histo1D((f"B Mass in D sideband region", f"B Mass in D sideband region", 75, bm_min, bm_max), 'B_M')
    hist_d1d2_sb = rdf_base.Histo2D((f"d1d2", f"d1d2", d1xbins, -d1_sb_max, d1_sb_max, d2ybins, -d2_sb_max, d2_sb_max), "D1_M_Plot", "D2_M_Plot")

    # ROOT.gStyle.SetOptStat("ne")

    # hist_string_d1 = f"{id_to_meson_string_dict[data_spec][1]} Mass in [{d1_mstart - 5*d}]"
    # hist_string_d2 = f"{id_to_meson_string_dict[data_spec][2]} Mass in "
    # entries = hist_bm_sb.Integral()
    #
    # dtpave = ROOT.TPaveText(0.60, 0.65, 0.80, 0.85, "NB NDC")
    # dtpave.SetFillStyle(0)
    # dtpave.SetTextSize(0.04)
    # dtpave.AddText(hist_string_d1)
    # dtpave.AddText(hist_string_d2)
    # dtpave.AddText(f"Number of Events: {entries}")


    # c1 = ROOT.TCanvas("c1","c1")
    # title = get_decay_string(id_to_meson_string_dict[data_spec], "Rec")
    # hist_bm_sb.GetXaxis().SetTitle(f"{title} [MeV]")
    # hist_bm_sb.Draw()
    # save_pdf(c1, f"charmless_plots", f"{data_spec}_sb_charmless")

    ROOT.gStyle.SetPalette(ROOT.kBird)

    xax = hist_d1d2_sb.GetXaxis()
    yax = hist_d1d2_sb.GetYaxis()

    xax.SetNdivisions(-10)
    yax.SetNdivisions(-10)

    for i in range(1, 12):
        xax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")
        yax.ChangeLabel(i,-1,0.05,-1,-1,-1,f"{i-6}#sigma")

    c2 = ROOT.TCanvas("c2","c2")
    c2.SetRightMargin(0.2)
    hist_d1d2_sb.GetXaxis().SetTitle(id_to_meson_string_dict[data_spec][1])
    hist_d1d2_sb.GetYaxis().SetTitle(id_to_meson_string_dict[data_spec][2])
    hist_d1d2_sb.GetZaxis().SetTitle(f"Candidates / {xbinw*ybinw:.2f} [MeV^{{2}}]")
    hist_d1d2_sb.GetYaxis().SetTitleOffset(1)
    hist_d1d2_sb.GetZaxis().SetTitleOffset(0.75)
    hist_d1d2_sb.Draw("COLZ")
    save_pdf(c2, f"charmless_plots", f"{data_spec}_d1d2_charmless")

# def build_tmc_mc(data_spec):
#     for mc_spec, list in id_to_b_values_dict.items():
#         if data_spec in mc_spec and "P_z_pst" not in mc_spec:
#             try:
#                 mcf = list[2]
#             except IndexError:
#                 mcf = "no"
#             if mcf == 0 or mcf == 1:
#
#                 mc_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/final_sample/{id_to_spec_dict[mc_spec]}.root"
#
#                 mc_tchain_sig = grab_files_and_chain(mc_path, ["DecayTreeTuple_SIG"])
#                 rdf_mc_sig = RDF(mc_tchain_sig)
#
#                 mc_tchain_sb = grab_files_and_chain(mc_path, ["DecayTreeTuple_SB"])
#                 rdf_mc_sb = RDF(mc_tchain_sb)
#
#                 truth_base_string = build_truth_base(mc_spec)
#
#                 print(rdf_mc_sig.Count().GetValue())
#
#                 rdf_mc_sig = rdf_mc_sig.Filter(f"({truth_base_string})", "truth_all")
#                 rdf_mc_sb = rdf_mc_sb.Filter(f"({truth_base_string})", "truth_all")
#
#                 tmc_sig_count = rdf_mc_sig.Count().GetValue()
#                 tmc_sb_count = rdf_mc_sb.Count().GetValue()
#
#                 print(tmc_sig_count, tmc_sb_count)
#
#                 outputfile = f"charm_mc/{mc_spec}.root"
#
#                 if os.path.exists(outputfile):
#                     os.remove(outputfile)
#                     print("First deleting old ", outputfile)
#                 else:
#                     print("making ", outputfile)
#
#                 rdf_mc_sig.Snapshot(f"DecayTreeTuple_SIG", outputfile, ["B_M"], opts)
#                 rdf_mc_sb.Snapshot(f"DecayTreeTuple_SB", outputfile, ["B_M"], opts)
#
#                 print(f"finished snapshot for {mc_spec}")

def mc_keys(mc_spec):

    sbws = ROOT.RooWorkspace("sbws")
    bmax = 5400
    bmin = 5000

    sbws.factory(f"B_M[{bmin},{bmax}]")
    b_m = sbws.var("B_M")
    data_args = ROOT.RooArgSet(b_m)

    mc_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/final_sample/{id_to_spec_dict[mc_spec]}.root"
    mc_tchain_sb = grab_files_and_chain(mc_path, ["DecayTreeTuple_SB"])
    rdf_mc_sb = RDF(mc_tchain_sb)

    data = ROOT.RooDataSet(f"{mc_spec}_data", f"{mc_spec}_data", mc_tchain_sb, data_args)

    sbws.Import(data)
    sbws.factory(f"RooKeysPdf::{mc_spec}_fit(B_M, {mc_spec}_data)")

    fit = sbws.pdf(f"{mc_spec}_fit")

    output_file = f"mc_keys/{mc_spec}.root"
    sbws.writeToFile(output_file)

    frame = b_m.frame(ROOT.RooFit.Title(f"{mc_spec}"), ROOT.RooFit.Bins(40))

    data.plotOn(frame, ROOT.RooFit.Name("data"))
    fit.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))

    p = residualPlot()
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

    title = get_decay_string(id_to_meson_string_dict[mc_spec], "Rec")
    pull.GetXaxis().SetTitle(f" m({title}) [MeV]")

    save_pdf(p, f"keys_pdfs", f"{mc_spec}_keys", rpflag = 1)

def get_ss_yield_and_fit(ws, data_spec, mc_spec, peak):

    mc_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/final_sample/{id_to_spec_dict[mc_spec]}.root"

    mc_tchain_sig = grab_files_and_chain(mc_path, ["DecayTreeTuple_SIG"])
    mc_tchain_sb = grab_files_and_chain(mc_path, ["DecayTreeTuple_SB"])

    sig_mc_events = mc_tchain_sig.GetEntries()
    sb_mc_events = mc_tchain_sb.GetEntries()

    sig_mc_events_ufloat = ufloat(sig_mc_events, np.sqrt(sig_mc_events))
    sb_mc_events_ufloat = ufloat(sb_mc_events, np.sqrt(sb_mc_events))

    data_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/data_files/{data_spec}_charm.root")
    data_ws = data_file.Get("fit_ws")

    data_yield = data_ws.var(f"{data_spec}_{peak}_yield")
    data_yield_err = data_yield.errorVar()
    data_yield_ufloat = ufloat(data_yield.getValV(), data_yield_err.getValV())

    ss_yield = (sb_mc_events_ufloat/sig_mc_events_ufloat)*data_yield_ufloat

    vars = data_ws.allVars()
    for i in vars:
        if i.GetName() != "B_M":
            temp = data_ws.var(i.GetName())
            temp.setConstant(True)
            print(f"{temp} is constant")

    shape_flag = ids_to_bestfit_dict[f"{data_spec}_{peak}"]

    fit = data_ws.pdf(f"{data_spec}_{peak}_fit_{shape_flag}")

    ws.Import(fit, ROOT.RooFit.RenameVariable(f"{data_spec}_{peak}_fit_{shape_flag}", f"{data_spec}_{peak}_fit"))
    # ws.factory(f"{data_spec}_{peak}_ss_yield[{ss_yield.n}]")

    ws.factory(f"{data_spec}_{peak}_ss_yield[{ss_yield.n}, {ss_yield.n - 5*ss_yield.s}, {ss_yield.n + 5*ss_yield.s}]")
    ws.factory(f"Gaussian::ss_{peak}_g({data_spec}_{peak}_ss_yield, {ss_yield.n}, {ss_yield.s})")

def fit_data_sb(data_spec, data_tchain_sb):

    ws = ROOT.RooWorkspace("fit_ws")
    bm_max = 5400
    bm_min = 5000
    ws.factory(f"B_M[{bm_min},{bm_max}]")

    if data_spec == "Z_m_p":

        get_ss_yield_and_fit(ws, data_spec, "01_Z_m_p", 0)
        get_ss_yield_and_fit(ws, data_spec, "02_Z_m_p", 1)

        ss_0_g_pdf = ws.pdf("ss_0_g")
        ss_1_g_pdf = ws.pdf("ss_1_g")

        glist = ROOT.RooArgSet(ss_0_g_pdf, ss_1_g_pdf)

        file_0 = ROOT.TFile(f"mc_keys/01_Z_m_p.root")
        file_1 = ROOT.TFile(f"mc_keys/02_Z_m_p.root")

        ws_0 = file_0.Get("sbws")
        ws_1 = file_1.Get("sbws")

        pdf_0 = ws_0.pdf(f"01_Z_m_p_fit")
        pdf_1 = ws_1.pdf(f"02_Z_m_p_fit")

        ws.Import(pdf_0)
        ws.Import(pdf_1)

        shape_flag_0 = ids_to_bestfit_dict[f"{data_spec}_0"]

        get_roofit_pdf(ws, f"{data_spec}_0", "G", "B_M")

        shape_flag_1 = ids_to_bestfit_dict[f"{data_spec}_1"]

        get_roofit_pdf(ws, f"{data_spec}_1", "G", "B_M")

        get_roofit_pdf(ws, f"{data_spec}", "BKG: Bernstein", "B_M")

        ws.factory(f"SUM::{data_spec}_sb_fit(\
                      Z_m_p_0_charmless_yield[55,0,1000]*Z_m_p_0_fit, \
                      Z_m_p_0_ss_yield*01_Z_m_p_fit, \
                      Z_m_p_1_charmless_yield[55,0,1000]*Z_m_p_1_fit, \
                      Z_m_p_1_ss_yield*02_Z_m_p_fit, \
                      Z_m_p_bkg_yield[100,0,10000]*Z_m_p_fit_bkg)")
    #
    ws.Print()
    b_m = ws.var("B_M")
    data_args = ROOT.RooArgSet(b_m)

    data = ROOT.RooDataSet(f"{data_spec}_sb", f"{data_spec}_sb", data_tchain_sb, data_args)
    model = ws.pdf(f"{data_spec}_sb_fit")

    fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(glist),ROOT.RooFit.Save())
    # fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0) ,ROOT.RooFit.Save())

    ws.Import(fit)
    fitresult = ws.obj(f"fitresult_{data_spec}_sb_fit_{data_spec}_sb")

    frame = b_m.frame(ROOT.RooFit.Title(f"{data_spec}_spectrum"), ROOT.RooFit.Bins(40))

    fit_list = ["Z_m_p_0_fit","01_Z_m_p_fit","Z_m_p_1_fit","02_Z_m_p_fit"]

    y_list = ["Z_m_p_0_charmless_yield","Z_m_p_0_ss_yield","Z_m_p_1_charmless_yield","Z_m_p_1_ss_yield"]

    data.plotOn(frame, ROOT.RooFit.Name("data"))
    model.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
    model.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_fit_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))

    color_list = [ROOT.kBlue, ROOT.kViolet, ROOT.kAzure, ROOT.kMagenta, ROOT.kCyan]

    for f, color in zip(fit_list, color_list):
        model.plotOn(frame, ROOT.RooFit.Components(f), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(color), ROOT.RooFit.Name(f"fit_{f}"))

    legend = ROOT.TLegend(id_to_legend_loc["rc"][0], id_to_legend_loc["rc"][1], id_to_legend_loc["rc"][2], id_to_legend_loc["rc"][3])
    legend.SetFillStyle(1001)

    d_nyield_bkg = ws.var(f"Z_m_p_bkg_yield")
    d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

    legend.AddEntry(frame.findObject("data"),"Run 2 Data","ep")
    legend.AddEntry(frame.findObject("pdf"),"Total Fit PDF","l")
    legend.AddEntry(frame.findObject("bkg"),f"Background PDF : {round(d_nyield_bkg.getValV(),1)} #pm {round(d_nyield_bkg_err,1)}","l")
    for fit_name, y_name in zip(fit_list, y_list):
        d_nyield = ws.obj(f"{y_name}")
        d_nyield_err = d_nyield.getPropagatedError(fitresult)
        legend.AddEntry(frame.findObject(f"fit_{fit_name}"), f"{y_name}: {round(d_nyield.getValV(),1)} #pm {round(d_nyield_err,1)}", "l")
    legend.SetTextSize(0.030)

    p = residualPlot()
    p.pt.cd()
    frame.SetMaximum(100)
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

    save_pdf(p, f"charm_fit", f"{data_spec}_sb", rpflag = 1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Charmless BKG EST")
    parser.add_argument('--Data_Spec_List', choices = ["Z_m_p","Z_z_z","P_z_p","M_m_z","norm7","norm8"],  nargs="+", help = 'Spec')
    parser.add_argument('--Plot_SB', action = 'store_true')
    parser.add_argument('--TMC_MC',action = 'store_true')
    parser.add_argument('--MC_KEYS',action = 'store_true')
    parser.add_argument('--FIT',action = 'store_true')

    args = parser.parse_args()

    data_spec_list = args.Data_Spec_List
    plot_sb = args.Plot_SB
    tmc_mc = args.TMC_MC
    mc_keys_flag = args.MC_KEYS
    fit = args.FIT

    if data_spec_list != None:
        for data_spec in data_spec_list:

            data_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{data_spec}.root"
            data_tchain_sig = grab_files_and_chain(data_path, ["DecayTreeTuple_SIG"])
            data_tchain_sb = grab_files_and_chain(data_path, ["DecayTreeTuple_SB"])

            if plot_sb:
                build_sb_plot(data_tchain_sb, data_spec)
                # get_ss_contribution(data_spec)
            if tmc_mc:
                build_tmc_mc(data_spec)

            if fit:
                fit_data_sb(data_spec, data_tchain_sb)

    if mc_keys_flag:
         for mc_spec, list in id_to_pid_dict.items():
             if data_spec_list != None and data_spec in mc_spec and "P_z_pst" not in mc_spec:
                 mc_keys(mc_spec)
