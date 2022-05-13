import sys
import os
import ROOT
import uproot
import pandas as pd
import argparse

basedir = os.getcwd().split('fits')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

def fit_data(version, data_spec, peak, varoi, tchain):

    if version == "sWeights":

        b_mean = id_to_b_values_dict[f"{data_spec}_{peak}"][0]
        b_window = id_to_b_values_dict[f"{data_spec}_{peak}"][1]

        bmin = b_mean - b_window
        bmax = b_mean + b_window

        dws = ROOT.RooWorkspace(f"fit_ws")
        dws.factory(f"{varoi}[{bmin},{bmax}]")

        shape_flag = ids_to_bestfit_dict[f"{data_spec}_{peak}"]

        get_roofit_pdf(dws, f"{data_spec}_{peak}", shape_flag, varoi)
        get_roofit_pdf(dws, f"{data_spec}", "BKG: Exponential", varoi)

        dws.Print()

        dws.factory(f"SUM::{data_spec}_{peak}_all_fit({data_spec}_{peak}_yield[500,0,10000]*{data_spec}_{peak}_fit_{shape_flag}, {data_spec}_bkg_yield[100,0,100000]*{data_spec}_fit_bkg)")

        fit_var = dws.var(varoi)
        data_args = ROOT.RooArgSet(fit_var)

        data_set = ROOT.RooDataSet(f"{data_spec}_{peak}_events", f"{data_spec}_{peak}_events", tchain, data_args)
        fit_pdf = dws.pdf(f"{data_spec}_{peak}_all_fit")
        fit_result = fit_pdf.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

        fit_ws = ROOT.RooWorkspace(f"fit_ws")

        fit_ws.Import(fit_pdf)
        fit_ws.Import(fit_result)
        fit_ws.Import(data_set)

        fit_output_file = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/data_files/{data_spec}_{peak}.root"
        fit_ws.writeToFile(fit_output_file)

        nyield_bkg = dws.var(f"{data_spec}_bkg_yield")
        nyield_1 = dws.var(f"{data_spec}_{peak}_yield")
        yields = ROOT.RooArgSet(nyield_1, nyield_bkg)
        MakeSWeights(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/sw_files/sw_{data_spec}_{peak}.root", "SW_tree", data_set, fit_pdf, yields)
    if version == "charmless":

        dws = ROOT.RooWorkspace("ws")

        bmin = 4800
        bmax = 5600

        if "norm" in data_spec:
            bmin = 5230
            bmax = 5330

        dws.factory(f"{varoi}[{bmin},{bmax}]")

        if data_spec not in ["M_m_z"]:
            shape_flag_0 = ids_to_bestfit_dict[f"{data_spec}_0"]
            get_roofit_pdf(dws, f"{data_spec}_0", shape_flag_0, varoi)
        if data_spec not in ["norm7","norm8"]:
            shape_flag_1 = ids_to_bestfit_dict[f"{data_spec}_1"]
            get_roofit_pdf(dws, f"{data_spec}_1", shape_flag_1, varoi)
        if data_spec not in ["P_z_pst","norm7","norm8"]:
            shape_flag_2 = ids_to_bestfit_dict[f"{data_spec}_2"]
            get_roofit_pdf(dws, f"{data_spec}_2", shape_flag_2, varoi)


        if data_spec in ["Z_m_p","Z_z_z","P_z_p"]:
            get_roofit_pdf(dws, data_spec, "BKG: Bernstein", varoi)
            dws.factory(f"SUM::{data_spec}_charm_all_fit({data_spec}_0_yield[500,0,10000]*{data_spec}_0_fit_{shape_flag_0}, \
                                                         {data_spec}_1_yield[500,0,10000]*{data_spec}_1_fit_{shape_flag_1}, \
                                                         {data_spec}_2_yield[500,0,10000]*{data_spec}_2_fit_{shape_flag_2}, \
                                                         {data_spec}_bkg_yield[100,0,100000]*{data_spec}_fit_bkg)")
        if data_spec in ["M_m_z"]:
            get_roofit_pdf(dws, data_spec, "BKG: Bernstein", varoi)
            dws.factory(f"SUM::{data_spec}_charm_all_fit({data_spec}_1_yield[500,0,10000]*{data_spec}_1_fit_{shape_flag_1}, \
                                                         {data_spec}_2_yield[500,0,10000]*{data_spec}_2_fit_{shape_flag_2}, \
                                                         {data_spec}_bkg_yield[100,0,100000]*{data_spec}_fit_bkg)")
        if data_spec in ["P_z_pst"]:
            get_roofit_pdf(dws, data_spec, "BKG: Bernstein", varoi)
            dws.factory(f"SUM::{data_spec}_charm_all_fit({data_spec}_0_yield[500,0,10000]*{data_spec}_0_fit_{shape_flag_0}, \
                                                         {data_spec}_1_yield[500,0,10000]*{data_spec}_1_fit_{shape_flag_1}, \
                                                         {data_spec}_bkg_yield[100,0,100000]*{data_spec}_fit_bkg)")
        if data_spec in ["norm7","norm8"]:
            get_roofit_pdf(dws, data_spec, "BKG: Exponential", varoi)
            dws.factory(f"SUM::{data_spec}_charm_all_fit({data_spec}_0_yield[500,0,10000]*{data_spec}_0_fit_{shape_flag_0}, \
                                                         {data_spec}_bkg_yield[100,0,100000]*{data_spec}_fit_bkg)")
        fit_var = dws.var(varoi)
        data_args = ROOT.RooArgSet(fit_var)

        data_set = ROOT.RooDataSet(f"{data_spec}_charm_events", f"{data_spec}_charm_events", tchain, data_args)
        fit_pdf = dws.pdf(f"{data_spec}_charm_all_fit")
        fit_result = fit_pdf.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

        fit_ws = ROOT.RooWorkspace(f"fit_ws")

        fit_ws.Import(fit_pdf)
        fit_ws.Import(fit_result)
        fit_ws.Import(data_set)

        fit_output_file = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/data_files/{data_spec}_charm.root"
        fit_ws.writeToFile(fit_output_file)

def plot_data(version, data_spec, peak, varoi):

    if version == "charmless":
        peak = "charm"

    base_fit_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/data_files/{data_spec}_{peak}.root")

    fit_ws = base_fit_file.Get(f"fit_ws")
    fit_var = fit_ws.var(f"{varoi}")

    fit_pdf = fit_ws.pdf(f"{data_spec}_{peak}_all_fit")
    data_set = fit_ws.data(f"{data_spec}_{peak}_events")

    fitresult = fit_ws.obj(f"fitresult_{data_spec}_{peak}_all_fit_{data_spec}_{peak}_events")
    frame = fit_var.frame(ROOT.RooFit.Title(f"{data_spec}_{peak}_spectrum"))

    data_set.plotOn(frame, ROOT.RooFit.Name("data"))
    fit_pdf.plotOn(frame, ROOT.RooFit.Name("pdf"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{data_spec}_pdf"])))
    fit_pdf.plotOn(frame, ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{data_spec}_bkg"])), ROOT.RooFit.Components(f"{data_spec}_fit_bkg"))

    if version == "sWeights":

        shape_flag = ids_to_bestfit_dict[f"{data_spec}_{peak}"]

        if shape_flag == "DG" or "_fr" in shape_flag:
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_{peak}_fit_{shape_flag}_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fa"))
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_{peak}_fit_{shape_flag}_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fb"))

    if version == "charmless":
        if data_spec in ["Z_m_p","Z_z_z","P_z_p"]:
            shape_flag_0 = ids_to_bestfit_dict[f"{data_spec}_0"]
            shape_flag_1 = ids_to_bestfit_dict[f"{data_spec}_1"]
            shape_flag_2 = ids_to_bestfit_dict[f"{data_spec}_2"]
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_0_fit_{shape_flag_0}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"p0_{shape_flag_0}"))
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_1_fit_{shape_flag_1}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"p1_{shape_flag_1}"))
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_2_fit_{shape_flag_2}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"p2_{shape_flag_2}"))

        if data_spec in ["M_m_z"]:
            shape_flag_1 = ids_to_bestfit_dict[f"{data_spec}_1"]
            shape_flag_2 = ids_to_bestfit_dict[f"{data_spec}_2"]
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_1_fit_{shape_flag_1}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"p1_{shape_flag_1}"))
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_2_fit_{shape_flag_2}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"p2_{shape_flag_2}"))

        if data_spec in ["P_z_pst"]:
            shape_flag_0 = ids_to_bestfit_dict[f"{data_spec}_0"]
            shape_flag_1 = ids_to_bestfit_dict[f"{data_spec}_1"]
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_0_fit_{shape_flag_0}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"p0_{shape_flag_0}"))
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_1_fit_{shape_flag_1}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"p1_{shape_flag_1}"))

        if data_spec in ["norm7","norm8"]:
            shape_flag_0 = ids_to_bestfit_dict[f"{data_spec}_0"]
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{data_spec}_0_fit_{shape_flag_0}"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"p0_{shape_flag_0}"))



    if data_spec == "Z_z_z":
        legend = ROOT.TLegend(id_to_legend_loc["rc"][0], id_to_legend_loc["rc"][1], id_to_legend_loc["rc"][2], id_to_legend_loc["rc"][3])
    else:
        legend = ROOT.TLegend(id_to_legend_loc["rc2"][0], id_to_legend_loc["rc2"][1], id_to_legend_loc["rc2"][2], id_to_legend_loc["rc2"][3])

    d_chi2 = frame.chiSquare(f"pdf", f"data")

    legend.SetFillStyle(1001)


    legend.AddEntry(frame.findObject("data"),"Run 2 Data","ep")
    legend.AddEntry(frame.findObject("pdf"),"Total Fit PDF","l")

    if version == "sWeights":
        if shape_flag == "DG":
            legend.AddEntry(frame.findObject("fa"),"Gauss 1","l")
            legend.AddEntry(frame.findObject("fb"),"Gauss 2","l")
        if "_fr" in shape_flag:
            shape_1 = shape_flag.split("Add")[0]
            shape_2 = shape_flag.split("Add")[1].split("_fr")[0]
            legend.AddEntry("fa",f"{shape_1}","l")
            legend.AddEntry("fb",f"{shape_2}","l")

    if data_spec != "M_m_z" and data_spec != "P_z_pst" and "norm" not in data_spec:
        y_list = [f"{data_spec}_0_yield", f"{data_spec}_1_yield", f"{data_spec}_2_yield", f"{data_spec}_bkg_yield"]
        fit_list = [f"{data_spec}_0_fit_{shape_flag_0}", f"{data_spec}_1_fit_{shape_flag_1}",  f"{data_spec}_2_fit_{shape_flag_2}", f"{data_spec}_fit_bkg"]
        name_list = [f"p0_{shape_flag_0}", f"p1_{shape_flag_1}", f"p2_{shape_flag_2}", "bkg"]
    if data_spec == "M_m_z":
        y_list = [f"{data_spec}_1_yield", f"{data_spec}_2_yield", f"{data_spec}_bkg_yield"]
        fit_list = [f"{data_spec}_1_fit_{shape_flag_1}",  f"{data_spec}_2_fit_{shape_flag_2}", f"{data_spec}_fit_bkg"]
        name_list = [f"p1_{shape_flag_1}", f"p2_{shape_flag_2}", "bkg"]
    if data_spec == "P_z_pst":
        y_list = [f"{data_spec}_0_yield", f"{data_spec}_1_yield", f"{data_spec}_bkg_yield"]
        fit_list = [f"{data_spec}_0_fit_{shape_flag_0}", f"{data_spec}_1_fit_{shape_flag_1}", f"{data_spec}_fit_bkg"]
        name_list = [f"p0_{shape_flag_0}", f"p1_{shape_flag_1}", "bkg"]
    if data_spec == "norm7" or data_spec == "norm8":
        y_list = [f"{data_spec}_0_yield", f"{data_spec}_bkg_yield"]
        fit_list = [f"{data_spec}_0_fit_{shape_flag_0}", f"{data_spec}_fit_bkg"]
        name_list = [f"p0_{shape_flag_0}", "bkg"]

    if version == "charmless":
        for y, fit, name in zip(y_list, fit_list, name_list):
            d_nyield = fit_ws.obj(y)
            d_nyield_err = d_nyield.getPropagatedError(fitresult)
            print(d_nyield.getValV())
            legend.AddEntry(frame.findObject(name), f"{name}: {round(d_nyield.getValV(),1)} #pm {round(d_nyield_err,1)}","l")


    legend.SetTextSize(0.035)

    p = residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)

    frame.Draw()
    legend.Draw()
    p.pb.cd()

    hpull = frame.pullHist(f"data", f"pdf")
    pull = fit_var.frame()
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

    title = get_decay_string(id_to_meson_string_dict[data_spec], "Rec")
    pull.GetXaxis().SetTitle(f" m({title}) [MeV]")

    save_pdf(p, f"data_fits", f"{data_spec}_{peak}_{version}", rpflag = 1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Fit a PDF to separate data samples")
    parser.add_argument('--Version', choices=["sWeights","charmless"], help = 'Version of Fit')
    parser.add_argument('--Data_Spec', choices=["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","norm7","norm8"], help = 'Data Spec')
    parser.add_argument('--Varoi', choices=["B_DTF_M","B_M"], default = "B_DTF_M", help = 'Variable to Fit')
    parser.add_argument('--FIT', action='store_true')
    parser.add_argument('--PLOT', action='store_true')
    parser.add_argument('--PEAK', type = int, choices=[0,1,2])

    args = parser.parse_args()

    version = args.Version
    data_spec = args.Data_Spec
    varoi = args.Varoi
    fit_flag = args.FIT
    plot_flag = args.PLOT
    peak_flag = args.PEAK


    file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{id_to_spec_dict[data_spec]}.root"
    DecayTree_List = ["DecayTreeTuple_SIG"]
    tchain = grab_files_and_chain(file_path, DecayTree_List)

    if fit_flag:
        fit_data(version, id_to_spec_dict[data_spec], peak_flag, varoi, tchain)
    if plot_flag:
        plot_data(version, id_to_spec_dict[data_spec], peak_flag, varoi)
