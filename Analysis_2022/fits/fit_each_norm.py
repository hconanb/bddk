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

def MakeSWeights(outfilename, outtreename, data, model, yields):

    """Determine s-weights from fit.

    arguments:
    outfilename -- name of .root file to create with `outtreename`
    outtreename -- name of TTree with s-weights to save in `outfilename`
    data -- RooDataSet to which `model` was fitted
    model -- fitted RooAbsPdf
    yields -- RooArgList of RooRealVars extracted from fitting `model` to `data`
    """
    from array import array

    print(f"using data '{data.GetName()}'")
    print(f"using model '{model.GetName()}'")
    print(f"using yields '{[x.GetName() for x in yields]}'")

    # ROOT.RooMsgService.instance().Print()
    # rme = ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.INFO)
    # ROOT.RooMsgService.instance().addStream(ROOT.RooFit.DEBUG, ROOT.RooFit.Topic(ROOT.RooFit.Contents), ROOT.RooFit.OutputFile("rfsmear_debug.log"))
    #projDeps = ROOT.RooArgSet(), useWeights = False, copyDataSet = True, newName = "test", fitToarg5 =  ROOT.RooFit.PrintLevel(3)

    sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields)
    print("npds")

    # print("Check SWeights:")
    # for y in yields:
    #     oval = y.getVal()
    #     sval = sData.GetYieldFromSWeight(y.GetName())
    #     print(f"Yield of {y.GetName()} is {oval}")
    #     print(f"from sWeights it is {sval}")
    #     if not (0.9995 < oval / sval < 1.0005):
    #         raise Exception("sWeighted yield should match")
    # for i in range(10):
    #     for y in yields:
    #         print(f"    {y.GetName()} Weight {sData.GetSWeight(i, y.GetName())}")
    #     totw = sData.GetSumOfEventSWeight(i)
    #     print(f"Total Weight {totw}")
    #     if not (0.9995 < totw < 1.0005):
    #         raise Exception("sum of sWeight should be 1")
    swnames = sorted([f"{x.GetName()}_sw" for x in yields])
    print(f"weights: {swnames}")
    # create output file
    nf = ROOT.TFile.Open(outfilename, "recreate")
    # create directory hierarchy
    nd = nf
    if len(outtreename.split("/")) > 1:
        for d in outtreename.split("/")[:-1]:
            nd = nd.mkdir(d)
    nd.cd()
    # create output TTree
    nt = ROOT.TTree(outtreename.split("/")[-1], outtreename.split("/")[-1])
    # declare output branches
    swvals = [array("f", [0]) for x in swnames]
    for val, nm in zip(swvals, swnames):
        nt.Branch(nm, val, f"{nm}/F")
    # loop data
    for i in range(data.numEntries()):
        # get vars
        swvars = sorted(
            [x for x in data.get(i) if x.GetName() in swnames],
            key=lambda x: x.GetName(),
        )
        assert [x.GetName() for x in swvars] == swnames  # check sorting worked
        # set values
        for val, var in zip(swvals, swvars):
            val[0] = var.getVal()
        # fill values
        nt.Fill()
    nt.Write()
    nf.Close()

def fit_norm(norm_spec, varoi, tchain, break_flag = False, year = "ALL", trigger = "ALL"):

    nws = ROOT.RooWorkspace("fit_ws")

    bmin = 5230
    bmax = 5330

    nws.factory(f"{varoi}[{bmin},{bmax}]")

    if not break_flag:
        nws.factory(f"Gaussian::{norm_spec}_signal_a({varoi}, mean_{norm_spec}[5279, 5270, 5290], width_{norm_spec}_a[10, 0.1, 15.0])")
        nws.factory(f"Gaussian::{norm_spec}_signal_b({varoi}, mean_{norm_spec}, width_{norm_spec}_b[5, 0.1, 10.0])")
        nws.factory(f"SUM::{norm_spec}_signal({norm_spec}_a_frac[0.5,0,1]*{norm_spec}_signal_a, {norm_spec}_signal_b)")
    if break_flag:
        nws.factory(f"Gaussian::{norm_spec}_signal({varoi}, mean_{norm_spec}[5279, 5270, 5290], width_{norm_spec}_a[10, 0.1, 15.0])")

    get_roofit_pdf(nws, norm_spec, "BKG: Exponential", varoi)
    nws.factory(f"SUM::{norm_spec}_spectrum_all_fit(n_{norm_spec}_signal[100,0,100000]*{norm_spec}_signal,n_{norm_spec}_bkg[100,0,1000000]*{norm_spec}_fit_bkg)")

    fit_var = nws.var(varoi)
    data_args = ROOT.RooArgSet(fit_var)

    data_set = ROOT.RooDataSet(f"{norm_spec}_events", f"{norm_spec}_events", tchain, data_args)
    fit_pdf = nws.pdf(f"{norm_spec}_spectrum_all_fit")
    fit_result = fit_pdf.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    fit_ws = ROOT.RooWorkspace(f"fit_ws")

    fit_ws.Import(fit_pdf)
    fit_ws.Import(fit_result)
    fit_ws.Import(data_set)

    if not break_flag:
        nyield_1 = nws.var(f"n_{norm_spec}_signal")
        nyield_2 = nws.var(f"n_{norm_spec}_bkg")
        yields = ROOT.RooArgSet(nyield_1, nyield_2)
        MakeSWeights(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/sw_files/sw_{norm_spec}.root", "SW_tree", data_set, fit_pdf, yields)

    fit_output_file = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/data_files/{norm_spec}_{year}_{trigger}.root"
    fit_ws.writeToFile(fit_output_file)


def plot_norm(norm_spec, varoi, break_flag, year, trigger):

    base_fit_file_path = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/data_files/{norm_spec}_{year}_{trigger}.root"

    base_fit_file = ROOT.TFile(base_fit_file_path)

    fit_ws = base_fit_file.Get(f"fit_ws")

    fit_var = fit_ws.var(f"{varoi}")
    fit_pdf = fit_ws.pdf(f"{norm_spec}_spectrum_all_fit")
    data_set = fit_ws.data(f"{norm_spec}_events")

    fitresult = fit_ws.obj(f"fitresult_{norm_spec}_spectrum_all_fit_{norm_spec}_events")
    frame = fit_var.frame(ROOT.RooFit.Title(f"{norm_spec}_spectrum"))

    data_set.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.Binning(40))

    fit_pdf.plotOn(frame, ROOT.RooFit.Name("pdf"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{norm_spec}_pdf"])))
    fit_pdf.plotOn(frame, ROOT.RooFit.Name("sig"), ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{norm_spec}_signal"])), ROOT.RooFit.Components(f"{norm_spec}_signal"))
    if not break_flag:
        fit_pdf.plotOn(frame, ROOT.RooFit.Name("sig_a"), ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{norm_spec}_signal_a"])), ROOT.RooFit.Components(f"{norm_spec}_signal_a"))
        fit_pdf.plotOn(frame, ROOT.RooFit.Name("sig_b"), ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{norm_spec}_signal_b"])), ROOT.RooFit.Components(f"{norm_spec}_signal_b"))
    fit_pdf.plotOn(frame, ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid),ROOT.RooFit.LineColor(ROOT.TColor.GetColor(id_to_color_values_dict[f"{norm_spec}_bkg"])), ROOT.RooFit.Components(f"{norm_spec}_fit_bkg"))

    legend = ROOT.TLegend(id_to_legend_loc[norm_spec][0],id_to_legend_loc[norm_spec][1],id_to_legend_loc[norm_spec][2],id_to_legend_loc[norm_spec][3])

    d_chi2 = frame.chiSquare(f"pdf", f"data")

    legend.SetFillStyle(1001)

    d_nyield_bkg = fit_ws.obj(f"n_{norm_spec}_bkg")
    d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

    d_nyield_sig = fit_ws.obj(f"n_{norm_spec}_signal")
    d_nyield_sig_err = d_nyield_sig.getPropagatedError(fitresult)

    legend.AddEntry(frame.findObject("data"),"Run 2 Data","ep")
    legend.AddEntry(frame.findObject("pdf"),"Total Fit PDF","l")
    legend.AddEntry(frame.findObject("sig"),f"Signal PDF : {round(d_nyield_sig.getValV(),1)} #pm {round(d_nyield_sig_err,1)}","l")
    if not break_flag:
        legend.AddEntry(frame.findObject("sig_a"),f"Signal PDF A","l")
        legend.AddEntry(frame.findObject("sig_b"),f"Signal PDF B","l")
    legend.AddEntry(frame.findObject("bkg"),f"Background PDF : {round(d_nyield_bkg.getValV(),1)} #pm {round(d_nyield_bkg_err,1)}","l")
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

    title = get_decay_string(id_to_meson_string_dict[norm_spec], "Rec")
    pull.GetXaxis().SetTitle(f" m({title}) [MeV]")

    save_pdf(p, f"norm_fits", f"{norm_spec}_{varoi}_{year}_{trigger}", rpflag = 1)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Fit a PDF to separate data samples")
    parser.add_argument('--Norm_Spec', choices=["norm7","norm8"], help = 'Norm Spec')
    parser.add_argument('--Varoi', choices=["B_DTF_M","B_M"], default = "B_DTF_M", help = 'Variable to Fit')
    parser.add_argument('--BREAK', action='store_true')
    parser.add_argument('--FIT', action='store_true')
    parser.add_argument('--PLOT', action='store_true')

    args = parser.parse_args()

    norm_spec = args.Norm_Spec
    varoi = args.Varoi
    break_flag = args.BREAK
    fit_flag = args.FIT
    plot_flag = args.PLOT

    if break_flag:
        for year in ["2016","2017","2018"]:
            for trigger in ["T","nTaT"]:
                file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/final_sample/{norm_spec}.root"
                DecayTree_List = [f"DecayTreeTuple_{trigger}"]
                tchain = grab_files_and_chain(file_path, DecayTree_List)
                if fit_flag:
                    fit_norm(norm_spec, varoi, tchain, break_flag, year, trigger)
                if plot_flag:
                    plot_norm(norm_spec, varoi, break_flag, year, trigger)
    else:
        file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{norm_spec}.root"
        DecayTree_List = [f"DecayTreeTuple_SIG"]
        tchain = grab_files_and_chain(file_path, DecayTree_List)
        if fit_flag:
            fit_norm(norm_spec, varoi, tchain)
        if plot_flag:
            plot_norm(norm_spec, varoi, year, trigger)
