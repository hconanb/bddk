import ROOT as ROOT
from Analysis_2021.essentials import get_shapes_bkg, get_free_shapes, save_png, MakeSWeights
import glob as glob
from Analysis_2021 import rootutils
from uncertainties import ufloat
from uncertainties import covariance_matrix, correlation_matrix
from uncertainties import unumpy, ufloat_fromstr
import math as math
import pandas
import numpy as np
RDF = ROOT.ROOT.RDataFrame


def get_shape(dws, spec, run_name, varoi):

    data_ws_base = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/{spec}_{run_name}_{varoi}.root"")
    dws = data_ws_base.Get(run_name)

    # Get Signal Yields
    y_01 = dws.var("Z_m_p_01_yield")
    y_02 = dws.var("Z_m_p_0203_yield")

    fit_01 = dws.pdf("Z_m_p_01_fit")
    fit_02 = dws.pdf("Z_m_p_0203_fit")

    # Get TMC sb ratio









        if fix_flag:
            vars = mcws.allVars()
            for i in vars:
                if not fix_mean_flag:
                    if i.GetName() != "B_DTF_M" and i.GetName() != "B_M" and "mean" not in i.GetName():
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")
                else:
                    if i.GetName() != "B_DTF_M" and i.GetName() != "B_M":
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")
        print(smear_flag)
        if smear_flag:
            print("in smear",tuple)
            print(mc_spec)
            mc_pdf = mcws.pdf(f"{mc_spec}")
            # dws.factory(
            #     f"Gaussian::{mc_spec}_gsmear(B_DTF_M, mean_{base_spec}_gsmear, width_{base_spec}_gsmear)"
            # )
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable(mc_spec, f"{mc_spec}_pc"), ROOT.RooFit.RenameVariable("B_DTF_M", "B_M"))
            dws.factory(
                f"FCONV::{mc_spec}({varoi}, {mc_spec}_pc, {base_spec}_gsmear)"
            )
        if not smear_flag and mc_spec not in ["Z_z_z_0710_fit","Z_z_z_040812_fit","P_z_p_020607_fit","P_z_p_0408_fit"]:
            print("in not smear",tuple)
            print(mc_spec)
            mc_pdf = mcws.pdf(f"{mc_spec}")
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable("B_DTF_M", "B_M"))

def build_charmless_fit(run_name, spec):

    dws = ROOT.RooWorkspace(run_name)
    bmin = 5000
    bmax = 5600
    dws.factory(f"B_M[{bmin},{bmax}]")


    get_charmless_fit()


    dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_ch_yield[500,0,10000]*Z_m_p_01_charmless_fit, \
                                             Z_m_p_01_lmb_yield*Z_m_p_01_lmb_fit, \
                                             Z_m_p_0203_ch_yield[500,0,10000]*Z_m_p_0203_charmless_fit, \
                                             Z_m_p_0203_lmb_yield*Z_m_p_0203_lmb_fit, \
                                             Z_m_p_bkg_yield[100,0,100000]*Z_m_p_spectrum_bkg)")




    b_m = dws.var("B_M")
    data_args = ROOT.RooArgSet(b_m)
    #############################

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{spec}.root")
    tchain = ROOT.TChain("SB")
    for file_name in file_list:
        tchain.Add(file_name)

    data = ROOT.RooDataSet(f"{spec}_final_data", f"{spec}_final_data", tchain, data_args)
    model = dws.pdf(f"{spec}_spectrum_all_fit")
    fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    frame = b_m.frame(ROOT.RooFit.Title(f"{data_spec}_spectrum"))

    if spec == "Z_m_p":
        ylist = ["01","0203","04"]
        title = "D^{-} D^{+} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{+}"

    data.plotOn(frame, ROOT.RooFit.Name("data"))
    model.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
    fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))

    color_list = [ROOT.kBlue, ROOT.kViolet, ROOT.kAzure, ROOT.kMagenta, ROOT.kCyan]

    for y, color in zip(ylist, color_list):
        fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_{y}_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(color), ROOT.RooFit.Name(f"fit_{y}"))


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
    pull = b_dtf_m.frame()
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

    pull.GetXaxis().SetTitle(f" m({title} {varoi}) [MeV]")
    # frame.GetXaxis().SetTitle(f" m({title}) [MeV]")

    save_png(p, f"fit_tests", f"{spec}_{run_name}", rpflag = 1)

    # dws.Import(data)
    # dws.Import(fit)
    #
    # output_base_data = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/{spec}_{run_name}.root"
    # dws.writeToFile(output_base_data)
    # print(f"Wrote dws to: {output_base_data}")
    # dws.Print()
