import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *


PE_d1d2 = "(pow(D1_PE + D2_PE,2))"
PX_d1d2 = "(pow(D1_PX + D2_PX,2))"
PY_d1d2 = "(pow(D1_PY + D2_PY,2))"
PZ_d1d2 = "(pow(D1_PZ + D2_PZ,2))"
M2_d1d2 = f"{PE_d1d2} - ({PX_d1d2} + {PY_d1d2} + {PZ_d1d2})"
M_d1d2 = f"sqrt({M2_d1d2})"

PE_d1kst = "(pow(D1_PE + KST_PE,2))"
PX_d1kst = "(pow(D1_PX + KST_PX,2))"
PY_d1kst = "(pow(D1_PY + KST_PY,2))"
PZ_d1kst = "(pow(D1_PZ + KST_PZ,2))"
M2_d1kst = f"{PE_d1kst} - ({PX_d1kst} + {PY_d1kst} + {PZ_d1kst})"
M_d1kst = f"sqrt({M2_d1kst})"

PE_d2kst = "(pow(D2_PE + KST_PE,2))"
PX_d2kst = "(pow(D2_PX + KST_PX,2))"
PY_d2kst = "(pow(D2_PY + KST_PY,2))"
PZ_d2kst = "(pow(D2_PZ + KST_PZ,2))"
M2_d2kst = f"{PE_d2kst} - ({PX_d2kst} + {PY_d2kst} + {PZ_d2kst})"
M_d2kst = f"sqrt({M2_d2kst})"
#########################################MCDecayTreeTuple
PE_d1d2_mcdtt = "(pow(C1_TRUEP_E + C2_TRUEP_E,2))"
PX_d1d2_mcdtt = "(pow(C1_TRUEP_X + C2_TRUEP_X,2))"
PY_d1d2_mcdtt = "(pow(C1_TRUEP_Y + C2_TRUEP_Y,2))"
PZ_d1d2_mcdtt = "(pow(C1_TRUEP_Z + C2_TRUEP_Z,2))"
M2_d1d2_mcdtt = f"{PE_d1d2_mcdtt} - ({PX_d1d2_mcdtt} + {PY_d1d2_mcdtt} + {PZ_d1d2_mcdtt})"
M_d1d2_mcdtt = f"sqrt({M2_d1d2_mcdtt})"

PE_d1kst_mcdtt = "(pow(C1_TRUEP_E + Kst_TRUEP_E,2))"
PX_d1kst_mcdtt = "(pow(C1_TRUEP_X + Kst_TRUEP_X,2))"
PY_d1kst_mcdtt = "(pow(C1_TRUEP_Y + Kst_TRUEP_Y,2))"
PZ_d1kst_mcdtt = "(pow(C1_TRUEP_Z + Kst_TRUEP_Z,2))"
M2_d1kst_mcdtt = f"{PE_d1kst_mcdtt} - ({PX_d1kst_mcdtt} + {PY_d1kst_mcdtt} + {PZ_d1kst_mcdtt})"
M_d1kst_mcdtt = f"sqrt({M2_d1kst_mcdtt})"

PE_d2kst_mcdtt = "(pow(C2_TRUEP_E + Kst_TRUEP_E,2))"
PX_d2kst_mcdtt = "(pow(C2_TRUEP_X + Kst_TRUEP_X,2))"
PY_d2kst_mcdtt = "(pow(C2_TRUEP_Y + Kst_TRUEP_Y,2))"
PZ_d2kst_mcdtt = "(pow(C2_TRUEP_Z + Kst_TRUEP_Z,2))"
M2_d2kst_mcdtt = f"{PE_d2kst_mcdtt} - ({PX_d2kst_mcdtt} + {PY_d2kst_mcdtt} + {PZ_d2kst_mcdtt})"
M_d2kst_mcdtt = f"sqrt({M2_d2kst_mcdtt})"


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

    sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields)

    print("Check SWeights:")
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

    tltxt = ROOT.TLatex()

def plotdalitz(k, title, xname, yname, histo, dim, rflag = 0):

    tltxt = ROOT.TLatex()

    if dim == 2 and rflag == 0:
        histo.GetXaxis().SetTitle(f"m^{{2}}({xname}) [GeV^{{2}}]")
        histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
        histo.GetYaxis().SetTitle(f"m^{{2}}({yname}) [GeV^{{2}}]")
        histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
    if dim == 1 and rflag == 0:
        histo.GetXaxis().SetTitle(f"m({xname}) [GeV]")
        histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000,histo.GetXaxis().GetXmax()/1000)
        bw = histo.GetXaxis().GetBinWidth(5)
        histo.GetYaxis().SetTitle(f"{yname} / ({bw} GeV)")
        histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000,histo.GetYaxis().GetXmax()/1000)
    if dim == 2 and rflag == 1:
        histo.GetXaxis().SetTitle(f"ratio")
        # histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
        histo.GetYaxis().SetTitle(f"ratio")
        # histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
    if dim == 1 and rflag == 1:
        histo.GetXaxis().SetTitle(f"ratio")
        # histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000,histo.GetXaxis().GetXmax()/1000)
        bw = histo.GetXaxis().GetBinWidth(5)
        histo.GetYaxis().SetTitle(f"{yname} / ({bw} GeV)")
        # histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000,histo.GetYaxis().GetXmax()/1000)

    c1 = ROOT.TCanvas("c1","c1")
    c1.SetRightMargin(0.15)
    c1.SetLeftMargin(0.9)
    c1.SetTopMargin(0.1)
    c1.SetBottomMargin(0.9)
    ht1 = histo.DrawCopy("COLZ")

    tltxt.DrawLatexNDC(.4,.925, title)

    pname = histo.GetTitle()
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}_pdf/'):
        os.makedirs(f'plots/{now.month}_{now.day}_pdf/')
    if not os.path.exists(f'plots/{now.month}_{now.day}_png/'):
        os.makedirs(f'plots/{now.month}_{now.day}_png/')
    c1.SaveAs(f"plots/{now.month}_{now.day}_pdf/{k}_{pname}.pdf")
    c1.SaveAs(f"plots/{now.month}_{now.day}_png/{k}_{pname}.png")

def prepdalitz(rdf, spec):

    dbins = 10

    if "sw" in spec:
        hdalitz_12_1k= rdf.Histo2D(("", f"d1d2_d1k", dbins, 13e6, 24e6, dbins, 5e6, 13e6), "M2_d1d2", "M2_d1k", "sig_sw")
        hdalitz_12_1k= rdf.Histo2D(("", f"d1d2_d1k", dbins, 13e6, 24e6, dbins, 5e6, 13e6), "M2_d1d2", "M2_d1k")
        hdalitz_1k_p = rdf.Histo1D(("", f"d1k_p", dbins, 1500, 4000), "M_d1k", "sig_sw")
        hdalitz_12_p = rdf.Histo1D(("", f"d1d2_p", dbins, 3500, 5500), "M_d1d2", "sig_sw")

    else:
        # hdalitz_12_2k = rdf.Histo2D(("", f"d1d2_d2kst", dbins, 14.25e6, 21e6, dbins, 7e6, 13e6), "M2_d1d2", "M2_d2k")
        hdalitz_12_1k = rdf.Histo2D(("", f"d1d2_d1k", dbins, 13e6, 24e6, dbins, 5e6, 13e6), "M2_d1d2", "M2_d1k")
        # hdalitz_1k_2k = rdf.Histo2D(("", f"d1kst_d2kst", dbins, 7e6, 13e6, dbins, 7e6, 13e6), "M2_d1kst", "M2_d2kst")
        hdalitz_1k_p = rdf.Histo1D(("", f"d1k_p", dbins, 1500, 4000), "M_d1k")
        # # hdalitz_2k_p = rdf.Histo1D(("", f"d2kst_p", dbins, 2500, 4000), "M_d2kst")
        hdalitz_12_p = rdf.Histo1D(("", f"d1d2_p", dbins, 3500, 5500), "M_d1d2")

    title = "B^{0} #rightarrow D^{-} D^{0} K^{+}"
    d1 = "D^{-}"
    d2 = "(D^{0} #rightarrow k#pi#pi#pi)"
    d3 = "K^{+}"

    plotdalitz(spec, title, f"{d1} {d2}", f"{d1} {d3}", hdalitz_12_1k, 2)
    plotdalitz(spec, title, f"{d1} {d3}", "Candidates", hdalitz_1k_p, 1)
    plotdalitz(spec, title, f"{d1} {d2}", "Candidates", hdalitz_12_p, 1)

# def Dalitz_Pre_SWeight(spec, year, trigger, data_tree, mc_rec_tree, mc_gen_tree):
def Dalitz_Pre_SWeight(spec, year, trigger, data_tree, mc_rec_tree, mc_gen_hist, mygentree):

    sw_file = ROOT.TFile(f"base_norm_files/{spec}_sw_file.root")
    sw_tree = sw_file.Get("norm8_sw")
    data_tree.AddFriend(sw_tree,  "sw_tree")

    rdf_data = RDF(data_tree)
    rdf_mc_rec = RDF(mc_rec_tree)
    rdf_mc_gen = RDF(mygentree)

    rdf_data_n1 = rdf_data.Define("sig_sw", f"sw_tree.n_{spec}_signal_sw") \
                         .Define("M2_d1d2", M2_d1d2.replace("KST","K")) \
                         .Define("M2_d1k", M2_d1kst.replace("KST","K")) \
                         .Define("M_d1d2",M_d1d2) \
                         .Define("M_d1k",M_d1kst.replace("KST","K")) \
                         .Filter("(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)")

    rdf_mc_rec_n1 = rdf_mc_rec.Define("M2_d1d2", M2_d1d2.replace("KST","K")) \
                              .Define("M2_d1k", M2_d1kst.replace("KST","K")) \
                              .Define("M_d1d2",M_d1d2) \
                              .Define("M_d1k",M_d1kst.replace("KST","K"))

    # rdf_mc_gen_n1 = rdf_mc_gen.Define("M2_d1d2", M2_d1d2_mcdtt.replace("Kst","K")) \
    #                           .Define("M2_d1k", M2_d1kst_mcdtt.replace("Kst","K")) \
    #                           .Define("M_d1d2",M_d1d2_mcdtt) \
    #                           .Define("M_d1k",M_d1kst_mcdtt.replace("Kst","K"))

    dbins = 10
    hdalitz_12_1k_dsw = rdf_data_n1.Histo2D(("", f"d1d2_d1k", dbins, 13e6, 24e6, dbins, 5e6, 13e6), "M2_d1d2", "M2_d1k", "sig_sw")
    hdalitz_12_1k_mcrec= rdf_mc_rec_n1.Histo2D(("", f"d1d2_d1k", dbins, 13e6, 24e6, dbins, 5e6, 13e6), "M2_d1d2", "M2_d1k")
    # hdalitz_12_1k_mcgen= rdf_mc_gen_n1.Histo2D(("", f"d1d2_d1k", dbins, 13e6, 24e6, dbins, 5e6, 13e6), "M2_d1d2", "M2_d1k")

    dws_clone = hdalitz_12_1k_dsw.Clone()
    mcrec_clone = hdalitz_12_1k_mcrec.Clone()
    mcgen_clone = mc_gen_hist.Clone()

    mcgen_clone.Scale(1/mcgen_clone.Integral())
    dws_clone.Multiply(mcgen_clone)
    dws_clone.Divide(mcrec_clone)

    print(f"old eff (no bootstrap) {hdalitz_12_1k_mcrec.Integral()}")

    print(f"data integral: {hdalitz_12_1k_dsw.Integral()}")
    print(f"dg/r integral: {dws_clone.Integral()}")
    print(f"new eff: {hdalitz_12_1k_dsw.Integral()/dws_clone.Integral()}")

    # dws_clone.Scale(1/dws_clone.Integral())
    # mcrec_clone.Scale(1/mcrec_clone.Integral())
    # mcgen_clone.Scale(1/mcgen_clone.Integral())

    # dws_clone.Divide(mcrec_clone)

    # for rdf, spec in zip([rdf_data_n1, rdf_data_n1, rdf_mc_rec_n1, rdf_mc_gen_n1],[f"{spec}_data",f"{spec}_data_sw", f"{spec}_mc_rec", f"{spec}_mc_gen"]):
    #     prepdalitz(rdf, spec)

    title = "B^{0} #rightarrow D^{-} D^{0} K^{+}"
    d1 = "D^{-}"
    d2 = "(D^{0} #rightarrow k#pi#pi#pi)"
    d3 = "K^{+}"

    plotdalitz("rtest", title, f"{d1} {d2}", f"{d1} {d3}", dws_clone, 2)


    print(dws_clone.Integral())
#     # dws_clone.SetName("test")
#     # dws_clone.Write()
#     # file_mcweights.Close()
#
#     # for i in range(mc_rec_tree.GetEntries()):
#     #     mc_rec_tree.GetEntry(i)
#     #     weight = dws_clone.GetBinContent(dws_clone.FindBin(mc_rec_tree.B_01, mc_rec_tree.B_02))
#     #     print(weight)
#     	# weight_all += weight
#     	# if tree.B_M>5300:
#     	# 	weight_pass += weight
#
#     ROOT.gInterpreter.ProcessLine("""
#     auto fweight = TFile::Open("MC_Test_Weights.root");
#     auto hweight = fweight->Get<TH2D>("test");
#     """)
# #hweight->GetBinContent
#     rdf_mc_rec_n2 = rdf_mc_rec_n1.Define("M2_d1k_ev", "M2_d1k/1000000") \
#                                  .Define("M2_d1d2_ev", "M2_d1d2/1000000") \
#                                  .Define("dalitz_weight", "hweight->GetBinContent(hweight->FindBin(M2_d1d2_ev, M2_d1k_ev))")
#     rdf_mc_rec_n2.Display(["dalitz_weight", "M2_d1k_ev", "M2_d1d2_ev"],50).Print()
#     weight_all = rdf_mc_rec_n2.Sum("dalitz_weight")
#     # print(weight_all.GetValue())
#     # # weight_pass = df2.Filter("B_M > 5300").Sum("weisght")
#     rdfsnap = rdf_mc_rec_n2.Snapshot(f"DecayTreeTuple", "tfinalw.root")

def build_norm_ws(spec, year, trigger, data_tree):

    norm_bmin = 5230
    norm_bmax = 5330
    nws = ROOT.RooWorkspace(f"{spec}")
    nws.factory(f"B_DTF_M[{norm_bmin},{norm_bmax}]")
    nws.factory(f"Gaussian::{spec}_signal(B_DTF_M, mean_{spec}[5279, 5270, 5290], width_{spec}[2, 0.0, 10.0])")
    nws.factory(f"Exponential:{spec}_bkg(B_DTF_M, c0_n[0, -5, 5])")
    nws.factory(f"SUM::{spec}_fit(n_{spec}_signal[100,0,100000]*{spec}_signal,n_{spec}_bkg[100,0,100000]*{spec}_bkg)")

    # nws.factory("Gaussian::{spec}_a(B_DTF_M, mean_norm[5279, 5270, 5290], width_norm_a[2, 0.0, 10.0])")
    # nws.factory("Gaussian::norm_b(B_DTF_M, mean_norm, width_norm_b[8, 0.0, 10.0])")
    # nws.factory("SUM::norm_signal(norm_a_frac[0.5,0,1]*norm_a, norm_b)")

    b_dtf_m = nws.var("B_DTF_M")
    m_args = ROOT.RooArgSet(b_dtf_m)

    data_set = ROOT.RooDataSet(f"{spec}_data",f"{spec}_data", data_tree, m_args)
    nws.Import(data_set)

    fit = nws.pdf(f"{spec}_fit")
    fit_Result = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    param_to_fix = fit.getVariables()

    # for param in param_to_fix:
    #     if param.GetName() != f"n_{spec}_signal" and param.GetName() != f"n_{spec}_bkg":
    #         print(param)
    #         param.setConstant(ROOT.kTRUE)
    #         print(param)

    nsig = nws.var(f"n_{spec}_signal")
    nbkg = nws.var(f"n_{spec}_bkg")

    MakeSWeights(f"base_norm_files/{spec}_sw_file.root", f"{spec}_sw", data_set, fit,  ROOT.RooArgList(nsig,nbkg))

    # for i in range(0,20):
    #     print(sdata.GetSWeight(i,f"n_{spec}_signal"))
    #     print(sdata.GetSWeight(i,f"n_{spec}_bkg"))
    #     print(sdata.GetSumOfEventSWeight(i))
    #
    # #
    # c1 = ROOT.TCanvas("c1","c1")
    # tframe = b_dtf_m.frame()
    # #DataError(RooAbsData::SumW2)
    # data_set.plotOn(tframe,ROOT.RooFit.Name("Hist"))
    # fit.plotOn(tframe)
    # tframe.Draw()
    # c1.SaveAs("t_old.png")
    # #
    # #
    # m_args.add(nws.var(f"n_{spec}_signal_sw"))

#     print(nws.var(f"n_{spec}_signal_sw"))
#
#     ROOT.gInterpreter.ProcessLine(f"auto histo = nws->var(n_{spec}_signal_sw);")
#
#     # dataw_sw = ROOT.RooDataSet("te","te", data_set, m_args, "", f"n_{spec}_signal_sw")
#     #
#     # c2 = ROOT.TCanvas("c2","c2")
#     # tframe_2 = b_dtf_m.frame()
#     # fit.plotOn(tframe_2, ROOT.RooFit.Components(f"{spec}_signal"))
#     # dataw_sw.plotOn(tframe_2)
#     #
#     # tframe_2.Draw()
#     # c2.SaveAs("t_new.png")
#     #
#     #
#     rdf = RDF(data_tree)
#     rdf.Define("sw",  "histo->GetBinContent(rdfentry_)")
#     # rdf.Print("sw")
#
#
#
#     # # dataw_sw.plotOn(tframe)
#     # nws.Print()
#     # # apply_dalitz_reweight()
#
#     # print(nsig.getVal())
#     # print(sdata.GetYieldFromSWeight(f"n_{spec}_signal"))
#
#
#     # fr_file = ROOT.TFile(f"base_norm_files/fitresult_{spec}_{year}_{trigger}.root","RECREATE") ;
    # fit_Result.Write("rfr")
#     # f.Close()
#     #
#     # nws.Print()
    nws.Import(fit)
    nws.Import(fit_Result)
#     #
#     #
#     #
    nws.writeToFile(f"base_norm_files/{spec}_{year}_{trigger}.root")
    print(f"base_norm_files/{spec}.root", "\n")
#
#
# for spec in ["norm7", "norm8"]:
    # for year in ["2016", "2017", "2018"]:
        # for trigger in ["T", "nTaT"]:
spec = "norm8"
year = "2016"
trigger = "T"
data_file = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/DATA_2021/2016/norm8.root")
data_tree = data_file.Get("DecayTreeTuple_norm8")

mc_file = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC/norm8_norm8_11198007_2016_T_spectrum_filtered.root")
mc_rec_tree = mc_file.Get("DecayTreeTuple")
mc_gen_tree = mc_file.Get("MCDecayTreeTuple")

matt_mc_file = ROOT.TFile(f"norm8_hists.root")
mc_gen_hist = matt_mc_file.Get("h2_true")

# build_norm_ws(spec, year, trigger, data_tree)
# Dalitz_Pre_SWeight(spec, year, trigger, data_tree, mc_rec_tree, mc_gen_tree)
Dalitz_Pre_SWeight(spec, year, trigger, data_tree, mc_rec_tree, mc_gen_hist, mc_gen_tree)
