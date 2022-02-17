PE_d1d2 = "(pow(D1_PE + D2_PE,2))"
PX_d1d2 = "(pow(D1_PX + D2_PX,2))"
PY_d1d2 = "(pow(D1_PY + D2_PY,2))"
PZ_d1d2 = "(pow(D1_PZ + D2_PZ,2))"
M2_d1d2 = f"{PE_d1d2} - ({PX_d1d2} + {PY_d1d2} + {PZ_d1d2})"
M_d1d2 = f"sqrt({M2_d1d2})"
def draw_distributions(MC, DATA, MC_weights, Data_weights, spec, mcid, savename):

    c1 = ROOT.TCanvas("c1","c1")

    for id, column in enumerate(dalitz_columns, 1):

        bins = 20

        data_np = DATA[column].to_numpy()
        data_df = ROOT.RDF.MakeNumpyDataFrame({column:data_np, "data_weights": Data_weights})

        MC_np = MC[column].to_numpy()
        MC_df = ROOT.RDF.MakeNumpyDataFrame({column:MC_np, "mc_weights": MC_weights})

        xmin = MC_df.Min(column).GetValue()
        xmax = MC_df.Max(column).GetValue()

        hist_data = data_df.Histo1D((f"data_{column}", f"data_{column}", bins, xmin, xmax), column, "data_weights")
        hist_MC = MC_df.Histo1D((f"MC_{column}", f"MC_{column}", bins, xmin, xmax), column,  "mc_weights")

        hist_MC.Scale(hist_data.Integral()/hist_MC.Integral())

        hist_data.SetMarkerColor(ROOT.kRed)
        hist_MC.SetMarkerColor(ROOT.kBlue)


        temp_data = hist_data.DrawCopy("E")
        temp_mc = hist_MC.DrawCopy("Same")

        temp_data.GetXaxis().SetTitle(f"m^{{2}}({column}) [GeV^{{2}}]")
        temp_data.GetXaxis().SetLimits(temp_data.GetXaxis().GetXmin()/1000000,temp_data.GetXaxis().GetXmax()/1000000)
        temp_mc.GetXaxis().SetLimits(temp_mc.GetXaxis().GetXmin()/1000000,temp_mc.GetXaxis().GetXmax()/1000000)

        legend = ROOT.TLegend(0.70, 0.60, 0.95, 0.95)
        legend.AddEntry(temp_data, f"Data_SWeighted","lp")
        legend.AddEntry(temp_mc, f"MC","lp")
        legend.Draw()

        save_png(c1, f"{savename}", f"{spec}_{mcid}_{column}", rpflag = 0)

        # print('KS over ', column, ' = ', ks_2samp_weighted(MC[column], DATA[column],
        #                          weights1=MC_weights, weights2=Data_weights))

        # xlim = numpy.percentile(numpy.hstack([target[column]]), [0.01, 99.99])
        # plt.subplot(2, 3, id)
        # plt.hist(original[column], weights=new_original_weights, range=xlim, **hist_settings)
        # plt.hist(target[column], weights=new_target_weights, range=xlim, **hist_settings)
        # plt.title(column)
        # print('KS over ', column, ' = ', ks_2samp_weighted(original[column], target[column],
        #                                  weights1=new_original_weights, weights2=new_target_weights))
        # plt.savefig(f'{savename}.png')
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
