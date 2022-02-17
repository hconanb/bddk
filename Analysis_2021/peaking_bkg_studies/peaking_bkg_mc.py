import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
    # rdf_data_all = rdf_data_all.Define("B_M013_fix","if (KSTH1_ID > 0) return B_M013; else return B_M023;") \
    #                  # .Define("B_VtxChi2_013_fix","if (KSTH1_ID > 0) return B_VtxChi2_013; else return B_VtxChi2_023;") \
    #                  # .Define("KSTPI_IPCHI2", "log(B_ENDVERTEX_CHI2 - B_VtxChi2_013_fix)") \


def sweight_data(infiles, outflag):

    run_name = f"sw_test_{outflag}"
    fit_strat = "dp_f_nn"

    specs = ["Z_z_z"]
    #["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]
    ids_and_shapes = {
                         "Z_m_p" : [["01", "G", 5280, 15],
                                   ["02", "BG", 5130, 15],
                                   ["04", "G", 4980, 15]],

                         "Z_z_z": [["09", "DG", 5280, 10],
                               ["0710", "BG", 5130, 15],
                               ["040812", "BG", 4975, 20]],

                         "P_z_p": [["05", "G", 5280, 10],
                               ["020607","BGEP", 5130, 15],
                               ["0408", "BG", 4975, 20]],

                         "M_m_z": [["03", "BG", 5130, 15],
                               ["04", "BG", 4975, 20]],

                         "P_z_pst": [["07", "BG", 5280, 10],
                               ["0408","DG", 5130, 15]],

                         "Zs_sm_p" : [["13","G", 5370, 10],
                                     ["14","G", 5225, 15],
                                     ["15","G", 5225, 15],
                                     ["16","BG", 5075, 15]]}

    #build fit for data -> Save
    #run fit for data -> Save S weighted Data
    #plot sweighted data for a peak vs relevant mc

    dws = ROOT.RooWorkspace(run_name)
    dws.factory(f"B_DTF_M[{bmin},{bmax}]")

    for spec in specs:
        get_shapes_bkg(spec, "Exponential", dws)
        for s_list in ids_and_shapes[spec]:
            get_free_shapes(dws, spec, fit_strat, s_list)

    for i in range(1, 18):
        dws.factory(f"nny_{i}[500,0,10000]")

    # dws.factory("SUM::Z_m_p_spectrum_all_fit(nny_1*Z_m_p_01_fit,nny_2*Z_m_p_02_fit,nny_3*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
    dws.factory("SUM::Z_z_z_spectrum_all_fit(nny_4*Z_z_z_09_fit,nny_5*Z_z_z_0710_fit,nny_6*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
    # dws.factory("SUM::P_z_p_spectrum_all_fit(nny_7*P_z_p_05_fit,nny_8*P_z_p_020607_fit,nny_9*P_z_p_0408_fit,n_P_z_P_z_p_bkg[100,0,10000]*P_z_p_spectrum_bkg)")
    # dws.factory("SUM::M_m_z_spectrum_all_fit(nny_10*M_m_z_03_fit,nny_11*M_m_z_04_fit,n_M_m_z_bkg[100,0,10000]*M_m_z_spectrum_bkg)")
    # dws.factory("SUM::P_z_pst_spectrum_all_fit(nny_12*P_z_pst_07_fit,nny_13*P_z_pst_0408_fit,n_P_z_st_bkg[100,0,1000]*P_z_pst_spectrum_bkg)")
    # dws.factory("SUM::s_spectrum_all_fit(nny_14*Zs_sm_p_13_fit,nny_15*Zs_sm_p_14_fit,nny_16*Zs_sm_p_15_fit,nny_17*Zs_sm_p_16_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")
    dws.Print()

    b_dtf_m = dws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_dtf_m)
    #
    for spec in specs:
        file_list = infiles
        tchain = ROOT.TChain("DecayTreeTuple")
        for file_name in file_list:
            tchain.Add(file_name)

        data_set = ROOT.RooDataSet(f"{spec}_data", f"{spec}_data", tchain, data_args)
        dws.Import(data_set)

        presw_data_set = dws.data(f"{spec}_data")
        fit_model = dws.pdf(f"{spec}_spectrum_all_fit")
        gofit = fit_model.fitTo(presw_data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

        if spec == "Z_m_p":
            nyield_1 = dws.obj(f"nny_1")
            nyield_2 = dws.obj(f"nny_2")
            nyield_3 = dws.obj(f"nny_3")
            n_bkg = dws.obj(f"n_{spec}_bkg")
        if spec == "Z_z_z":
            nyield_1 = dws.obj(f"nny_4")
            nyield_2 = dws.obj(f"nny_5")
            nyield_3 = dws.obj(f"nny_6")
            n_bkg = dws.obj(f"n_{spec}_bkg")
        if spec == "P_z_p":
            nyield_1 = dws.obj(f"nny_7")
            nyield_2 = dws.obj(f"nny_8")
            nyield_3 = dws.obj(f"nny_9")
            n_bkg = dws.obj(f"n_{spec}_bkg")
        if spec == "M_m_z":
            nyield_1 = dws.obj(f"nny_10")
            nyield_2 = dws.obj(f"nny_11")
            n_bkg = dws.obj(f"n_{spec}_bkg")
        if spec == "P_z_pst":
            nyield_1 = dws.obj(f"nny_12")
            nyield_2 = dws.obj(f"nny_13")
            n_bkg = dws.obj(f"n_{spec}_bkg")
        if spec == "Zs_sm_p":
            nyield_1 = dws.obj(f"nny_14")
            nyield_2 = dws.obj(f"nny_15")
            nyield_3 = dws.obj(f"nny_16")
            nyield_4 = dws.obj(f"nny_17")
            n_bkg = dws.obj(f"n_{spec}_bkg")

        cs = ROOT.TCanvas("cs","cs")
        frame = b_dtf_m.frame()
        presw_data_set.plotOn(frame)
        fit_model.plotOn(frame)

        if spec == "Z_m_p":
            fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
        if spec == "Z_z_z":
            fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_09_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
            fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_0710_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
            fit_model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_040812_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))

        # vars = dws.allVars()
        # mean_vars = ROOT.RooArgSet()
        # for i in vars:
        #     if f"{spec}" in i.GetName() and "mean" in i.GetName():
        #         mean_vars.add(i)
        #
        # dtpave = ROOT.TPaveText(0.50, 0.65, 0.70, 0.85, "NB NDC")
        # dtpave.SetFillStyle(0)
        # dtpave.AddText(f"Mean 0: {round(d_chi2, 3)}")
        # dtpave.AddText(f"Mean 1: {round(d_chi2, 3)}")
        # dtpave.AddText(f"Mean 2: {round(d_chi2, 3)}")



        frame.Draw()
        save_png(cs, f"fit_tests", f"{run_name}", rpflag = 0)
        if spec == "Z_m_p" or spec == "P_z_p" or spec == "Z_z_z":
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, n_bkg)

        sData = ROOT.RooStats.SPlot("sData","An SPlot", presw_data_set, fit_model, yields)
        # dws.Import(presw_data_set, ROOT.RooFit.Rename(f"swdata_{spec}"))
        # dws.writeToFile(f"root_files/sw_{run_name}.root")

        swnames = sorted([f"{x.GetName()}_sw" for x in yields])
        print(f"weights: {swnames}")
        # create output file
        nf = ROOT.TFile.Open(f"root_files/sw_{outflag}_{spec}.root", "recreate")
        outtreename = "testtree"
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
        for i in range(presw_data_set.numEntries()):
            # get vars
            swvars = sorted(
                [x for x in presw_data_set.get(i) if x.GetName() in swnames],
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

    # dws.writeToFile(f"root_files/{run_name}.root")

def plot_vtx_chi2(datatree, mc_spec_list, spec):

    rdf_data_all = RDF(datatree)
    rdf_data_all = rdf_data_all.Define("nextnny_1_sw", f"sw_tree.nny_4_sw")
    rdf_data_all = rdf_data_all.Define("nextnny_2_sw", f"sw_tree.nny_5_sw")
    rdf_data_all = rdf_data_all.Define("nextnny_3_sw", f"sw_tree.nny_6_sw")
    rdf_mc_base_list = []
    hist_list = []

    for mc_spec in mc_spec_list:
        mc_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_MC/*/post_veto/{mc_spec}*")
        mc_tree_chain = ROOT.TChain(f"DecayTreeTuple")
        for mcfile in mc_file_list:
            mc_tree_chain.Add(mcfile)

        min_ddk = 4200
        max_ddk = 5300

        rdf_mc_base = RDF(mc_tree_chain)
        hist_temp = rdf_mc_base.Histo1D((f"KSTPI_IPCHI2_{mc_spec}", f"KSTPI_IPCHI2_{mc_spec}", 100, -5, 5), 'KSTPI_IPCHI2')
        hist_list.append(hist_temp.GetPtr())

    hist_bmddkc_hi2_data = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", 100, -5, 5), 'KSTPI_IPCHI2')

    c1 = ROOT.TCanvas("c1","c1")

    hist_bmddkc_hi2_data.SetTitle("B - DDK")
    hist_bmddkc_hi2_data.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    # hist_bmddkc_hi2_data.GetYaxis().SetRangeUser(0., 200.)

    hist_list[0].Scale(hist_bmddkc_hi2_data.Integral()/hist_list[0].Integral())
    hs_list = [hist_bmddkc_hi2_data.GetPtr(), hist_list[0]]
    hs = DrawStack(ROOT.gPad, hs_list, legend=(0.20, 0.5, 0.40, 0.9), drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"comp_b_ddk_log_{spec}", rpflag = 0)


    # hist_bmddkc_hi2_mc2.Scale(hist_bmddkc_hi2_data.Integral()/hist_bmddkc_hi2_mc2.Integral())
    # hist_bmddkc_hi2_mc2.SetLineColor(ROOT.kBlue)
    # hist_bmddkc_hi2_mc2.Draw("SAME HIST")
    #
    # hist_bmddkc_hi2_mc3.Scale(hist_bmddkc_hi2_data.Integral()/hist_bmddkc_hi2_mc3.Integral())
    # hist_bmddkc_hi2_mc3.SetLineColor(ROOT.kGreen)
    # hist_bmddkc_hi2_mc3.Draw("SAME HIST")

    # # if spec == "Z_m_p":
    # #     legend = ROOT.TLegend(0.25, 0.4, 0.35, 0.70)
    # #     legend.AddEntry("KSTPI_IPCHI2_data", "Data","l")
    # #     legend.AddEntry("KSTPI_IPCHI2_mc1", "D-D+ MC","l")
    # #     legend.AddEntry("KSTPI_IPCHI2_mc2", "D*-D+/D-D*+ MC","l")
    # #     legend.AddEntry("KSTPI_IPCHI2_mc3", "D*-D*+","l")
    # #     legend.SetTextSize(0.040)
    # #     legend.Draw()
    # if spec == "Z_z_z":
    #     legend = ROOT.TLegend(0.25, 0.4, 0.35, 0.70)
    #     legend.AddEntry("KSTPI_IPCHI2_data", "Data","l")
    #     legend.AddEntry("KSTPI_IPCHI2_mc1", "D0D0 MC","l")
    #     legend.AddEntry("KSTPI_IPCHI2_mc2", "D*0D0/D0D*0 MC","l")
    #     legend.AddEntry("KSTPI_IPCHI2_mc3", "D*0D*0","l")
    #     legend.SetTextSize(0.040)
    #     legend.Draw()

def submass_harris(datatree, spec, mc_spec_list):

    print(f"reading trees for {spec}")

    rdf_data_all = RDF(datatree)
    rdf_data_all = rdf_data_all.Define("DSTM_DM",'D1H1D1H2KSTH2 - D1H1D1H2')
    rdf_data_cut = rdf_data_all.Filter("DSTM_DM > 150")
    rdf_data_cut2 = rdf_data_cut.Filter("D1H1D1H2D2H1D2H2KSTH1 < 5100")

    rdf_mc_base_list = []
    rdf_mc_wcut_list = []
    rdf_mc_wcut2_list = []
    hist_list = []
    hist_list2 = []

    min_ddk = 4200
    max_ddk = 5300

    # for mc_spec in mc_spec_list:
    #     mc_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_MC/*/pre_d/{mc_spec}*")
    #     print(mc_file_list)
    #     mc_tree_chain = ROOT.TChain(f"DecayTreeTuple_{mc_spec}")
    #     for mcfile in mc_file_list:
    #         mc_tree_chain.Add(mcfile)
    #
    #     rdf_mc_base_1 = RDF(mc_tree_chain)
    #     rdf_mc_base = rdf_mc_base_1.Define("DSTM_DM",'D1H1D1H2KSTH2 - D1H1D1H2')
    #     rdf_mc_base_list.append(rdf_mc_base)
    #     rdf_mc_wcut = rdf_mc_base.Filter("DSTM_DM > 150", f"{mc_spec}_1")
    #     rdf_mc_wcut_list.append(rdf_mc_wcut)
    #     rdf_mc_wcut2 = rdf_mc_wcut.Filter("D1H1D1H2D2H1D2H2KSTH1 < 5100",f"{mc_spec}_2")
    #     rdf_mc_wcut2_list.append(rdf_mc_wcut2)
    #
    #     # print('All stats:')
    #     # allCutsReport = rdf_mc_base_1.Report()
    #     # allCutsReport.Print()
    #     if spec == "Z_z_z":
    #         hist_temp = rdf_mc_wcut.Histo1D((f"DDK_{mc_spec}", f"DDK_{mc_spec}", 100, min_ddk, max_ddk), 'D1H1D1H2D2H1D2H2KSTH1')
    #     if spec == "P_z_p":
    #         hist_temp = rdf_mc_wcut.Histo1D((f"DDK_{mc_spec}", f"DDK_{mc_spec}", 100, min_ddk, max_ddk), 'D1H1D1H2D2H1D2H2KSTH1D2H3')
    #
    #     hist_temp2 = rdf_mc_wcut2.Histo1D((f"DDKpi_{mc_spec}", f"DDKpi_{mc_spec}", 100, 4200, 5400), 'B_DTF_M')
    #
    #     hist_list.append(hist_temp.GetPtr())
    #     hist_list2.append(hist_temp2.GetPtr())


    dmd_bins = 100
    dmd_min = 140
    dmd_max = 180
    dmd_bw = (dmd_max- dmd_min)/dmd_bins

    hist_dstmd_data = rdf_data_all.Histo1D((f"Kpi(D1)pi(KST)-Kpi", f"Kpi(D1)pi(KST)-Kpi", dmd_bins, dmd_min, dmd_max), 'DSTM_DM')

    bins = 70
    bw = (max_ddk - min_ddk)/100

    if spec == "Z_z_z":

        hist_ddk_pre_cut = rdf_data_all.Histo1D((f"DDK_pc", f"DDK_pc", bins, min_ddk, max_ddk), 'D1H1D1H2D2H1D2H2KSTH1')
        hist_ddk_post_cut = rdf_data_cut.Histo1D((f"DDK_pc", f"DDK_pc", bins, min_ddk, max_ddk), 'D1H1D1H2D2H1D2H2KSTH1')
        ddk_title = "M(#bar{D^{0}}D^{0}K^{+}) (MeV/c^{2})"

    if spec == "P_z_p":
        hist_ddk_pre_cut = rdf_data_all.Histo1D((f"DDK_pc", f"DDK_pc", bins, min_ddk, max_ddk), 'D1H1D1H2D2H1D2H2KSTH1D2H3')
        hist_ddk_post_cut = rdf_data_cut.Histo1D((f"DDK_pc", f"DDK_pc", bins, min_ddk, max_ddk), 'D1H1D1H2D2H1D2H2KSTH1D2H3')
        ddk_title = "M(#bar{D^{0}}D^{+}K^{+}) (MeV/c^{2})"

    ### Plot of dst-d precut
    c1 = ROOT.TCanvas("c1","c1")
    hist_dstmd_data.GetXaxis().SetTitle("M(#bar{D^{0}}#pi^{-}) - M(#bar{D^{0}}) (MeV/c^{2})")
    hist_dstmd_data.GetYaxis().SetTitle(f"Events/({dmd_bw})")
    hist_dstmd_data.Draw()
    c1.Update()

    pmax = ROOT.gPad.GetUymax()
    print(pmax)

    line_cut = ROOT.TLine(150, 0, 150, pmax)
    line_cut.SetLineColor(ROOT.kRed)
    line_cut.SetLineWidth(5)

    line_cut.Draw()
    c1.Update()
    save_png(c1, "his_test", f"precut_{spec}_dstdm", rpflag = 0)

    ### plot of ddk post cut with MC samples
    c2 = ROOT.TCanvas("c2","c2")
    hist_ddk_post_cut.GetXaxis().SetTitle(ddk_title)
    hist_ddk_post_cut.GetYaxis().SetTitle(f"Events/({bw})")

    hist_ddk_pre_cut.SetLineColor(ROOT.kGreen)
    hist_ddk_post_cut.SetLineColor(ROOT.kRed)

    hist_ddk_post_cut.Draw("")
    hist_ddk_pre_cut.Draw("SAME")

    # for next_hist in hist_list:
    #     next_hist.Scale(0.5*(hist_ddk_post_cut.Integral()/next_hist.Integral()))
    # hist_list.insert(0,hist_ddk_post_cut.GetPtr())
    # hs = DrawStack(ROOT.gPad, hist_list, legend=(0.70, 0.7, 0.90, 0.9), drawopts="nostack plc pmc")
    # col_list = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kYellow, ROOT.kViolet, ROOT.kOrange]

    legend = ROOT.TLegend(0.75, 0.6, 0.85, 0.80)
    legend.AddEntry(hist_ddk_pre_cut.GetPtr(), "Before Cut", "l")
    legend.AddEntry(hist_ddk_post_cut.GetPtr(), "After Cut", "l")
    legend.SetTextSize(0.040)
    legend.Draw()

    save_png(c2, "his_test", f"postcut_{spec}_", rpflag = 0)

    # c3 = ROOT.TCanvas("c3","c3")
    # hist_ddk_post_cut.GetXaxis().SetTitle(ddk_title)
    # hist_ddk_post_cut.Draw()
    # save_png(c3, "his_test", f"postcut_nomc_{spec}", rpflag = 0)

    # ############
    # c3 = ROOT.TCanvas("c3","c3")
    # hist_ddk_post_cut.GetXaxis().SetTitle("M(D^{0}#bar{D^{0}}K^{+}) (MeV/c^{2})")
    # hist_ddk_post_cut.Draw()
    # # col_list = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kYellow, ROOT.kViolet, ROOT.kOrange]
    # # legend = ROOT.TLegend(0.65, 0.6, 0.75, 0.90)
    # hist_list[3].Scale(0.7*(hist_ddk_post_cut.Integral()/hist_list[3].Integral()))
    # hist_list[3].SetLineColor(ROOT.kRed)
    # hist_list[3].Draw("SAME HIST")
    # legend.SetTextSize(0.040)
    # # legend.Draw()
    # save_png(c3, "his_test", f"1mc_postcut_{spec}", rpflag = 0)
    #
    #
    # c4 = ROOT.TCanvas("c4","c4")
    # hist_ddkpi_post_cut.GetXaxis().SetTitle("m(D^{0}#bar{D^{0}}K^{*0}) DTF (MeV/c^{2})")
    # hist_ddkpi_post_cut.Draw()
    # save_png(c4, "his_test", f"ddkpi_postcut2_{spec}", rpflag = 0)


    # c3 = ROOT.TCanvas("c3","c3")
    #
    # hist_ddk_post_cut.SetTitle("M(#bar{D^{0}}#pi^{-})")
    # hist_ddk_post_cut.GetXaxis().SetTitle("m(D^{0}#bar{D^{0}}K^{+}) (MeV/c^{2})")
    # hist_ddk_post_cut.Draw()
    # col_list = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kYellow, ROOT.kViolet, ROOT.kOrange]
    # legend = ROOT.TLegend(0.65, 0.6, 0.75, 0.90)
    # for next_hist, col, mcs in zip(hist_list, col_list, mc_spec_list):
    # # next_hist = hist_list[0]
    # # # temp_hist = next_hist.GetPtr().Clone()
    #
    #     next_hist.Scale(0.7*(hist_ddk_post_cut.Integral()/next_hist.Integral()))
    #     next_hist.SetLineColor(col)
    #     next_hist.Draw("SAME HIST")
    #     legend.AddEntry(next_hist, mcs, "l")
    #
    # legend.SetTextSize(0.040)
    # legend.Draw()
    #
    # save_png(c3, "his_test", f"ddk_pc_{spec}", rpflag = 0)

def plot_new_mass(data_tree):
    rdf_data_all = RDF(data_tree)
    # rdf_data_all = rdf_data_all.Define("DSTM_DM",'D1H1D1H2KSTH2 - D1H1D1H2')
    rdf_data_all = rdf_data_all.Define("DSTM_DM",'D2H1D2H2KSTH2 - D2H1D2H2')

    hist_dstmd_data = rdf_data_all.Histo1D((f"Kpi(D1)pi(KST)-Kpi", f"Kpi(D1)pi(KST)-Kpi", 100, 0, 1000), 'DSTM_DM')

    c1 = ROOT.TCanvas("c1","c1")
    hist_dstmd_data.GetXaxis().SetTitle("M(#bar{D^{0}}#pi^{-}) - M(#bar{D^{0}}) (MeV/c^{2})")
    hist_dstmd_data.GetYaxis().SetTitle(f"Events/(100)")
    hist_dstmd_data.Draw()

    save_png(c1, "his_test", f"precut_dstdm", rpflag = 0)


# zzz_mc_spec_list = ["04_Z_z_z_11198023","07_Z_z_z_12197045","08_Z_z_z_12197423","09_Z_z_z_11196019","10_Z_z_z_11196413","12_Z_z_z_11196414"]

# pzp_mc_spec_list = ["02_P_z_p_11198005","04_P_z_p_11198410","05_P_z_p_12197023","06_P_z_p_12197410","07_P_z_p_12197400","08_P_z_p_12197401"]
#
# pzp_data_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/pre_d/P_z_p.root")
# pzp_data_tree_chain = ROOT.TChain(f"DecayTreeTuple_P_z_p")
# for datafile in pzp_data_file_list:
#     pzp_data_tree_chain.Add(datafile)
#
# zzz_data_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/pre_d/Z_z_z.root")
# zzz_data_tree_chain = ROOT.TChain(f"DecayTreeTuple_Z_z_z")
# for datafile in zzz_data_file_list:
#     zzz_data_tree_chain.Add(datafile)

# submass_harris(zzz_data_tree_chain, "Z_z_z", zzz_mc_spec_list)
# submass_harris(pzp_data_tree_chain, "P_z_p", pzp_mc_spec_list)

mmz_data_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/final_sample/M_m_z.root")
mmz_data_tree_chain = ROOT.TChain(f"DecayTreeTuple")
for datafile in mmz_data_file_list:
    mmz_data_tree_chain.Add(datafile)
plot_new_mass(mmz_data_tree_chain)

# sweight_data(glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_veto/Z_z_z_postveto.root"), "noecut")
# sweight_data(glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_veto/Z_z_z_postvetoecut.root"), "ecut")
# zzz_data_file_list_pv = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_veto/Z_z_z_postveto.root")
# zzz_data_tree_chain_pv = ROOT.TChain(f"DecayTreeTuple")
# for datafile in zzz_data_file_list_pv:
#     zzz_data_tree_chain_pv.Add(datafile)
#
# sw_file = ROOT.TFile(f"root_files/sw_sw_test_Z_z_z.root")
# sw_tree = sw_file.Get(f"testtree")
# zzz_data_tree_chain_pv.AddFriend(sw_tree, "sw_tree")
#
# plot_vtx_chi2(zzz_data_tree_chain_pv, zzz_mc_spec_list, "Z_z_z")
# plot_vtx_chi2(zzz_data_tree_chain_pv, zzz_mc_spec_list, "Z_z_z")
