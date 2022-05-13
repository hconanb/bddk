def build_test_data_sw(run_name, spec, split_flag, fix_flag, smear_flag):
    dws = ROOT.RooWorkspace(run_name)
    bmin = 5200
    bmax = 5350
    dws.factory(f"B_DTF_M[{bmin},{bmax}]")

    dws.factory(f"width_{spec}_gsmear[10, 0.01, 30]")
    dws.factory(f"mean_{spec}_gsmear[0]")
    dws.factory(f"Gaussian::{spec}_gsmear(B_DTF_M, mean_{spec}_gsmear, width_{spec}_gsmear)")

    fix_flag = False
    get_mc_shape_nn(dws, spec, split_flag, fix_flag, smear_flag)
    get_shapes_bkg(spec, "Exponential", dws)

    dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_yield[500,0,10000]*Z_m_p_01_fit, Z_m_p_bkg_yield[100,0,100000]* Z_m_p_spectrum_bkg)")

    b_dtf_m = dws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_dtf_m)
    #############################

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/final_sample/{spec}.root")
    tchain = ROOT.TChain("DecayTreeTuple")
    for file_name in file_list:
        tchain.Add(file_name)

    data = ROOT.RooDataSet(f"{spec}_final_data", f"{spec}_final_data", tchain, data_args)
    model = dws.pdf(f"{spec}_spectrum_all_fit")
    fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    dws.Import(data)
    dws.Import(fit)

    output_base_data = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/{spec}_test_sw_{run_name}.root"
    dws.writeToFile(output_base_data)
    print(f"Wrote dws to: {output_base_data}")
    dws.Print()

    nyield_1 = dws.var(f"Z_m_p_01_yield")
    nyield_bkg = dws.var(f"Z_m_p_bkg_yield")
    yields = ROOT.RooArgSet(nyield_1, nyield_bkg)

    MakeSWeights(f"{analysis_path}/data_fit/sw_files/sw_{spec}_swcomp.root", "SW_tree", data, model, yields)

    frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}_spectrum"))

    title = "D^{-} D^{+} K^{*0}"

    data.plotOn(frame, ROOT.RooFit.Name("data"))
    model.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
    model.plotOn(frame, ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))

    p = ROOT.TCanvas("p1","p1")
    p.cd()

    frame.GetXaxis().SetTitle(f" m({title}) [MeV]")
    frame.Draw()

    save_png(p, f"sw_fit_tests", f"{spec}_rw_test_{run_name}", rpflag = 0)
