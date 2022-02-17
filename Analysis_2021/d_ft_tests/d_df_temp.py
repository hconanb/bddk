import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

def temp_norm_store(type):
    spec = "norm8"
    for year in ["2016","2017","2018"]:
        if type == "DATA":
            file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_data/{year}/*/*/ntuple.root")
            tree_name = f"data_{spec}_Data/DecayTreeTuple"
        if type == "MC":
            file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/norm8_norm8_11198007/11198007_{year}_*/*/ntuple.root")
            tree_name = f"norm8_norm8_11198007/DecayTreeTuple"
        tree_chain = ROOT.TChain(tree_name)
        new_file_list = []
        for file_name in file_list:
            tempo = ROOT.TFile(file_name)
            dirlist = [b.GetTitle() for b in tempo.GetListOfKeys()]
            if tree_name.split("/")[0] in dirlist:
                new_file_list.append(file_name)
            tempo.Close()
        for next_file_name in new_file_list:
            tree_chain.Add(next_file_name)

        temp_rdf = RDF(tree_chain)
        temp_rdf = temp_rdf.Define("B_DTF_M", "B_dtf_M[0]")
        temp_rdf = temp_rdf.Filter("B_DTF_M >= 5200 && B_DTF_M <= 5360")
        temp_rdf = temp_rdf.Filter("abs(D1_M - 1869.62) < 60 && abs(D2_M - 1864.840) < 60")
        # hist_2016 = temp_rdf.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", 100, 5200, 5360), 'B_DTF_M')
        #
        # cs = ROOT.TCanvas("cs","cs")
        # hist_2016.Draw()
        # save_png(cs, "dfd_test", f"testbm", rpflag = 0)
        # break
# temp_norm_store()
    #     frame = b_m.frame()
    #     data_set.plotOn(frame)
    #     fit.plotOn(frame)

        clist = temp_rdf.GetColumnNames()
        outputfile = f"{spec}_{year}_{type}_test.root"
        print(f"Starting snapshot for {outputfile}")
        rdfsnap =  temp_rdf.Snapshot(f"{spec}_{year}_{type}_test", outputfile, clist)
        print("adsfh")
# temp_norm_store("DATA")
#
tree_list = []
# #
file_2016 = ROOT.TFile("norm8_2016_DATA_test.root")
tree_2016 = file_2016.Get("norm8_2016_DATA_test")

file_2017 = ROOT.TFile("norm8_2017_DATA_test.root")
tree_2017 = file_2017.Get("norm8_2017_DATA_test")

file_2018 = ROOT.TFile("norm8_2018_DATA_test.root")
tree_2018 = file_2018.Get("norm8_2018_DATA_test")
#
for tree, year in zip([tree_2016, tree_2017, tree_2018],["2016", "2017", "2018"]):
    spec = "norm8"
    ##############################################
    ws = ROOT.RooWorkspace(f"sw_norm8_{year}")
    ws.factory(f"B_DTF_M[5200,5360]")
    # ws.factory(f"Chebychev:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[0.,-3,3],c1_{spec}[0.,-3,3], c2_{spec}[0.,-3,3]}})")
    ws.factory(f"Exponential:{spec}_spectrum_bkg(B_DTF_M, c0_{spec}[0, -5, 5])")
    ws.factory(f"Gaussian::{spec}_fit_1(B_DTF_M, mean_1[5280,5270,5290], width_1[5.0,0.0,30.0])")
    # ws.factory(f"Gaussian::{spec}_fit_2(B_DTF_M, mean_1, width_3[5.0,0.0,100.0])")
    ws.factory(f"SUM::{spec}_all_fit(n_1[100,0,100000]*{spec}_fit_1, n_bkg[100,0,100000]*{spec}_spectrum_bkg)")
    # ws.factory(f"SUM::{spec}_all_fit(n_1[100,0,1000]*{spec}_fit_1,n_2[100,0,1000]*{spec}_fit_2,n_3[100,0,1000]*{spec}_fit_3,n_bkg[100,0,10000]*{spec}_spectrum_bkg)")
    b_m = ws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_m)
    fit = ws.pdf(f"{spec}_all_fit")
    data_set = ROOT.RooDataSet(f"{spec}_data_{year}", f"{spec}_data_{year}", tree, data_args)
    all_fit = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    ws.Import(all_fit)
    fitresult = ws.obj(f"fitresult_{spec}_all_fit_{spec}_data")

    cs = ROOT.TCanvas("cs","cs")
    frame = b_m.frame()
    data_set.plotOn(frame, ROOT.RooFit.Binning(80))
    fit.plotOn(frame)
    fit.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit_1"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
    # fit.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit_2"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
    fit.plotOn(frame, ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
    frame.Draw()
    save_png(cs, f"fit_tests", f"{spec}_{year}", rpflag = 0)


    nyield_1 = ws.obj(f"n_1")
    print(nyield_1.getVal())

    n_bkg = ws.obj(f"n_bkg")
    yields = ROOT.RooArgSet(nyield_1, n_bkg)

    sData = ROOT.RooStats.SPlot("sData","An SPlot", data_set, fit, yields)

    swnames = sorted([f"{x.GetName()}_sw" for x in yields])
    print(f"weights: {swnames}")
    # create output file
    nf = ROOT.TFile.Open(f"sw_{year}_{spec}.root", "recreate")
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
    for i in range(data_set.numEntries()):
        # get vars
        swvars = sorted(
            [x for x in data_set.get(i) if x.GetName() in swnames],
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

    sw_file = ROOT.TFile(f"sw_{year}_{spec}.root")
    sw_tree = sw_file.Get(f"testtree")
    tree.AddFriend(sw_tree, "sw_tree")
    tree_list.append(tree)
#############################################

file_2016_mc = ROOT.TFile("norm8_2016_mc_test.root")
tree_2016_mc = file_2016_mc.Get("norm8_2016_mc_test")

file_2017_mc = ROOT.TFile("norm8_2017_mc_test.root")
tree_2017_mc = file_2017_mc.Get("norm8_2017_mc_test")

file_2018_mc = ROOT.TFile("norm8_2018_mc_test.root")
tree_2018_mc = file_2018_mc.Get("norm8_2018_mc_test")

rdf_mc_2016 = RDF(tree_2016_mc)
rdf_mc_2017 = RDF(tree_2017_mc)
rdf_mc_2018 = RDF(tree_2017_mc)

rdf_2016 = RDF(tree_list[0])
rdf_2017 = RDF(tree_list[1])
rdf_2018 = RDF(tree_list[2])

rdf_data_list = [rdf_2016, rdf_2017, rdf_2018]
rdf_mc_list = [rdf_mc_2016, rdf_mc_2017, rdf_mc_2018]

for rdf_data, rdf_mc, year in zip(rdf_data_list, rdf_mc_list, ["2016","2017","2018"]):

    rdf_data = rdf_data.Define("log_D1_x2","log(D1_FDCHI2_ORIVX)")
    rdf_mc = rdf_mc.Define("log_D1_x2","log(D1_FDCHI2_ORIVX)")

    rdf_data = rdf_data.Define("log_D2_x2","log(D2_FDCHI2_ORIVX)")
    rdf_mc = rdf_mc.Define("log_D2_x2","log(D2_FDCHI2_ORIVX)")

    rdf_data = rdf_data.Define("ny_1_sw", f"sw_tree.n_1_sw")

    for dname in ["D1","D2"]:

        nbins = 25
        min = -1
        max = 15

        nbinsx2 = 25
        minx2 = -6
        maxx2 = 15

        hist_data = rdf_data.Histo1D((f"{dname}_FD_ORIVX_{spec}_data", f"{dname}_FD_ORIVX_{spec}_data", nbins, min, max), f"{dname}_FD_ORIVX", "ny_1_sw")
        hist_mc = rdf_mc.Histo1D((f"{dname}_FD_ORIVX_{spec}_mc", f"{dname}_FD_ORIVX_{spec}_mc", nbins, min, max), f"{dname}_FD_ORIVX")

        hist_datax2  = rdf_data.Histo1D((f"{dname}_FDCHI2_ORIVX_{spec}_data", f"{dname}_FDCHI2_ORIVX_{spec}_data", nbinsx2, minx2, maxx2), f"log_{dname}_x2", "ny_1_sw")
        hist_mcx2  = rdf_mc.Histo1D((f"{dname}_FDCHI2_ORIVX_{spec}_mc", f"{dname}_FDCHI2_ORIVX_{spec}_mc", nbinsx2, minx2, maxx2), f"log_{dname}_x2")

        hist_list_1 = [hist_data.GetPtr(), hist_mc.GetPtr()]
        histx2_list_1 = [hist_datax2.GetPtr(), hist_mcx2.GetPtr()]

        pname_list = ["x","x2"]

        for list,pname in zip([hist_list_1, histx2_list_1],pname_list):
            list[1].Scale(list[0].Integral()/list[1].Integral())
            list[0].Sumw2()
            c1 = ROOT.TCanvas("c1","c1")
            hs = DrawStack(ROOT.gPad, list, drawopts="nostack plc pmc")
            save_png(c1, "dfd_test", f"{dname}_{spec}_{year}_{pname}", rpflag = 0)
