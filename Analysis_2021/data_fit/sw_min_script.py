import ROOT as ROOT
import glob as glob
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")

analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"
ids_and_shapes = {
                    "Z_m_p" : [("Z_m_p_01_fit", "DG"), ("Z_m_p_0203_fit" , "BGEP"), ("Z_m_p_04_fit" , "BGEP")],
}

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
    ROOT.RooMsgService.instance().addStream(ROOT.RooFit.DEBUG, ROOT.RooFit.Topic(ROOT.RooFit.Contents), ROOT.RooFit.OutputFile("rfsmear_debug.log"))
    #projDeps = ROOT.RooArgSet(), useWeights = False, copyDataSet = True, newName = "test", fitToarg5 =  ROOT.RooFit.PrintLevel(3)
    sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields)

    print("Check SWeights:")
    for y in yields:
        oval = y.getVal()
        sval = sData.GetYieldFromSWeight(y.GetName())
        print(f"Yield of {y.GetName()} is {oval}")
        print(f"from sWeights it is {sval}")
        if not (0.9995 < oval / sval < 1.0005):
            raise Exception("sWeighted yield should match")
    for i in range(10):
        for y in yields:
            print(f"    {y.GetName()} Weight {sData.GetSWeight(i, y.GetName())}")
        totw = sData.GetSumOfEventSWeight(i)
        print(f"Total Weight {totw}")
        if not (0.9995 < totw < 1.0005):
            raise Exception("sum of sWeight should be 1")
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
        assert [x.GetName() for x in swvars] == swnames
        # print(swnames)
        # print([x.GetName() for x in swvars]) # check sorting worked
        # set values
        for val, var in zip(swvals, swvars):
            val[0] = var.getVal()
        # fill values
        nt.Fill()
    nt.Write()
    nf.Close()

def get_shapes_bkg(spec, flag, ws):
    if flag == "Exponential":
        ws.factory(f"Exponential:{spec}_spectrum_bkg(B_DTF_M, c0_{spec}[0, -0.01, 0.01])")

def get_mc_shape(dws, base_spec, split_flag, fix_flag, smear_flag):
    if split_flag:
        slist = ids_and_shapes[f"{base_spec}_split"]
    if not split_flag:
        slist = ids_and_shapes[f"{base_spec}"]

    for tuple in slist:
        mc_spec = tuple[0]
        strat = tuple[1]
        mcws_base = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_fit/fit_mc_files/{mc_spec}_{strat}.root")
        mcws = mcws_base.Get(f"fit_ws")
        if fix_flag:
            vars = mcws.allVars()
            for i in vars:
                if i.GetName() != "B_DTF_M" and "mean" not in i.GetName():
                # if "mean" not in i.GetName():# and "width" not in i.GetName():
                    temp = mcws.var(i.GetName())
                    temp.setConstant(True)
                    print(f"{temp} is constant")
        if smear_flag:
            mc_pdf = mcws.pdf(f"{mc_spec}")
            dws.factory(
                f"Gaussian::{mc_spec}_gsmear(B_DTF_M, mean_{base_spec}_gsmear, width_{base_spec}_gsmear)"
            )
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable(mc_spec, f"{mc_spec}_pc"))
            dws.factory(
                f"FCONV::{mc_spec}(B_DTF_M, {mc_spec}_pc, {mc_spec}_gsmear)"
            )
        if not smear_flag:
            mc_pdf = mcws.pdf(f"{mc_spec}")
            dws.Import(mc_pdf)

dws = ROOT.RooWorkspace("ws")

bmin = 4800
bmax = 5600
spec = "Z_m_p"
dws.factory(f"B_DTF_M[{bmin},{bmax}]")
tv = dws.var("B_DTF_M")
# tv.setBins(10000, "cache")

dws.factory(f"width_{spec}_gsmear[22, 1, 30]")
dws.factory(f"mean_{spec}_gsmear[0]")
split_flag = False
fix_flag = True
smear_flag = True

get_mc_shape(dws, spec, split_flag, fix_flag, smear_flag)
get_shapes_bkg(spec, "Exponential", dws)

dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_yield[500,0,10000]*Z_m_p_01_fit, Z_m_p_0203_yield[500,0,10000]*Z_m_p_0203_fit, Z_m_p_04_yield[500,0,10000]*Z_m_p_04_fit, Z_m_p_bkg_yield[1000,0,100000]*Z_m_p_spectrum_bkg)")

dws.Print()
model = dws.pdf("Z_m_p_spectrum_all_fit")
# data = model.generate(tv, 10000)

file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/final_sample/{spec}.root")
tchain = ROOT.TChain("DecayTreeTuple")
for file_name in file_list:
    tchain.Add(file_name)

data = ROOT.RooDataSet(f"{spec}_final_data", f"{spec}_final_data", tchain, ROOT.RooArgSet(tv))
# dws.Import(temp_ds)

fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

# c1 = ROOT.TCanvas("c1","c1")
# frame = tv.frame()
# data.plotOn(frame)
# model.plotOn(frame)
# frame.Draw()
# c1.SaveAs("t.png")

nyield_1 = dws.obj(f"Z_m_p_01_yield")
nyield_2 = dws.obj(f"Z_m_p_0203_yield")
nyield_3 = dws.obj(f"Z_m_p_04_yield")
nyield_bkg = dws.obj(f"Z_m_p_bkg_yield")
yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, nyield_bkg)
# sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields, projDeps = ROOT.RooArgSet(), useWeights = False, copyDataSet = True, newName = "test", fitToarg5 =  ROOT.RooFit.PrintLevel(3))

MakeSWeights("t.root", "tree", data, model, yields)

# t = ROOT.RooRealVar("t", "t", 5200, 5350)
#
# # Construct landau(t,ml,sl)
# m1 = ROOT.RooRealVar("m1", "mean g1", 5280, 5270, 5290)
# s1 = ROOT.RooRealVar("s1", "sigma g1", 5, 0, 30)
# g1 = ROOT.RooGaussian("g1", "g1", t, m1, s1)
#
# # Construct gauss(t,mg,sg)
# m2 = ROOT.RooRealVar("m2", "mg2", 0)
# s2 = ROOT.RooRealVar("s2", "sg2", 3, 0.1, 10)
# g2 = ROOT.RooGaussian("g2", "g2", t, m2, s2)
#
# m3 = ROOT.RooRealVar("m3", "mean g3", 5230, 5220, 5240)
# s3 = ROOT.RooRealVar("s3", "sigma g3", 10, 0, 30)
# g3 = ROOT.RooGaussian("g3", "g3", t, m3, s3)
# #
# # # Construct gauss(t,mg,sg)
# # m4 = ROOT.RooRealVar("m4", "mg4", 0)
# # s4 = ROOT.RooRealVar("s4", "sg4", 1)
# # g4 = ROOT.RooGaussian("g4", "g4", t, m4, s4)
#
# # Construct convolution pdf
# # ---------------------------------------
#
# # Set #bins to be used for FFT sampling to 10000
# t.setBins(10000, "cache")
#
# dws.Import(t)
# dws.factory(f"Exponential:bkg(t, c0[0.005, -0.01, 0.01])")
# bkg3 = dws.pdf("bkg")
#
# # Construct landau (x) gauss
# g1g2 = ROOT.RooFFTConvPdf("g1g2", "g1g2", t, g1, g2)
# g3g4 = ROOT.RooFFTConvPdf("g3g2", "g3g2", t, g3, g2)
#
# y1 = ROOT.RooRealVar("y1", "y1", 1000, 0, 10000)
# y2 = ROOT.RooRealVar("y2", "y2", 2000, 0, 10000)
# bkgy = ROOT.RooRealVar("y3", "y4", 1500, 0, 10000)
#
# shapes = ROOT.RooArgList(g1g2, g3g4, bkg3)
# yields = ROOT.RooArgList(y1,y2,bkgy)
# model = ROOT.RooAddPdf("model", "model", shapes, yields)
# # Sample, fit and plot convoluted pdf
# # ----------------------------------------------------------------------
#
# # Sample 1000 events in x from gxlx
# data = model.generate(t, 10000)
# model.fitTo(data)
# sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields, projDeps = ROOT.RooArgSet(), useWeights = False, copyDataSet = True, newName = "test", fitToarg5 =  ROOT.RooFit.PrintLevel(3))

# c1 = ROOT.TCanvas("c1","c1")
# frame = t.frame()
# data.plotOn(frame)
# frame.Draw()
# c1.SaveAs("t.png")
