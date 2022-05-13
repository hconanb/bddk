import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021")
from essentials import *
from rootutils import residualPlot
ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")
ROOT.gStyle.SetPalette(ROOT.kBird)


def build_norm_ws(spec, year, trigger, data_tree, strat):

    norm_mean = 5280
    norm_bmin = norm_mean - 40
    norm_bmax = norm_mean + 40

    nws = ROOT.RooWorkspace(f"{spec}")
    nws.factory(f"B_DTF_M[{norm_bmin},{norm_bmax}]")
    if strat == "DG":
        nws.factory(f"Gaussian::{spec}_signal_a(B_DTF_M, mean_{spec}[5279, 5270, 5290], width_{spec}_a[10, 0.1, 15.0])")
        nws.factory(f"Gaussian::{spec}_signal_b(B_DTF_M, mean_{spec}, width_{spec}_b[5, 0.1, 10.0])")
        nws.factory(f"SUM::{spec}_signal({spec}_a_frac[0.5,0,1]*{spec}_signal_a, {spec}_signal_b)")
    if strat == "G":
        nws.factory(f"Gaussian::{spec}_signal(B_DTF_M, mean_{spec}[5279, 5270, 5290], width_{spec}_a[2, 0.1, 30.0])")
    nws.factory(f"Exponential:{spec}_bkg(B_DTF_M, c0_n[0, -5, 5])")
    nws.factory(f"SUM::{spec}_fit(n_{spec}_signal[100,0,100000]*{spec}_signal,n_{spec}_bkg[100,0,100000]*{spec}_bkg)")

    b_dtf_m = nws.var("B_DTF_M")
    m_args = ROOT.RooArgSet(b_dtf_m)

    if year == "all":
        data_set = ROOT.RooDataSet(f"{spec}_data",f"{spec}_data", data_tree[0], m_args)
        data_set_nTaT = ROOT.RooDataSet(f"{spec}_data",f"{spec}_data", data_tree[1], m_args)
        data_set.append(data_set_nTaT)
    else:
        data_set = ROOT.RooDataSet(f"{spec}_data",f"{spec}_data", data_tree, m_args)

    nws.Import(data_set)

    model = nws.pdf(f"{spec}_fit")
    fit_Result = model.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    nws.Import(fit_Result)

    if year == "all":
        nyield_1 = nws.var(f"n_{spec}_signal")
        nyield_2 = nws.var(f"n_{spec}_bkg")
        yields = ROOT.RooArgSet(nyield_1, nyield_2)

        MakeSWeights(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/normalization_fit/sw_files/sw_{spec}.root", "SW_tree", data_set, model, yields)

    output_base_data = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/normalization_fit/fit_files/{spec}_{year}_{trigger}.root"
    nws.writeToFile(output_base_data)
    print(f"Wrote nws to: {output_base_data}")
    # dws.Print()

    frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}"), ROOT.RooFit.Bins(40))

    data_set.plotOn(frame, ROOT.RooFit.Name(f"{spec}_data"))
    model.plotOn(frame,ROOT.RooFit.Components(f"{spec}_fit"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kGreen),ROOT.RooFit.Name(f"{spec}_fit"))
    model.plotOn(frame,ROOT.RooFit.Components(f"{spec}_signal"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kBlue),ROOT.RooFit.Name(f"{spec}_signal"))
    model.plotOn(frame,ROOT.RooFit.Components(f"{spec}_bkg"),ROOT.RooFit.LineStyle(ROOT.kDashed),ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Name(f"{spec}_bkg"))
    if spec == "norm8":
        fancy_string = "D^{-} (D^{0} #rightarrow K#pi#pi#pi) K^{+}"
    if spec == "norm7":
        fancy_string = "#bar{D^{0}} (D^{0} #rightarrow K#pi#pi#pi) K^{+}"

    d_chi2 = frame.chiSquare(f"{spec}_fit", f"{spec}_data")
    d_mean = nws.obj(f"mean_{spec}")
    d_nyield = nws.var(f"n_{spec}_signal")
    d_nbkg = nws.var(f"n_{spec}_bkg")
    d_std = nws.var(f"width_{spec}_a")

    d_nyield_err = d_nyield.getPropagatedError(fit_Result)

    dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.AddText(f"mean: {round(d_mean.getValV(),3)}")
    dtpave.AddText(f"nsig: {round(d_nyield.getValV(), 3)} #pm {round(d_nyield_err,3)}")
    dtpave.AddText(f"nbkg: {round(d_nbkg.getValV(),3)}")
    dtpave.AddText(f"width: {round(d_std.getValV(),3)}")

    frame.addObject(dtpave)

    p = residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)

    leg = ROOT.TLegend(0.75,0.5,0.95,0.87)
    # leg.SetFillColor(kWhite)
    # leg.SetLineColor(kWhite)
    leg.AddEntry(f"{spec}_data","Data", "P");
    leg.AddEntry(frame.findObject(f"{spec}_fit"),"Signal + background","L");
    leg.AddEntry(frame.findObject(f"{spec}_bkg"),"Background only", "L");
    leg.AddEntry(frame.findObject(f"{spec}_signal"),"Signal only", "L");

    frame.Draw()
    leg.Draw()
    p.pb.cd()

    hpull = frame.pullHist(f"{spec}_data", f"{spec}_fit")
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

    pull.GetXaxis().SetTitle(f"{fancy_string} Mass [MeV]")

    output_name = f"{spec}_{year}_{trigger}_final_fit"
    folder = "norm_fits"
    save_png(p, folder, output_name, rpflag = 1)
# Split for yields
for spec in ["norm8"]:
    for year in ["2016", "2017", "2018"]:
        for trigger in ["T", "nTaT"]:
            data_file = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/final_sample/{spec}.root")
            data_tree_base = data_file.Get(f"DecayTreeTuple_{trigger}")
            build_norm_ws(spec, year, trigger, data_tree_base, "G")

# #Combined for rw
# for spec in ["norm7", "norm8"]:
#
#     data_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{spec}.root")
#
#     data_tree_T_base = ROOT.TChain(f"DecayTreeTuple_T")
#     data_tree_nTaT_base = ROOT.TChain(f"DecayTreeTuple_nTaT")
#
#     for file_name in data_file_list:
#         data_tree_T_base.Add(file_name)
#         data_tree_nTaT_base.Add(file_name)
#
#     build_norm_ws(spec, "all","all", [data_tree_T_base, data_tree_nTaT_base], "G")
