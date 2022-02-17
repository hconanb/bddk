import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *


def plot_fit(var, frame, folder, name, fancy_string):

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
    #
    hpull = frame.pullHist(f"{spec}_data", f"{spec}_fit")
    pull = var.frame()
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
    #
    pull.GetXaxis().SetTitle(f"{fancy_string} Mass [MeV]")

    save_png(p, folder, spec, rpflag = 1)

for spec in ["norm7","norm8"]:

    base_file = ROOT.TFile(f"base_norm_files/{spec}_2016_nTaT.root")
    nws = base_file.Get(f"{spec}")

    b_dtf_m = nws.var("B_DTF_M")
    frame = b_dtf_m.frame()

    data_set = nws.data(f"{spec}_data")
    model = nws.pdf(f"{spec}_fit")
    data_set.plotOn(frame, ROOT.RooFit.Name(f"{spec}_data"))

    model.plotOn(
        frame,
        ROOT.RooFit.Components(f"{spec}_fit"),
        ROOT.RooFit.LineStyle(ROOT.kDashed),
        ROOT.RooFit.LineColor(ROOT.kGreen),
        ROOT.RooFit.Name(f"{spec}_fit"),
    )
    model.plotOn(
        frame,
        ROOT.RooFit.Components(f"{spec}_signal"),
        ROOT.RooFit.LineStyle(ROOT.kDashed),
        ROOT.RooFit.LineColor(ROOT.kBlue),
        ROOT.RooFit.Name(f"{spec}_signal"),
    )
    model.plotOn(
        frame,
        ROOT.RooFit.Components(f"{spec}_bkg"),
        ROOT.RooFit.LineStyle(ROOT.kDashed),
        ROOT.RooFit.LineColor(ROOT.kRed),
        ROOT.RooFit.Name(f"{spec}_bkg"),
    )
    if spec == "norm8":
        fancy_string = "D^{-} (D^{0} #rightarrow K#pi#pi#pi) K^{+}"
    if spec == "norm7":
        fancy_string = "D^{0} (D^{0} #rightarrow K#pi#pi#pi) K^{+}"
    d_chi2 = frame.chiSquare(f"{spec}_fit", f"{spec}_data")
    d_mean = nws.obj(f"mean_{spec}")
    d_nyield = nws.var(f"n_{spec}_signal")
    d_nbkg = nws.var(f"n_{spec}_bkg")
    d_std = nws.var(f"width_{spec}")

    dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.AddText(f"mean: {round(d_mean.getValV(),3)}")
    dtpave.AddText(f"nsig: {round(d_nyield.getValV(), 3)}")
    dtpave.AddText(f"nbkg: {round(d_nbkg.getValV(),3)}")
    dtpave.AddText(f"width: {round(d_std.getValV(),3)}")

    frame.addObject(dtpave)

    plot_fit(b_dtf_m, frame, "norm_fit", f"{spec}_fit", fancy_string)
