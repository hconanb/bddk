import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

def plot_norm(nws, name, trigger_flag):

    print(nws)
    b_dtf_m = nws.var("B_DTF_M")
    ds = nws.data(f"{name}_{trigger_flag}")
    norm_fit = nws.pdf("norm_fit")

    frame = b_dtf_m.frame(ROOT.RooFit.Range("myrange"))
    ds.plotOn(frame, ROOT.RooFit.CutRange("myrange"), ROOT.RooFit.Name("data"))
    nData = ds.sumEntries("", "myrange")

    norm_fit.plotOn(frame, ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Components(f"norm_bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("bkg"))
    norm_fit.plotOn(frame, ROOT.RooFit.NormRange("myrange"), ROOT.RooFit.Normalization(nData, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Range("myrange"), ROOT.RooFit.Components(f"norm_fit"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))

    c1 = ROOT.TCanvas("c1","c1")
    frame.Draw()
    c1.SaveAs(f"plots/{name}_{trigger_flag}.png")
def save_yields(nws, name, p):
    n_norm_signal = nws.var("n_norm_signal")
    n = n_norm_signal.getValV()
    e = n_norm_signal.getError()
    f = open(p, "a")
    f.write(f"n_{name.split('_')[0]} {name.split('_')[1]} {n} {e}\n")
    f.close()

ws_list = ["norm8_ToT","norm8_T","norm8_nTaT","norm7_ToT","norm7_T","norm7_nTaT"]

p = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/norm_yield.txt"
if os.path.exists(p):
    os.remove(p)

for name in ws_list:
    pfile= ROOT.TFile(f"normalization/{name}.root")
    ws = pfile.Get(name)
    print(ws)
    plot_norm(ws, name.split("_")[0],  name.split("_")[1])
    save_yields(ws, name, p)
