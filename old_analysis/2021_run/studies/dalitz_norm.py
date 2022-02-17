import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *

# for spec in ["norm7", "norm8"]:
#     for year in ["2016", "2017", "2018"]:
#         for trigger in ["T", "nTaT"]:

data_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/normalization_fit/base_norm_files/norm8_2016_T.root")
ws = data_file.Get("norm8")

spec = "norm8"

pdf = ws.pdf(f"{spec}_fit")
sigYield = ws.var(f"n_{spec}_signal")
bkgYield = ws.var(f"n_{spec}_bkg")
dataset = ws.data(f"{spec}_data")

param_to_fix = pdf.getVariables()

for param in param_to_fix:
    print(param)
    param.setConstant(ROOT.kTRUE)
    print(param)

ws.Print()


sdata = ROOT.RooStats.SPlot("sData","An SPlot", dataset, pdf, ROOT.RooArgList(sigYield,bkgYield))


# data_file = ROOT.TFile(f"{root_basepath}DATA/{spec}_{year}_{trigger}_spectrum_filtered.root")
# data_tree = data_file.Get("DecayTreeTuple")
