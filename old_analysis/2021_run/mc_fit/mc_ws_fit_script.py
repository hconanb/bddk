import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *

def run_mc_fit(ws_name, shape, wflag = 0):

        mcws_base_file = ROOT.TFile(f"{analysis_path}/mc_fit/base_mc_files/{ws_name}.root")
        mcws = mcws_base_file.Get(f"{ws_name}")
        b_dtf_m = mcws.var("B_DTF_M")
        data_set = mcws.data(f"{ws_name}_events")

        fitws = ROOT.RooWorkspace(f"fit_{ws_name}")
        fitws.Import(data_set)

        fit = mcws.pdf(f"{ws_name}_{shape}_fit")
        if wflag == 0:
            fit_pdf = fit.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
            fitws.Import(fit_pdf)
            fitws.Import(fit)
            fitws.Print()
            fitws.writeToFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_{ws_name}_{shape}.root")


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--ws_name')
parser.add_argument('--shape')

args = parser.parse_args()
ws_name = args.ws_name
shape = args.shape

run_mc_fit(ws_name, shape)
