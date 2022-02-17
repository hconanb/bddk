import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *

def run_data_fit(run_name, gc_onflag, specs):

    lp = 0
    if "z" in specs:
        lp = lp + 4
    if "zz" in specs:
        lp = lp + 6
    if "p" in specs:
        lp = lp + 6
    if "m" in specs:
        lp = lp + 2
    if "st" in specs:
        lp = lp + 3
    if "s" in specs:
        lp = lp + 4

    dws_base_file = ROOT.TFile(f"{analysis_path}/data_fit/base_data_files/{run_name}.root")
    dws = dws_base_file.Get(f"{run_name}")
    if gc_onflag == 1:
        glist = ROOT.RooArgSet()
        for i in range(0, lp):
            gvar = dws.pdf(f"{run_name}_nu{i}_g")
            glist.add(gvar)

    all_data_sets = dws.data("all_data_sets")
    all_cats = dws.cat("all_cats")
    b_dtf_m = dws.var("B_DTF_M")
    all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)

    if "z" in specs:
        z_model = dws.pdf("z_spectrum_all_fit")
        all_fit.addPdf(z_model, "z_spectrum")
    if "zz" in specs:
        zz_model = dws.pdf("zz_spectrum_all_fit")
        all_fit.addPdf(zz_model, "zz_spectrum")
    if "p" in specs:
        p_model = dws.pdf("p_spectrum_all_fit")
        all_fit.addPdf(p_model, "p_spectrum")
    if "m" in specs:
        m_model = dws.pdf("m_spectrum_all_fit")
        all_fit.addPdf(m_model, "m_spectrum")
    if "st" in specs:
        st_model = dws.pdf("st_spectrum_all_fit")
        all_fit.addPdf(st_model, "st_spectrum")
    if "s" in specs:
        s_model = dws.pdf("s_spectrum_all_fit")
        all_fit.addPdf(s_model, "s_spectrum")
    if gc_onflag == 1:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(glist), ROOT.RooFit.Save())
    if gc_onflag == 0:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    fitws = ROOT.RooWorkspace(f"fit_{run_name}")
    fitws.Import(all_data_sets)
    fitws.Import(all_cats)
    fitws.Import(all_fit)
    fitws.Import(nall_fit)
    fitws.Print()
    fitws.writeToFile(f"{analysis_path}/data_fit/fit_data_files/fit_{run_name}.root")

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--run_name')
parser.add_argument('--gc_onflag', type=int)
parser.add_argument('--specs', nargs='+')

args = parser.parse_args()
run_name = args.run_name
gc_onflag = args.gc_onflag
specs  = args.specs

run_data_fit(run_name, gc_onflag, specs)
