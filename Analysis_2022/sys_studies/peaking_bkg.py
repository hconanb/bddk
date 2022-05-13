import sys
import os

basedir = os.getcwd().split('sys_studies')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

from hep_ml import reweight
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from hep_ml.metrics_utils import ks_2samp_weighted
import matplotlib.pyplot as plt
import warnings

def convertTuple(tup, spec):
    pname = ""
    for p in tup:
        pname += p
    return pname
def invmass(plist):
    """arguments:
    plist -- list of particles, e.g., ["K1", "K2", ...]
    """
    sflag = 0
    for p in plist:
        if sflag == 0:
            e = f"{p}_PE"
            x = f"{p}_PX"
            y = f"{p}_PY"
            z = f"{p}_PZ"
        if sflag == 1:
            e = f"{e}+{p}_PE"
            x = f"{x}+{p}_PX"
            y = f"{y}+{p}_PY"
            z = f"{z}+{p}_PZ"
        sflag = 1
    e2 = f"pow({e},2)"
    x2 = f"pow({x},2)"
    y2 = f"pow({y},2)"
    z2 = f"pow({z},2)"
    s = f"(sqrt({e2}-{x2}-{y2}-{z2}))"
    return s
def plot_mass(data_spec, mc_spec_list, ver):

    if ver == "post_trigger":
        DecayTree_List = ["DecayTreeTuple"]
        path = ver
    if ver == "SIG":
        DecayTree_List = [f"DecayTreeTuple_{ver}"]
        path = "post_signal_and_sb"
    if ver == "SB":
        DecayTree_List = [f"DecayTreeTuple_{ver}"]
        path = "post_signal_and_sb"

    data_file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/{path}/{id_to_spec_dict[data_spec]}.root"
    data_tchain = grab_files_and_chain(data_file_path, DecayTree_List)

    mc_tchain_list = []
    for mc_spec in mc_spec_list:
        mc_file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/{path}/{id_to_spec_dict[mc_spec]}.root"
        mc_tchain = grab_files_and_chain(mc_file_path, DecayTree_List)
        mc_tchain_list.append(mc_tchain)

    data_rdf = RDF(data_tchain)
    # window = id_to_b_values_dict[mc_spec][1]
    # mean = id_to_b_values_dict[mc_spec][0]
    # data_rdf = data_rdf.Filter(f"B_DTF_M < {mean + window} && B_DTF_M > {mean - window}")
    # mc_rdf = mc_rdf.Filter(f"B_DTF_M < {mean + window} && B_DTF_M > {mean - window}")

    nlist = [2, 3, 4, 5]
    if "norm" not in data_spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "KSTH1", "KSTH2"]
        nlist = [2, 3, 4, 5]
    if "norm7" in data_spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "D2H3", "D2H4","K"]
        nlist = nlist + [6]
    if "norm8" in data_spec:
        tl = ["D1H1", "D1H2", "D1H3", "D2H1", "D2H2", "D2H3", "D2H4","K"]
        nlist = nlist + [6, 7]
    if "Z_m_p" in data_spec:
        tl = tl + ["D1H3", "D2H3"]
        nlist = nlist + [6, 7]
    if "P_z_p" in data_spec:
        tl = tl + ["D2H3"]
        nlist = nlist + [6]
    if "M_m_z" in data_spec:
        tl = tl + ["D1H3"]
        nlist = nlist + [6]
    for n in nlist:
        # Final_n_List = []
        all_n_combinations = itertools.combinations(tl, n)
        tupflag = 0
        for tup in all_n_combinations:
            tname = convertTuple(tup, data_spec)

            if len(tup) == 2:
                print(tname)
                if tname in ["D1H1KSTH1","D1H1D2H1","D2H1KSTH1"]:
                    bins = 20
                    xmin = 985
                    xmax = 1055
                else:
                    continue

            elif len(tup) == 3:
                return
                # xmin = 800
                # xmax = 5000

            else:
                return

            c1 = ROOT.TCanvas("c1","c1")


            data_histo = data_rdf.Histo1D((tname, tname, bins, xmin, xmax), tname)
            yh = data_histo.GetBinContent(data_histo.GetMaximumBin())*1.25
            temp_data = data_histo.DrawCopy("E")
            temp_data.SetMarkerColor(ROOT.kRed)

            mc_hist_list = []
            color_list = [ROOT.kBlue, ROOT.kGreen, ROOT.kViolet]
            for chain, color in zip(mc_tchain_list, color_list):
                mc_rdf = RDF(chain)
                mc_histo = mc_rdf.Histo1D((tname, tname, bins, xmin, xmax), tname)
                temp_rec = mc_histo.DrawCopy("E SAME")
                temp_rec.Scale(0.5*temp_data.Integral()/temp_rec.Integral())
                temp_rec.SetMarkerColor(color)
                mc_hist_list.append(temp_rec)

            temp_data.SetMaximum(yh)
            temp_data.GetXaxis().SetTitle(tname)

            legend = ROOT.TLegend(0.70, 0.65, 0.95, 0.90)

            legend.AddEntry(temp_data,f"Data {ver}","lp")
            for mc_spec, temp_rec in zip(mc_spec_list, mc_hist_list):
                legend.AddEntry(temp_rec, f"{mc_spec} {ver}","lp")

            legend.SetTextSize(0.0300)
            legend.SetFillStyle(1001)

            legend.Draw()

            save_pdf(c1, f"peaking_{ver}/{len(tup)}/{data_spec}", tname)
def plot_cut(data_spec, mc_spec_list, ver):
    if ver == "post_trigger":
        DecayTree_List = ["DecayTreeTuple"]
        path = ver
    if ver == "SIG":
        DecayTree_List = [f"DecayTreeTuple_{ver}"]
        path = "post_signal_and_sb"
    if ver == "SB":
        DecayTree_List = [f"DecayTreeTuple_{ver}"]
        path = "post_signal_and_sb"

    data_file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/{path}/{id_to_spec_dict[data_spec]}.root"
    data_tchain = grab_files_and_chain(data_file_path, DecayTree_List)
    data_rdf = RDF(data_tchain)

    if data_spec == "Z_m_p":
        cut_string_1 = "D1H1D2H1 > 1050"
        cut_string_2 = "D2H1KSTH1 > 1050"

        data_rdf_1 = data_rdf.Filter(cut_string_1, "D1H1D2H1")
        data_rdf_2 = data_rdf.Filter(cut_string_2, "D2H1KSTH1")

        print('All stats Data:')
        allCutsReport = data_rdf.Report()
        allCutsReport.Print()

    for mc_spec in mc_spec_list:
        mc_file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/{path}/{id_to_spec_dict[mc_spec]}.root"
        mc_tchain = grab_files_and_chain(mc_file_path, DecayTree_List)
        mc_rdf = RDF(mc_tchain)

        mc_rdf_1 = mc_rdf.Filter(cut_string_1, "D1H1D2H1")
        mc_rdf_2 = mc_rdf.Filter(cut_string_2, "D2H1KSTH1")

        print(f'All stats {mc_spec}:')
        allCutsReport = mc_rdf.Report()
        allCutsReport.Print()



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Apply D Window Cuts Conditions")
    parser.add_argument("--Data_Spec", choices = id_to_spec_dict.keys(), help = 'Spec')
    parser.add_argument("--Version", choices = ["post_trigger","SIG","SB"])
    parser.add_argument('--Comp', action = 'store_true')
    parser.add_argument('--Cut', action = 'store_true')

    args = parser.parse_args()
    data_spec = args.Data_Spec
    mc_spec_list = data_to_mc_dict[data_spec]
    ver = args.Version
    comp = args.Comp
    cut = args.Cut

    if comp:
        plot_mass(data_spec, mc_spec_list, ver)
    if cut:
        plot_cut(data_spec, mc_spec_list, ver)
