import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run')
from essentials import *
from uncertainties import correlated_values
from uncertainties.umath import sqrt
import root_pandas as rp
import itertools

branch_dict = {

"D1H1D1H2D2H1D2H2KSTH1" : "K+pi-(D0bar) K-pi+(D0) K+(Kst0)",
"D1H1D1H2D2H1D2H2KSTH2" : "K+pi-(D0bar) K-pi+(D0) pi-(Kst0)",
"D1H1D1H2D2H1D2H2KSTH1KSTH2" : "K+pi-(D0bar) K-pi+(D0) K+pi-(Kst0)",
"D1H1D1H2KSTH1" : "K+pi-(D0bar)  K+(Kst0)",
"D2H1D2H2KSTH1" : "K-pi+(D0)  K+(Kst0)",
"D1H1D1H2KSTH2" : "K+pi-(D0bar)  pi-(Kst0)",
"D2H1D2H2KSTH2" : "K-pi+(D0)  pi-(Kst0)",

"B_M01_Fixed" : "K+pi-(D0bar)  K-pi+(D0)",
"B_M02_Fixed" : "K+pi-(D0bar)  K+(Kst0)",
"B_M03_Fixed" : "K+pi-(D0bar)  pi-(Kst0)",
"B_M12_Fixed" : "K-pi+(D0)  K+(Kst0)",
"B_M13_Fixed" : "K-pi+(D0)  pi-(Kst0)",
"B_M23_Fixed" : "K+(Kst0)  pi-(Kst0)",

"B_M012_Fixed" : "K+pi-(D0bar) K-pi+(D0) K+(Kst0)",
"B_M013_Fixed" : "K+pi-(D0bar) K-pi+(D0) pi-(Kst0)",
"B_M0123_Fixed" : "K+pi-(D0bar) K-pi+(D0) K+(Kst0) pi-(Kst0)",

"04" : "B0 -> D*- ->D0bar pi- D*+ -> D0 pi+ Kst0",
"07" : "B+ -> D0bar D*+->D0 pi+ Kst0",
"08" : "B+ -> D*0bar D*+->D0 pi+ Kst0",
"09" : "B0 -> D0bar D0 Kst0",
"10" : "B0 -> D*0bar D0 K*0 + D0bar D*0 Kst0",
"12" : "B0 -> D*0bar D*0 Kst0",
"Data" : "D0bar D0 Kst0",

}

def rdf_plot(rdf, column, nbins, min, max, title):
    id_name = title.split("_")[0]
    hist = rdf.Histo1D((f"{column}", f"{column}", nbins, min, max), f"{column}")
    bincan = ((max-min)/nbins)
    c1 = ROOT.TCanvas("c1","c1")
    hist_plot = hist.DrawCopy()
    if "IP" in column:
        column = column.replace("IPChi2_","M")
    hist_plot.GetXaxis().SetTitle(branch_dict[column])
    hist_plot.GetYaxis().SetTitle(f"Candidates / ({bincan} MeV)")
    # ROOT.gStyle.SetOptTitle(0);
    print(branch_dict[id_name])
    pt = ROOT.TPaveText(.70,.85,1,1, "BNDC")
    pt.AddText(branch_dict[id_name])
    pt.SetFillColorAlpha (0, 1)
    # pt.SetFillStyle(4050)
    # ltitle = ROOT.TPaveLabel(.11,.95,.35,.99, branch_dict[id_name],"brndc")
    pt.Draw()
    c1.Update()
    save_png(c1, "InvM_SecondCut", title)


def plot_invm():
    zz_invm_file_list = glob.glob(f"rootfiles/*.root")
    for file_name in zz_invm_file_list:
        id_name = file_name.split("/")[1].split("_")[0]
        file = ROOT.TFile(file_name)
        dtt = file.Get("DTT")
        rdf_base = RDF(dtt)
        rdf = rdf_base.Filter("D1H1D1H2KSTH2 > 2050")
        rdf = rdf.Filter("D1H1D1H2D2H1D2H2KSTH1 < 5080")

        print(id_name)
        if "Data" in file_name:
            rdf_plot(rdf, "D1H1D1H2KSTH1", 100, 1800, 3000, f"{id_name}_d1K")
            rdf_plot(rdf, "D2H1D2H2KSTH1", 100, 1800, 3000, f"{id_name}_d2K")
            rdf_plot(rdf, "D1H1D1H2KSTH2", 100, 1800, 3000, f"{id_name}_d1pi")
            rdf_plot(rdf, "D2H1D2H2KSTH2", 100, 1800, 3000, f"{id_name}_d2pi")
        rdf_plot(rdf, "D1H1D1H2D2H1D2H2KSTH1", 100, 3700, 5400, f"{id_name}_K_invM")
        rdf_plot(rdf, "D1H1D1H2D2H1D2H2KSTH2", 100, 3700, 5400, f"{id_name}_pi_invM")
        rdf_plot(rdf, "D1H1D1H2D2H1D2H2KSTH1KSTH2", 100, 4800, 5400, f"{id_name}_6_track_invM")

def plot_subm():
    zz_id_list = ["04", "07", "08", "09", "10", "12"]
    for id in zz_id_list:
        tc = ROOT.TChain(f"DecayTreeTuple")
        zz_subm_file_list = glob.glob(f"../submass_fix/rootfiles/{id}*.root")
        for file in zz_subm_file_list:
            tc.Add(file)
        rdf_base = RDF(tc)
        rdf = rdf_base.Filter("B_M03_Fixed > 2050 && B_M012_Fixed < 5080")
        # rdf_plot(rdf, "B_M02_Fixed", 100, 1800, 3000, f"{id}_D0barK")
        # rdf_plot(rdf, "B_M03_Fixed", 100, 1800, 3000, f"{id}_D0barpi")
        # rdf_plot(rdf, "B_M12_Fixed", 100, 1800, 3000, f"{id}_D0K")
        # rdf_plot(rdf, "B_M13_Fixed", 100, 1800, 3000, f"{id}_D0pi")
        # rdf_plot(rdf, "B_M012_Fixed", 100, 3700, 5400, f"{id}_K_subM")
        # rdf_plot(rdf, "B_M013_Fixed", 100, 3700, 5400, f"{id}_pi_subM")
        oldn = rdf_base.Count().GetValue()
        newn = rdf.Count().GetValue()
        print(f"{id} -  Old Rec N: {oldn}  New Rec N: {newn}")
def plot_ipchi2():
    zz_id_list = ["04", "07", "08", "09", "10", "12"]
    for id in zz_id_list:
        tc = ROOT.TChain(f"DecayTreeTuple")
        zz_subm_file_list = glob.glob(f"../submass_fix/rootfiles/{id}*.root")
        for file in zz_subm_file_list:
            tc.Add(file)
        rdf_base = RDF(tc)
        rdf = rdf_base.Filter("B_M03_Fixed > 2050")
        ipmin = 0
        ipmax = 600
        # rdf_plot(rdf, "B_IPChi2_01_Fixed", 100, ipmin, ipmax,  f"{id}_D0barD0")
        # rdf_plot(rdf, "B_IPChi2_02_Fixed", 100, ipmin, ipmax, f"{id}_D0barK")
        # rdf_plot(rdf, "B_IPChi2_03_Fixed", 100, ipmin, ipmax, f"{id}_D0barpi")
        # rdf_plot(rdf, "B_IPChi2_12_Fixed", 100, ipmin, ipmax, f"{id}_D0K")
        # rdf_plot(rdf, "B_IPChi2_13_Fixed", 100, ipmin, ipmax, f"{id}_D0pi")
        # rdf_plot(rdf, "B_IPChi2_23_Fixed", 100, ipmin, ipmax, f"{id}_Kpi")
        #
        # rdf_plot(rdf, "B_IPChi2_012_Fixed", 100, ipmin, ipmax, f"{id}_D0barD0K")
        # rdf_plot(rdf, "B_IPChi2_013_Fixed", 100, ipmin, ipmax, f"{id}_D0barD0pi")
        # rdf_plot(rdf, "B_IPChi2_023_Fixed", 100, ipmin, ipmax, f"{id}_D0barKpi")
        # rdf_plot(rdf, "B_IPChi2_123_Fixed", 100, ipmin, ipmax, f"{id}_D0barKpi")

        # rdf_plot(rdf, "B_IPChi2_0123_Fixed", 100, ipmin, ipmax, f"{id}_D0barpi")

plot_invm()
plot_subm()
# plot_ipchi2()
