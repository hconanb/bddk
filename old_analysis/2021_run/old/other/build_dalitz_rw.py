import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run')
from essentials import *

####### Dalitz Weights
PE_d1d2 = "(pow(D1_PE + D2_PE,2))"
PX_d1d2 = "(pow(D1_PX + D2_PX,2))"
PY_d1d2 = "(pow(D1_PY + D2_PY,2))"
PZ_d1d2 = "(pow(D1_PZ + D2_PZ,2))"
M2_d1d2 = f"{PE_d1d2} - ({PX_d1d2} + {PY_d1d2} + {PZ_d1d2})"
M_d1d2 = f"sqrt({M2_d1d2})"

PE_d1kst = "(pow(D1_PE + KST_PE,2))"
PX_d1kst = "(pow(D1_PX + KST_PX,2))"
PY_d1kst = "(pow(D1_PY + KST_PY,2))"
PZ_d1kst = "(pow(D1_PZ + KST_PZ,2))"
M2_d1kst = f"{PE_d1kst} - ({PX_d1kst} + {PY_d1kst} + {PZ_d1kst})"
M_d1kst = f"sqrt({M2_d1kst})"

PE_d2kst = "(pow(D2_PE + KST_PE,2))"
PX_d2kst = "(pow(D2_PX + KST_PX,2))"
PY_d2kst = "(pow(D2_PY + KST_PY,2))"
PZ_d2kst = "(pow(D2_PZ + KST_PZ,2))"
M2_d2kst = f"{PE_d2kst} - ({PX_d2kst} + {PY_d2kst} + {PZ_d2kst})"
M_d2kst = f"sqrt({M2_d2kst})"
#########################################MCDecayTreeTuple
PE_d1d2_mcdtt = "(pow(C0_TRUEP_E + C1_TRUEP_E,2))"
PX_d1d2_mcdtt = "(pow(C0_TRUEP_X + C1_TRUEP_X,2))"
PY_d1d2_mcdtt = "(pow(C0_TRUEP_Y + C1_TRUEP_Y,2))"
PZ_d1d2_mcdtt = "(pow(C0_TRUEP_Z + C1_TRUEP_Z,2))"
M2_d1d2_mcdtt = f"{PE_d1d2_mcdtt} - ({PX_d1d2_mcdtt} + {PY_d1d2_mcdtt} + {PZ_d1d2_mcdtt})"
M_d1d2_mcdtt = f"sqrt({M2_d1d2_mcdtt})"

PE_d1kst_mcdtt = "(pow(C0_TRUEP_E + Kst_TRUEP_E,2))"
PX_d1kst_mcdtt = "(pow(C0_TRUEP_X + Kst_TRUEP_X,2))"
PY_d1kst_mcdtt = "(pow(C0_TRUEP_Y + Kst_TRUEP_Y,2))"
PZ_d1kst_mcdtt = "(pow(C0_TRUEP_Z + Kst_TRUEP_Z,2))"
M2_d1kst_mcdtt = f"{PE_d1kst_mcdtt} - ({PX_d1kst_mcdtt} + {PY_d1kst_mcdtt} + {PZ_d1kst_mcdtt})"
M_d1kst_mcdtt = f"sqrt({M2_d1kst_mcdtt})"

PE_d2kst_mcdtt = "(pow(C1_TRUEP_E + Kst_TRUEP_E,2))"
PX_d2kst_mcdtt = "(pow(C1_TRUEP_X + Kst_TRUEP_X,2))"
PY_d2kst_mcdtt = "(pow(C1_TRUEP_Y + Kst_TRUEP_Y,2))"
PZ_d2kst_mcdtt = "(pow(C1_TRUEP_Z + Kst_TRUEP_Z,2))"
M2_d2kst_mcdtt = f"{PE_d2kst_mcdtt} - ({PX_d2kst_mcdtt} + {PY_d2kst_mcdtt} + {PZ_d2kst_mcdtt})"
M_d2kst_mcdtt = f"sqrt({M2_d2kst_mcdtt})"
##########################################################
dbins = 20
d1kstmin, d1kstmax = 2600, 3800
d2kstmin, d2kstmax = 2600, 3800
d1d2min, d1d2max = 3600, 4500

def plotdalitz(title, xname, yname, histo, dim):

    tltxt = ROOT.TLatex()

    if dim == 2:
        histo.GetXaxis().SetTitle(f"m^{{2}}({xname}) [GeV^{{2}}]")
        histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
        histo.GetYaxis().SetTitle(f"m^{{2}}({yname}) [GeV^{{2}}]")
        histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)

    # if dim == 2:
    #     histo.GetXaxis().SetTitle(f"m^{{2}}({xname}) [GeV^{{2}}]")
    #     histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
    #     histo.GetYaxis().SetTitle(f"m^{{2}}({yname}) [GeV^{{2}}]")
    #     histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
    #
    # if dim == 1:
    #     histo.GetXaxis().SetTitle(f"m({xname}) [GeV]")
    #     histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000,histo.GetXaxis().GetXmax()/1000)
    #     bw = histo.GetXaxis().GetBinWidth(5)
    #     histo.GetYaxis().SetTitle(f"{yname} / ({bw} GeV)")
    #     histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000,histo.GetYaxis().GetXmax()/1000)
    # if dim == 2 and rflag == 1:
    #     histo.GetXaxis().SetTitle(f"ratio")
    #     # histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
    #     histo.GetYaxis().SetTitle(f"ratio")
    #     # histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
    # if dim == 1 and rflag == 1:
    #     histo.GetXaxis().SetTitle(f"ratio")
    #     # histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000,histo.GetXaxis().GetXmax()/1000)
    #     bw = histo.GetXaxis().GetBinWidth(5)
    #     histo.GetYaxis().SetTitle(f"{yname} / ({bw} GeV)")
    #     # histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000,histo.GetYaxis().GetXmax()/1000)

    c1 = ROOT.TCanvas("c1","c1")
    c1.SetRightMargin(0.15)
    c1.SetLeftMargin(0.9)
    c1.SetTopMargin(0.1)
    c1.SetBottomMargin(0.9)
    ht1 = histo.DrawCopy("COLZ")
    tltxt.DrawLatexNDC(.4,.925, title)
    pname = histo.GetTitle()
    save_png(c1, "daliz_plots", pname, rpflag=0)

def make_dalitz(file, type_flag):
    nfile = ROOT.TFile(file)
    tree = nfile.Get("DecayTreeTuple")
    if type_flag == "MC":
        spec = file.split("/")[-1].split("_spectrum_filtered.root")[0]
        spec_name = file.split("/")[-1].split("_")[1]
    if type_flag == "DATA":
        spec = file.split("/")[-1].split("_spectrum_filtered.root")[0]
        spec_name = spec
        spec = "DATA_" + spec
    rdf_base = RDF(tree)
    rdf_next = rdf_base.Define("M2_d1d2", M2_d1d2) \
                       .Define("M2_d2kst", M2_d2kst) \
                       .Define("M2_d1kst", M2_d1kst) \
                       .Define("M_d1d2",M_d1d2) \
                       .Define("M_d2kst",M_d2kst) \
                       .Define("M_d1kst",M_d1kst) \


    # hdalitz_1k_p = rdf_next.Histo1D(("", f"{spec}_d1kst_projection", dbins, d1kstmin, d1kstmax), "M_d1kst")
    # hdalitz_2k_p = rdf_next.Histo1D(("", f"{spec}_d2kst_projection", dbins, d2kstmin, d2kstmax), "M_d2kst")
    # hdalitz_12_p = rdf_next.Histo1D(("", f"{spec}_d1d2_projection", dbins, d1d2min, d1d2max), "M_d1d2")

    hdalitz_12_2k = rdf_next.Histo2D(("", f"{spec}_d1d2_d2kst", dbins,  d1d2min**2, d1d2max**2, dbins, d2kstmin**2, d2kstmax**2), "M2_d1d2", "M2_d2kst")
    hdalitz_12_1k = rdf_next.Histo2D(("", f"{spec}_d1d2_d1kst", dbins,  d1d2min**2, d1d2max**2, dbins, d1kstmin**2, d1kstmax**2), "M2_d1d2", "M2_d1kst")
    hdalitz_1k_2k = rdf_next.Histo2D(("", f"{spec}_d1kst_d2kst", dbins, d1kstmin**2,d1kstmax**2, dbins, d2kstmin**2, d2kstmax**2), "M2_d1kst", "M2_d2kst")

    d3 = "K^{*0}"
    if spec_name == "z" or spec_name == "z_mc_rec" or spec_name == "z_mc_gen":
        title = "B^{0} #rightarrow D^{-} D^{+} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{+}"
    if spec_name == "zz":
        title = "B^{0} #rightarrow #bar{D^{0}} D^{0} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{0}"
    if spec_name == "p":
        title = "B^{+} #rightarrow #bar{D^{0}} D^{+} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{+}"
    if spec_name == "m":
        title = "B^{-} #rightarrow D^{-} D^{0} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{0}"
    if spec_name == "st":
        title = "B^{+} #rightarrow #bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
        d1 = "D^{-}"
        d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
    if spec_name == "s0" or spec_name == "s1":
        title = "B^{0}_{s} #rightarrow D_{s}^{-} D^{+} K^{*0}"
        d1 = "D_{s}^{-}"
        d2 = "D^{+}"

    plotdalitz(title, f"{d1} {d2}", f"{d2} {d3}", hdalitz_12_2k, 2)
    plotdalitz(title, f"{d1} {d2}", f"{d1} {d3}", hdalitz_12_1k, 2)
    plotdalitz(title, f"{d1} {d3}", f"{d2} {d3}", hdalitz_1k_2k, 2)
    # plotdalitz(title, f"{d1} {d3}", "Candidates", hdalitz_1k_p, 1)
    # plotdalitz(title, f"{d2} {d3}", "Candidates", hdalitz_2k_p, 1)
    # plotdalitz(title, f"{d1} {d2}", "Candidates", hdalitz_12_p, 1)

data_path_list = glob.glob(f"{root_basepath}DATA/*.root")
for file in data_path_list:
    make_dalitz(file, "DATA")


mc_path_list = glob.glob(f"{root_basepath}MC/*.root")
for file in mc_path_list:
    make_dalitz(file, "MC")
