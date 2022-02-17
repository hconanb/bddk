import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

tla = ROOT.TLatex()
tla.SetTextFont(32)
tla.SetTextColor(1)
tla.SetTextSize(0.03)
tla.SetTextAlign(12)

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")
m_data_file = ROOT.TFile(data_basepath+"m_spectrum.root")
st_data_file = ROOT.TFile(data_basepath+"st_spectrum.root")
s_data_file = ROOT.TFile(data_basepath+"s_spectrum.root")

z_tree = z_data_file.Get("DecayTreeTuple")
zz_tree = zz_data_file.Get("DecayTreeTuple")
p_tree = p_data_file.Get("DecayTreeTuple")
m_tree = m_data_file.Get("DecayTreeTuple")
st_tree = st_data_file.Get("DecayTreeTuple")
s_tree = s_data_file.Get("DecayTreeTuple")

z0_mc_file = ROOT.TFile("/mnt/c/Users/Harris/Google Drive/LHCb/rootfiles/mc/1_z_11198000.root")
zz0_mc_file = ROOT.TFile("/mnt/c/Users/Harris/Google Drive/LHCb/rootfiles/mc/9_zz_11196000.root")
z0_mc_rec_tree = z0_mc_file.Get("DecayTreeTuple")
z0_mc_gen_tree = z0_mc_file.Get("MCDecayTreeTuple")

zz0_mc_rec_tree = zz0_mc_file.Get("DecayTreeTuple")
zz0_mc_gen_tree = zz0_mc_file.Get("MCDecayTreeTuple")

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
#################################################3


peak = 0

z_cut = f"(abs(D1_M - {dpmass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow}) && abs(B_DTF_M - {5279} + {peak}) < 50"
zz_cut = f"(abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow}) && abs(B_DTF_M - {5279} + {peak}) < 50"
p_cut = f"(abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow}) && abs(B_DTF_M - {5279} + {peak}) < 50"
m_cut = f"(abs(D1_M - {dpmass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow}) && abs(B_DTF_M - {5279} + {peak}) < 50"
st_cut = f"(abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow}) && abs(B_DTF_M - {5279} + {peak}) < 50"

z_mc_cut = f"(abs(D1_M - {dpmass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow}) && abs(B_DTF_M - {5279} + {peak}) < 50 && abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313"

s0_cut = f"(abs(D1_M - {dsmass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow}) && abs(B_DTF_M - {5366} + {0}) < 50"
s1_cut = f"(abs(D1_M - {dsmass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow}) && abs(B_DTF_M - {5366} + {145}) < 50"


RDF = ROOT.ROOT.RDataFrame

rdf_Z = RDF(z_tree)
rdf_ZZ = RDF(zz_tree)
rdf_P = RDF(p_tree)
rdf_M = RDF(m_tree)
rdf_ST = RDF(st_tree)
rdf_S = RDF(s_tree)

rdf_z0_rec = RDF(z0_mc_rec_tree)
rdf_zz0_rec = RDF(z0_mc_rec_tree)

rdf_z0_gen = RDF(z0_mc_gen_tree)
rdf_zz0_gen = RDF(z0_mc_gen_tree)

#rdf_list = [rdf_Z, rdf_ZZ, rdf_P, rdf_M, rdf_ST, rdf_S]
rdf_list = [rdf_S, rdf_S]

#rdf_cuts = [z_cut, zz_cut, p_cut, m_cut, st_cut, s_cut]
rdf_cuts = [s0_cut, s1_cut]

#rdf_id = ["z", "zz", "p", "m", "st", "s"]
rdf_id = ["s0", "s1"]

rdf_mc_list = [rdf_z0_rec, rdf_z0_gen]
rdf_mc_cuts = [z_mc_cut, z_mc_cut]
rdf_mc_id = ["z_mc_rec", "z_mc_gen"]
mcf = ["r", "g"]

def plotdalitz(k, title, xname, yname, histo, dim, rflag = 0):

    tltxt = ROOT.TLatex()

    if dim == 2 and rflag == 0:
        histo.GetXaxis().SetTitle(f"m^{{2}}({xname}) [GeV^{{2}}]")
        histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
        histo.GetYaxis().SetTitle(f"m^{{2}}({yname}) [GeV^{{2}}]")
        histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
    if dim == 1 and rflag == 0:
        histo.GetXaxis().SetTitle(f"m({xname}) [GeV]")
        histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000,histo.GetXaxis().GetXmax()/1000)
        bw = histo.GetXaxis().GetBinWidth(5)
        histo.GetYaxis().SetTitle(f"{yname} / ({bw} GeV)")
        histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000,histo.GetYaxis().GetXmax()/1000)
    if dim == 2 and rflag == 1:
        histo.GetXaxis().SetTitle(f"ratio")
        # histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
        histo.GetYaxis().SetTitle(f"ratio")
        # histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
    if dim == 1 and rflag == 1:
        histo.GetXaxis().SetTitle(f"ratio")
        # histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000,histo.GetXaxis().GetXmax()/1000)
        bw = histo.GetXaxis().GetBinWidth(5)
        histo.GetYaxis().SetTitle(f"{yname} / ({bw} GeV)")
        # histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000,histo.GetYaxis().GetXmax()/1000)

    c1 = ROOT.TCanvas("c1","c1")
    c1.SetRightMargin(0.15)
    c1.SetLeftMargin(0.9)
    c1.SetTopMargin(0.1)
    c1.SetBottomMargin(0.9)
    ht1 = histo.DrawCopy("COLZ")

    tltxt.DrawLatexNDC(.4,.925, title)

    pname = histo.GetTitle()
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}_pdf/'):
        os.makedirs(f'plots/{now.month}_{now.day}_pdf/')
    if not os.path.exists(f'plots/{now.month}_{now.day}_png/'):
        os.makedirs(f'plots/{now.month}_{now.day}_png/')
    c1.SaveAs(f"plots/{now.month}_{now.day}_pdf/{k}_{pname}.pdf")
    c1.SaveAs(f"plots/{now.month}_{now.day}_png/{k}_{pname}.png")


rec_list = []
gen_list = []

for i, j, k in zip(rdf_list, rdf_cuts, rdf_id):
#for i, j, k, z in zip(rdf_mc_list, rdf_mc_cuts, rdf_mc_id, mcf):
    # if z == "r":
    rdf_next = i.Define("M2_d1d2", M2_d1d2) \
                .Define("M2_d2kst", M2_d2kst) \
                .Define("M2_d1kst", M2_d1kst) \
                .Define("M_d1d2",M_d1d2) \
                .Define("M_d2kst",M_d2kst) \
                .Define("M_d1kst",M_d1kst) \
                .Filter(j)
    # if z == "g":
    #     rdf_next = i.Define("M2_d1d2", M2_d1d2_mcdtt) \
    #                 .Define("M2_d2kst", M2_d2kst_mcdtt) \
    #                 .Define("M2_d1kst", M2_d1kst_mcdtt) \
    #                 .Define("M_d1d2",M_d1d2_mcdtt) \
    #                 .Define("M_d2kst",M_d2kst_mcdtt) \
    #                 .Define("M_d1kst",M_d1kst_mcdtt)

    dbins = 10

    hdalitz_12_2k = rdf_next.Histo2D(("", f"d1d2_d2kst", dbins, 14.25e6, 21e6, dbins, 7e6, 13e6), "M2_d1d2", "M2_d2kst")
    hdalitz_12_1k = rdf_next.Histo2D(("", f"d1d2_d1kst", dbins, 14.25e6, 21e6, dbins, 7e6, 13e6), "M2_d1d2", "M2_d1kst")
    hdalitz_1k_2k = rdf_next.Histo2D(("", f"d1kst_d2kst", dbins, 7e6, 13e6, dbins, 7e6, 13e6), "M2_d1kst", "M2_d2kst")
    hdalitz_1k_p = rdf_next.Histo1D(("", f"d1kst_p", dbins, 2500, 4000), "M_d1kst")
    hdalitz_2k_p = rdf_next.Histo1D(("", f"d2kst_p", dbins, 2500, 4000), "M_d2kst")
    hdalitz_12_p = rdf_next.Histo1D(("", f"d1d2_p", dbins, 3500, 5500), "M_d1d2")

    # if z == "r":
    # rec_list.append(hdalitz_12_2k)
    # rec_list.append(hdalitz_12_1k)
    # rec_list.append(hdalitz_1k_2k)
    # rec_list.append(hdalitz_1k_p)
    # rec_list.append(hdalitz_2k_p)
    # rec_list.append(hdalitz_12_p)
    # if z == "g":
    #     gen_list.append(hdalitz_12_2k)
        # gen_list.append(hdalitz_12_1k)
        # gen_list.append(hdalitz_1k_2k)
        # gen_list.append(hdalitz_1k_p)
        # gen_list.append(hdalitz_2k_p)
        # gen_list.append(hdalitz_12_p)

    d3 = "K^{*0}"
    if k == "z" or k == "z_mc_rec" or k == "z_mc_gen":
        title = "B^{0} #rightarrow D^{-} D^{+} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{+}"
    if k == "zz":
        title = "B^{0} #rightarrow #bar{D^{0}} D^{0} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{0}"
    if k == "p":
        title = "B^{+} #rightarrow #bar{D^{0}} D^{+} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{+}"
    if k == "m":
        title = "B^{-} #rightarrow D^{-} D^{0} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{0}"
    if k == "st":
        title = "B^{+} #rightarrow #bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
        d1 = "D^{-}"
        d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
    if k == "s0" or k == "s1":
        title = "B^{0}_{s} #rightarrow D_{s}^{-} D^{+} K^{*0}"
        d1 = "D_{s}^{-}"
        d2 = "D^{+}"
    plotdalitz(k, title, f"{d1} {d2}", f"{d2} {d3}", hdalitz_12_2k, 2)
    plotdalitz(k, title, f"{d1} {d2}", f"{d1} {d3}", hdalitz_12_1k, 2)
    plotdalitz(k, title, f"{d1} {d3}", f"{d2} {d3}", hdalitz_1k_2k, 2)
    plotdalitz(k, title, f"{d1} {d3}", "Candidates", hdalitz_1k_p, 1)
    plotdalitz(k, title, f"{d2} {d3}", "Candidates", hdalitz_2k_p, 1)
    plotdalitz(k, title, f"{d1} {d2}", "Candidates", hdalitz_12_p, 1)

# for i in range(0,1):
#     rec_list[i].GetPtr().Divide(gen_list[i].GetPtr())
#     rec_list[i].GetPtr().GetZaxis().SetRangeUser(0,0.002)
#     title = "B^{0} #rightarrow D^{-} D^{+} K^{*0}"
#     d1 = "D^{-}"
#     d2 = "D^{+}"
#
#     plotdalitz("mcr", title, f"{d1} {d2}", f"{d2} {d3}", rec_list[i], 2, 1)
    # plotdalitz("mcr", title, f"{d1} {d2}", f"{d1} {d3}", hdalitz_12_1k, 2, 1)
    # plotdalitz("mcr", title, f"{d1} {d3}", f"{d2} {d3}", hdalitz_1k_2k, 2, 1)
    # plotdalitz("mcr", title, f"{d1} {d3}", "Candidates", hdalitz_1k_p, 1, 1)
    # plotdalitz("mcr", title, f"{d2} {d3}", "Candidates", hdalitz_2k_p, 1, 1)
    # plotdalitz("mcr", title, f"{d1} {d2}", "Candidates", hdalitz_12_p, 1, 1)
