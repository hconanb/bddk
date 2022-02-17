import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *


norm8_mc_file = ROOT.TFile(mc_basepath+"norm8_norm8_11198030.root")
norm8_mc_tree = norm8_mc_file.Get("DecayTreeTuple")

norm8_data_file = ROOT.TFile(data_basepath+"norm8_data.root")
norm8_data_tree = norm8_data_file.Get("DecayTreeTuple")

norm7_mc_file = ROOT.TFile(mc_basepath+"norm7_norm7_12197009.root")
norm7_mc_tree = norm7_mc_file.Get("DecayTreeTuple")

norm7_data_file = ROOT.TFile(data_basepath+"norm7_data.root")
norm7_data_tree = norm8_data_file.Get("DecayTreeTuple")

name = "ws_norm"

ws = ROOT.RooWorkspace(name)
ws.factory("B_DTF_M[-1000,7000]")
ws.factory("D1_M[0,6000]")
ws.factory("D2_M[0,6000]")

ws.factory("Gaussian::norm8_mc_a(B_DTF_M, mean_norm8[5279, 5270, 5290], width_norm8_a[9, 0.0, 10])")
ws.factory("Gaussian::norm8_mc_b(B_DTF_M, mean_norm8, width_norm8_b[1, 0.0, 10])")
ws.factory("SUM::norm8_signal(norm8_b_frac[0.5,0,1]*norm8_mc_a, norm8_mc_b)")

ws.factory("Gaussian::norm7_mc_a(B_DTF_M, mean_norm7[5279, 5270, 5290], width_norm7_a[9, 0.0, 10])")
ws.factory("Gaussian::norm7_mc_b(B_DTF_M, mean_norm7, width_norm7_b[1, 0.0, 10])")
ws.factory("SUM::norm7_signal(norm7_b_frac[0.5,0,1]*norm7_mc_a, norm7_mc_b)")

d1_mass = ws.var("D1_M")
d2_mass = ws.var("D2_M")
b_dtf_m = ws.var("B_DTF_M")

norm8_cut = f"abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow}"
norm7_cut = f"abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow}"

m_args = ROOT.RooArgSet(d1_mass, d2_mass, b_dtf_m)
emptyargs = ROOT.RooArgSet()

eff_norm8 = eff_calc(norm8_mc_file, d0mass, dpmass)
ws.factory(f"eff_norm8[{eff_norm8}]")
eff_norm7 = eff_calc(norm7_mc_file, d0mass, d0mass)
ws.factory(f"eff_norm7[{eff_norm7}]")

e11 = get_mc_error(norm8_mc_file, d0mass, dpmass)
e22 = get_mc_error(norm7_mc_file, d0mass, d0mass)

print(e11, e22)

ws.Print()

# ws.factory("Chebychev:norm8_bkg(B_DTF_M,{c0[0.,-2,20],c1[0.,-2,2],c2[0.,-2,2]})")
ws.factory("Exponential:norm8_bkg(B_DTF_M, c0_n8[0, -2, 2])")
ws.factory("SUM::norm8_fit(n_norm8[100,0,10000]*norm8_signal,n_norm8_bkg[100,0,10000]*norm8_bkg)")

# ws.factory("Chebychev:norm7_bkg(B_DTF_M,{c3[0.,-2,20],c4[0.,-2,2],c5[0.,-2,2]})")
ws.factory("Exponential:norm7_bkg(B_DTF_M, c0_n7[0, -2, 2])")
ws.factory("SUM::norm7_fit(n_norm7[100,0,10000]*norm7_signal,n_norm7_bkg[100,0,10000]*norm7_bkg)")
#
n8sig = ws.var("n_norm8")
n8bkg = ws.var("n_norm8_bkg")
#
n7sig = ws.var("n_norm7")
n7bkg = ws.var("n_norm7_bkg")
#
params8 = ROOT.RooArgSet(n8sig, n8bkg, "p8")
params7 = ROOT.RooArgSet(n7sig, n7bkg, "p7")
#
fit_plots(ws, "norm8_fit", norm8_data_tree, m_args, params8, norm8_cut, b_dtf_m,  5230, 5330, 5230, 5330, 1)
fit_plots(ws, "norm7_fit", norm7_data_tree, m_args, params7, norm7_cut, b_dtf_m,  5230, 5330, 5230, 5330, 1)

ws.Print()

ws.writeToFile(f"../workspaces/{name}"+".root")
print("Saved to:", f"../workspaces/{name}"+".root")
