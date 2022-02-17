import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *


name = "mc_ws"

m_1z_file = ROOT.TFile(basepath+"mc/1_z_11198000.root")
m_1z_tree = m_1z_file.Get("DecayTreeTuple")

m_2az_file = ROOT.TFile(basepath+"mc/2a_z_11198400.root")
m_2az_tree = m_2az_file.Get("DecayTreeTuple")

m_4az_file = ROOT.TFile(basepath+"mc/4a_z_11198401.root")
m_4az_tree = m_4az_file.Get("DecayTreeTuple")

m_7zz_file = ROOT.TFile(basepath+"mc/7b_zz_12197022.root")
m_7zz_tree = m_7zz_file.Get("DecayTreeTuple")

m_9zz_file = ROOT.TFile(basepath+"mc/9_zz_11196000.root")
m_9zz_tree = m_9zz_file.Get("DecayTreeTuple")

m_10gzz_file = ROOT.TFile(basepath+"mc/10_g_zz_11196200.root", "Update")
m_10gzz_tree = m_10gzz_file.Get("DecayTreeTuple")

m_10pizz_file = ROOT.TFile(basepath+"mc/10_pi_zz_11196410.root", "Update")
m_10pizz_tree = m_10pizz_file.Get("DecayTreeTuple")

m_12gzz_file = ROOT.TFile(basepath+"mc/12_g_zz_11196210.root", "Update")
m_12gzz_tree = m_12gzz_file.Get("DecayTreeTuple")

m_12pizz_file = ROOT.TFile(basepath+"mc/12_pi_zz_11196420.root", "Update")
m_12pizz_tree = m_12pizz_file.Get("DecayTreeTuple")

w1  = array( 'd', [ 0.343 ] )
w1Branch_10 = m_10gzz_tree.Branch("w1", w1, "w1/D")
w1Branch_12 = m_12gzz_tree.Branch("w1", w1, "w1/D")

w2  = array( 'd', [ 0.647 ] )
w2Branch_10 = m_10pizz_tree.Branch("w2", w2, "w2/D")
w2Branch_12 = m_12pizz_tree.Branch("w2", w2, "w2/D")

for i in range(m_10gzz_tree.GetEntries()):
    m_10gzz_tree.Fill()

for j in range(m_10pizz_tree.GetEntries()):
    m_10pizz_tree.Fill()

for i in range(m_12gzz_tree.GetEntries()):
    m_12gzz_tree.Fill()

for j in range(m_12pizz_tree.GetEntries()):
    m_12pizz_tree.Fill()

# 10_g_zz_11196200.root   12_pig_zz_11196620.root  16_s_13198600.root  4a_z_11198401.root   8b_pi_zz_12197420.root
# 10_pi_zz_11196410.root  13_s_13198040.root       1_z_11198000.root   7a_p_12197400.root   9_zz_11196000.root
# 12_g_zz_11196210.root   14_s_13198200.root       2a_z_11198400.root  7b_zz_12197022.root
# 12_pi_zz_11196420.root  15_s_13198400.root       2b_p_11198005.root  8a_p_12197401.root

ws = ROOT.RooWorkspace(name)
ws.factory("B_DTF_M[4800,5600]")
ws.factory("D1_M[0,6000]")
ws.factory("D2_M[0,6000]")
ws.factory("w1[0,1]")
ws.factory("w2[0,1]")
ws.factory("B_TRUEID[-1000,1000]")
ws.factory("KST_TRUEID[-1000,1000]")
ws.factory("TRUE_BM_REC[0,6500]")
ws.factory("Res[-100,100]")

#fit for z1,z2,z3 bdtfm
ws.factory("Gaussian::1az(B_DTF_M, mean_1az[5285,5270,5290], width_1az[7.0, 0, 30.0])")
ws.factory("Gaussian::1bz(B_DTF_M, mean_1az, width_1bz[5.0, 0, 30.0])")
ws.factory("SUM::z0_fit(frac_1za[0.5, 0.0, 1.0]*1az, 1bz)")

ws.factory(f"BifurGauss::z1_fit(B_DTF_M,mean_z1[5130,5100,5160],width_2az_l[10,0.01,1000],width_2az_r[10,0.01,1000])")
ws.factory(f"BifurGauss::z2_fit(B_DTF_M,mean_z2[4985,4960,5000],width_z2_l[10,0.01,1000],width_z2_r[10,0.01,1000])")


# ws.factory("Gaussian::9azz(B_DTF_M, mean_9azz[5285,5270,5290], width_9azz[7.0, 0, 30.0])")
# ws.factory("Gaussian::9bzz(B_DTF_M, mean_9azz, width_9bzz[5.0, 0, 30.0])")
# ws.factory("SUM::9zz_fit(frac_9zza[0.5, 0.0, 1.0]*9azz, 9bzz)")
#
# ws.factory("Gaussian::10gazz(B_DTF_M, mean_10gazz[5100,5050,5260], width_10gazz[5.0, 0, 30.0])")
# ws.factory("Gaussian::10gbzz(B_DTF_M, mean_10gazz, width_10gbzz[5.0, 0, 30.0])")
# ws.factory("SUM::10gzz_fit(frac_10gzza[0.5, 0.0, 1.0]*10gazz, 10gbzz)")
#
# ws.factory("BifurGauss::10pizz_fit(B_DTF_M,mean_10pizz[5100,4900,5260],width_10pizz_l[5,0.01,1000],width_10pizz_r[10,0.01,1000])")
# ws.factory("BifurGauss::10gzz_fit(B_DTF_M,mean_10gzz[5100,4900,5260],width_10gzz_l[5,0.01,1000],width_10gzz_r[10,0.01,1000])")
# ws.factory("Landau::10gzz_fit(B_DTF_M,mean_10gzz[5100,4900,5260],width_10gzz_l[5,0.01,1000])")

d1_mass = ws.var("D1_M")
d2_mass = ws.var("D2_M")
b_trueid = ws.var("B_TRUEID")
kst_trueid = ws.var("KST_TRUEID")
b_dtf_m = ws.var("B_DTF_M")
w1 = ws.var("w1")
w2 = ws.var("w2")
# tbr = ws.var("TRUE_BM_REC")
# res = ws.var("Res")

z_mcut = f"(abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313 && abs(D1_M - {dpmass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow})"
zz_mcut = f"(abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313 && abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow})"
# p_mcut = f"(abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313 && abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow})"

m_args = ROOT.RooArgSet(d1_mass, d2_mass, b_trueid, kst_trueid, w1, w2, b_dtf_m)
# tbr_args = ROOT.RooArgSet(d1_mass, d2_mass, b_trueid, kst_trueid, tbr)
# res_args = ROOT.RooArgSet(d1_mass, d2_mass, b_trueid, kst_trueid, res)
emptyargs = ROOT.RooArgSet()

# newds = add_mc([m_1z_tree, m_2az_tree, m_4az_tree], m_args, z_mcut)
# getattr(ws, 'import')(newds)

newds = add_mc_2(m_10gzz_tree, m_10pizz_tree, m_args, zz_mcut)
display_data("test", newds, b_dtf_m, 5000, 5250)

# #fit the b_dtm_plots
# fit_plots(ws, "z0_fit", m_1z_tree, m_args, emptyargs, z_mcut, b_dtf_m, 5250, 5310, 5250, 5310, 1)
# fit_plots(ws, "z1_fit", m_2az_tree, m_args, emptyargs, z_mcut, b_dtf_m, 5090, 5170, 5090, 5170, 1)
# fit_plots(ws, "z2_fit", m_4az_tree, m_args, emptyargs, z_mcut, b_dtf_m, 4920, 5050, 4920, 5050, 1)
# # # fit_plots(ws, "9zz_fit", m_9zz_tree, m_args, emptyargs, zz_mcut, b_dtf_m, 5250, 5310, 5250, 5310, 0)
# # # fit_plots(ws, "10gzz_fit", m_10gzz_tree, m_args, emptyargs, zz_mcut, b_dtf_m, 5000, 5250, 5000, 5250, 0)
# # # fit_plots(ws, "10pizz_fit", m_10pizz_tree, m_args, emptyargs, zz_mcut, b_dtf_m, 5000, 5250, 5000, 5250, 0)
# #
# # #truth fits for second fit strategy
# data_set_11198400 = ROOT.RooDataSet("11198400_ds","11198400_ds", m_2az_tree, tbr_args, z_mcut)
# getattr(ws, 'import')(data_set_11198400)
# ws.factory("KeysPdf::2az_t(TRUE_BM_REC,11198400_ds,MirrorBoth,2)")
# ws.factory("Gaussian::4az_t(TRUE_BM_REC, mean_t[4990,4900,5200], width_t[5.0,0.01,100])")
#
# fit_plots(ws, "2az_t", m_2az_tree, tbr_args, emptyargs, z_mcut, tbr, 5100, rmax1, 5100, rmax1, 1)
# # fit_plots(ws, "3az_t", m_2az_tree, tbr_args, emptyargs, z_mcut, tbr, 5100, rmax1, 5100, rmax1, 1)
# fit_plots(ws, "4az_t", m_4az_tree, tbr_args, emptyargs, z_mcut, tbr, 4945, 5025, 4945, 5025, 1)
#
#
# ws.writeToFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/workspaces/{name}"+".root")
# print("Saved to:", f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/workspaces/{name}"+".root")
