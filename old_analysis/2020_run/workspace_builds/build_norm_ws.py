import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

dwindow_file = ROOT.TFile(f"dmass/d_mass_fits.root","READ")
dwindow_ws = dwindow_file.Get("d_mass_fits")

def build_norm_ws(name, trigger_flag, data_tree):

    nws = ROOT.RooWorkspace(f"{name}_{trigger_flag}")
    nws.factory("B_DTF_M[-1000,7000]")
    nws.factory("D1_M[0,6000]")
    nws.factory("D2_M[0,6000]")
    nws.factory("D1_DIRA_ORIVX[-5,5]")
    nws.factory("D2_DIRA_ORIVX[-5,5]")
    nws.factory("D1_FDCHI2_ORIVX[0,1000]")
    nws.factory("D2_FDCHI2_ORIVX[0,1000]")

    # nws.factory("B_L0Global_TOS[0,1]")
    # nws.factory("B_L0HadronDecision_TIS[0,1]")
    # nws.factory("B_L0MuonDecision_TIS[0,1]")
    # nws.factory("B_L0ElectronDecision_TIS[0,1]")
    # nws.factory("B_L0PhotonDecision_TIS[0,1]")
    # nws.factory("Gaussian::norm_a(B_DTF_M, mean_norm[5279, 5270, 5290], width_norm_a[2, 0.0, 10.0])")
    # nws.factory("Gaussian::norm_b(B_DTF_M, mean_norm, width_norm_b[8, 0.0, 10.0])")
    # nws.factory("SUM::norm_signal(norm_a_frac[0.5,0,1]*norm_a, norm_b)")

    d1_mass = nws.var("D1_M")
    d2_mass = nws.var("D2_M")
    b_dtf_m = nws.var("B_DTF_M")
    d1_dira = nws.var("D1_DIRA_ORIVX")
    d2_dira = nws.var("D2_DIRA_ORIVX")

    b_dtf_m.setRange("myrange", 5230, 5330)

    m_cut = get_dwindow_values(dwindow_ws, name)
    norm_cut = f"{m_cut}"
    #&& {trigger_cuts[trigger_flag]}"

    m_args = ROOT.RooArgSet(d1_mass, d2_mass, b_dtf_m, d1_dira, d2_dira)

    norm_data_set = ROOT.RooDataSet(f"{name}_{trigger_flag}",f"{name}_{trigger_flag}", data_tree, m_args, norm_cut)

    nws.factory("Gaussian::norm_signal(B_DTF_M, mean_norm[5279, 5270, 5290], width_norm_a[2, 0.0, 10.0])")
    nws.factory("Exponential:norm_bkg(B_DTF_M, c0_n[0, -5, 5])")
    nws.factory("SUM::norm_fit(n_norm_signal[100,0,100000]*norm_signal,n_norm_bkg[100,0,100000]*norm_bkg)")

    norm_fit = nws.pdf("norm_fit")
    norm_fit.fitTo(norm_data_set, ROOT.RooFit.Range("myrange"), ROOT.RooFit.PrintLevel(0))

    nws.Print()
    getattr(nws,'import')(norm_data_set)
    nws.writeToFile(f"normalization/{name}_{trigger_flag}.root")
    print(f"wrote normalization/{name}_{trigger_flag}.root", "\n")

norm7_tot_data_file = ROOT.TFile(data_basepath+"norm7_data.root")
norm7_tot_data_tree = norm7_tot_data_file.Get("DecayTreeTuple")

norm7_tt_data_file = ROOT.TFile(TT_data_basepath + "TT_norm7_data.root")
norm7_t_data_tree = norm7_tt_data_file.Get("DecayTreeTuple_T")
norm7_ntat_data_tree = norm7_tt_data_file.Get("DecayTreeTuple_nTaT")

norm8_tot_data_file = ROOT.TFile(data_basepath+"norm8_data.root")
norm8_tot_data_tree = norm8_tot_data_file.Get("DecayTreeTuple")

norm8_tt_data_file = ROOT.TFile(TT_data_basepath + "TT_norm8_data.root")
norm8_t_data_tree = norm8_tt_data_file.Get("DecayTreeTuple_T")
norm8_ntat_data_tree = norm8_tt_data_file.Get("DecayTreeTuple_nTaT")

build_norm_ws("norm8", "ToT", norm8_tot_data_tree)
build_norm_ws("norm8", "T", norm8_t_data_tree)
build_norm_ws("norm8", "nTaT", norm8_ntat_data_tree)

build_norm_ws("norm7", "ToT", norm7_tot_data_tree)
build_norm_ws("norm7", "T", norm7_t_data_tree)
build_norm_ws("norm7", "nTaT", norm7_ntat_data_tree)
