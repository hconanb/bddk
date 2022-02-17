import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

def build_data_fits(dws, bkg_flag, peak_0_flag, peak_1_flag, peak_2_flag, nflag):

    # dws.factory(f"Exponential:s_spectrum_bkg(B_DTF_M, c0_s[0, -2, 2])")
    dws.factory(f"Bernstein:s_spectrum_bkg(B_DTF_M,{{c0_s[1,0,10], c1_s[1,0,10], c2_s[1,0,10], c3_s[1,0,10]}})")
    dws.factory(f"Gaussian::s0_a(B_DTF_M, mean_s0[5366,5350,5375], width_a_s0[8.0,0.01,20])")
    dws.factory(f"Gaussian::s0_b(B_DTF_M, mean_s0, width_b_s0[7.0,0.01,20])")
    dws.factory(f"SUM::s0_fit(s0_a_frac[0.5,0,1]*s0_a, s0_b)")
    # dws.factory(f"BifurGauss::s1_fit(B_DTF_M,mean_s1[5230,5200,5250],width_1_sl[15,0.01,100],width_1_sr[10,0.01,100])")
    # dws.factory(f"BifurGauss::s2_fit(B_DTF_M,mean_s2[5090,5050,5130],width_2_sl[4,0.01,100],width_2_sr[5,0.01,100])")
    dws.factory(f"CBShape::s1_fit(B_DTF_M,mean_s1[5230,5200,5250],width_1_s[10,0,50],alpha_1_s[1,0.1,100], n1_s[3,1,100])")
    dws.factory(f"CBShape::s2_fit(B_DTF_M,mean_s2[5075,5050,5100],width_2_s[10,0,50],alpha_5_s[1,0.1,100], n5_s[3,1,100])")

    speclist = [
    "z",
    "zz",
    "p",
    "m",
    "st",
    ]

    for spec in speclist:

        if bkg_flag == "Exponential":
            dws.factory(f"Exponential:{spec}_spectrum_bkg(B_DTF_M, c0_{spec}[0, -2, 2])")
        if bkg_flag == "Chebychev":
            dws.factory(f"Chebychev:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[0.,-3,3],c1_{spec}[0.,-3,3]}})")
        if bkg_flag == "Bernstein":
            dws.factory(f"Bernstein:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[1,0,10], c1_{spec}[1,0,10], c2_{spec}[1,0,10], c3_{spec}[1,0,10]}})")

        if peak_0_flag == "DG":
            dws.factory(f"Gaussian::{spec}0_a(B_DTF_M, mean_{spec}0[5279,5270,5290], width_a_{spec}0[8.0,0.01,20])")
            dws.factory(f"Gaussian::{spec}0_b(B_DTF_M, mean_{spec}0, width_b_{spec}0[7.0,0.01,20])")
            dws.factory(f"SUM::{spec}0_fit({spec}0_a_frac[0.5,0,1]*{spec}0_a, {spec}0_b)")
        if peak_0_flag == "G":
            dws.factory(f"Gaussian::{spec}0_fit(B_DTF_M, mean_{spec}0[9,-20,20], width_a_{spec}0[400,200,1000])")

        if peak_1_flag == "GEP":
            dws.factory(f"RooGaussExp::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}l[15,0.01,100],alphae_{spec}[10,0.01,100])")
        if peak_1_flag == "BGEP":
            dws.factory(f"RooBifurGaussExp::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}l[15,0.01,100],width_1_{spec}R[15,0.01,100],alphae_{spec}[10,0.01,100],alphae2_{spec}[10,0.01,100])")
        if peak_1_flag == "BG":
            dws.factory(f"BifurGauss::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}l[15,0.01,100],width_1_{spec}r[10,0.01,100])")
        if peak_1_flag == "cb2":
            dws.factory(f"CBShape::{spec}1_a_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}[10,0,50],alpha_1_{spec}[-5,-100,-0.1], n1_{spec}[3,1,100])")
            dws.factory(f"CBShape::{spec}1_b_fit(B_DTF_M,mean_{spec}1,                width_1_{spec},         alpha_2_{spec}[11,0.01,100], n2_{spec}[6,1,100])")
            dws.factory(f"SUM::{spec}1_fit({spec}1_a_frac[0.5,0,1]*{spec}1_a_fit, {spec}1_b_fit)")
        if peak_1_flag == "cb1R":
            dws.factory(f"CBShape::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}[10,0,50],alpha_1_{spec}[-5,-100,-0.1], n1_{spec}[3,1,100])")
        if peak_1_flag == "cb1L":
            dws.factory(f"CBShape::{spec}1_fit(B_DTF_M,mean_{spec}1[5130,5100,5150],width_1_{spec}[10,0,50],alpha_1_{spec}[1,0.1,100], n1_{spec}[3,1,100])")

        if peak_2_flag == "GEP":
            dws.factory(f"RooGaussExp::{spec}2_fit(B_DTF_M,mean_{spec}2[4999,4950,5020],width_2_{spec}l[15,0.01,100],alpha2e_{spec}[10,0.01,100])")
        if peak_2_flag == "BGEP":
            dws.factory(f"RooBifurGaussExp::{spec}2_fit(B_DTF_M,mean_{spec}2[4999,4950,5020],width_2_{spec}l[15,0.01,100],width_2_{spec}R[15,0.01,100],alpha2e_{spec}[10,0.01,100],alpha2e2_{spec}[10,0.01,100])")
        if peak_2_flag == "BG":
            dws.factory(f"BifurGauss::{spec}2_fit(B_DTF_M,mean_{spec}2[4999,4950,5020],width_2_{spec}l[4,0.01,100],width_2_{spec}r[5,0.01,100])")
        if peak_2_flag == "cb2":
            dws.factory(f"CBShape::{spec}2_a_fit(B_DTF_M,mean_{spec}2[4995,4950,5020],width_2_{spec}[10,0,50],alpha_3_{spec}[-5,-100,-0.1], n3_{spec}[3,1,100])")
            dws.factory(f"CBShape::{spec}2_b_fit(B_DTF_M,mean_{spec}2,                width_2_{spec},         alpha_4_{spec}[11,0.01,100], n4_{spec}[6,1,100])")
            dws.factory(f"SUM::{spec}2_fit({spec}2_a_frac[0.5,0,1]*{spec}2_a_fit, {spec}2_b_fit)")
        if peak_2_flag == "cb1R":
            dws.factory(f"CBShape::{spec}2_fit(B_DTF_M,mean_{spec}2[4995,4950,5020],width_2_{spec}[10,0,50],alpha_5_{spec}[-5,-100,-0.1], n5_{spec}[3,1,100])")
        if peak_2_flag == "cb1L":
            dws.factory(f"CBShape::{spec}2_fit(B_DTF_M,mean_{spec}2[4995,4950,5020],width_2_{spec}[10,0,50],alpha_5_{spec}[1,0.1,100], n5_{spec}[3,1,100])")

    if nflag == "calc":
        dws.factory("SUM::z_spectrum_all_fit(n_01_z*z0_fit,n_23_z*z1_fit,n_04_z*z2_fit,n_z_bkg[100,0,10000]*z_spectrum_bkg)")
        dws.factory("SUM::zz_spectrum_all_fit(n_09_zz*zz0_fit,n_71011_zz*zz1_fit,n_4812_zz*zz2_fit,n_zz_bkg[100,0,100000]*zz_spectrum_bkg)")
        dws.factory("SUM::p_spectrum_all_fit(n_05_p*p0_fit,n_267_p*p1_fit,n_48_p*p2_fit,n_p_bkg[100,0,10000]*p_spectrum_bkg)")
        dws.factory("SUM::m_spectrum_all_fit(n_03_m*m1_fit,n_04_m*m2_fit,n_m_bkg[100,0,10000]*m_spectrum_bkg)")
        dws.factory("SUM::st_spectrum_all_fit(n_07_st*st0_fit,n_48_st*st1_fit,n_st_bkg[100,0,1000]*st_spectrum_bkg)")
        dws.factory("SUM::s_spectrum_all_fit(n_13_s*s0_fit,n_1415_s*s1_fit,n_16_s*s2_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")

    dws.Print()
    return dws
def fit_ws(ttflag):

    dws_base_file = ROOT.TFile(f"signal/dws_{ttflag}.root")
    dws = dws_base_file.Get(f"dws_{ttflag}")

    if gc_onflag == 1:
        clist = dws.allPdfs()

    all_cats = dws.cat("all_cats")
    b_dtf_m = dws.var("B_DTF_M")
    b_dtf_m.setRange("myrange", bmin, bmax)

    all_data_sets = dws.data("all_data_sets")
    # fd_cut = 2/10*(i - 1)
    # test_fd_cut = f"D1_FDCHI2_ORIVX > {fd_cut} && D2_FDCHI2_ORIVX > {fd_cut}"

    # if i != 0:
    #     all_data_sets = tds.reduce(ROOT.RooArgSet(all_cats,b_dtf_m), f"{dira_cut} && {test_fd_cut}")
    # else:
    #     all_data_sets = tds

    dws = build_data_fits(dws, fbkg, f0, f1, f2, nflag)

    z_model = dws.pdf("z_spectrum_all_fit")
    zz_model = dws.pdf("zz_spectrum_all_fit")
    p_model = dws.pdf("p_spectrum_all_fit")
    m_model = dws.pdf("m_spectrum_all_fit")
    st_model = dws.pdf("st_spectrum_all_fit")
    # s_model = dws.pdf("s_spectrum_all_fit")

    all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)
    all_fit.addPdf(z_model, "z_spectrum")
    all_fit.addPdf(zz_model, "zz_spectrum")
    all_fit.addPdf(p_model, "p_spectrum")
    all_fit.addPdf(m_model, "m_spectrum")
    all_fit.addPdf(st_model, "st_spectrum")
    # all_fit.addPdf(s_model, "s_spectrum")

    if gc_onflag == 0:
        fr = all_fit.fitTo(all_data_sets, ROOT.RooFit.Range("myrange"), ROOT.RooFit.PrintLevel(0))
        # fit_file = ROOT.TFile(f"fits/fr_{ttflag}.root","RECREATE")
        # fr.Write(f"rf_{ttflag}")
        # fit_file.Close()
    if gc_onflag == 1:
        fr = all_fit.fitTo(all_data_sets, ROOT.RooFit.Range("myrange"), ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(clist))
        # fit_file = ROOT.TFile(f"fits/fr_{ttflag}.root","RECREATE")
        # fr.Write(f"rf_{ttflag}")
        # fit_file.Close()

    fit_dws = ROOT.RooWorkspace(f"fit_{ttflag}")

    getattr(fit_dws,'import')(all_data_sets)
    getattr(fit_dws,'import')(all_fit)

    fit_dws.writeToFile(f"fits/fit_{ttflag}.root")


fit_ws("ToT")
# fit_ws("T")
# fit_ws("nTaT")
