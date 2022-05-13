import sys
import os

basedir = os.getcwd().split('sys_studies')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

def grab_count(spec, name):
    print(spec)
    file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/post_d_no_trigger/{spec}_hh.root"
    DecayTree_List = [f"DecayTreeTuple_{name}"]
    tchain = grab_files_and_chain(file_path, DecayTree_List)
    print(tchain)
    rdf_base = RDF(tchain)
    count = rdf_base.Count().GetValue()
    print(count)
    return count

def calc_trigger_mc_TISTOS(spec):

    # trigger_TOS_ALL_L0 = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)"
    # trigger_TIS_ALL_L0 = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"
    #
    # trigger_TOS_HLT1 = "(B_Hlt1TrackMVADecision_TOS == 1 || B_Hlt1TwoTrackMVADecision_TOS == 1)"
    # trigger_TIS_HLT1 = "(B_Hlt1TrackMVADecision_TIS == 1 || B_Hlt1TwoTrackMVADecision_TIS == 1)"
    #
    # trigger_TOS_HLT2 = "(B_Hlt2Topo2BodyDecision_TOS == 1 || B_Hlt2Topo3BodyDecision_TOS == 1 || B_Hlt2Topo4BodyDecision_TOS == 1)"
    # trigger_TIS_HLT2 = "(B_Hlt2Topo2BodyDecision_TIS == 1 || B_Hlt2Topo3BodyDecision_TIS == 1 || B_Hlt2Topo4BodyDecision_TIS == 1)"
    #
    # rdf_base = RDF(tchain)
    # #####################################
    # mc_only_tos_rdf =  rdf_base.Filter(f"({trigger_TOS_ALL_L0})")
    # mc_only_tis_rdf =  rdf_base.Filter(f"({trigger_TIS_ALL_L0})")
    #
    # gen_den_ufloat = ufloat(rdf_base.Count().GetValue(), np.sqrt(rdf_base.Count().GetValue()))
    #
    # mc_nb = ufloat(mc_only_tos_rdf.Count().GetValue(), np.sqrt(mc_only_tos_rdf.Count().GetValue()))
    # mc_nc = ufloat(mc_only_tis_rdf.Count().GetValue(), np.sqrt(mc_only_tis_rdf.Count().GetValue()))
    #
    # e_tos_real = mc_nb / gen_den_ufloat
    # e_tis_real = mc_nc / gen_den_ufloat
    # #################################
    # a_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && {trigger_TOS_ALL_L0}"
    # b_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && !{trigger_TOS_ALL_L0}"
    #
    # c_line = f"({a_line}) && {trigger_TOS_HLT1} && {trigger_TOS_HLT2}"
    # d_line = f"({b_line}) && {trigger_TOS_HLT1} && {trigger_TOS_HLT2}"
    #
    # m_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && {trigger_TOS_HLT1} && {trigger_TOS_HLT2}"
    # n_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && !({trigger_TOS_HLT1} && {trigger_TOS_HLT2})"
    #
    # x_line = f"{trigger_TOS_ALL_L0} && {trigger_TOS_HLT1} && {trigger_TOS_HLT2} && {trigger_TIS_ALL_L0}"
    # y_line = f"{trigger_TOS_ALL_L0} && {trigger_TOS_HLT1} && {trigger_TOS_HLT2} && !{trigger_TIS_ALL_L0}"
    #
    a = grab_count(spec, "a")
    b = grab_count(spec, "b")
    c = grab_count(spec, "c")
    # d = rdf_base.Filter(d_line).Count().GetValue()
    m = grab_count(spec, "m")
    n = grab_count(spec, "n")
    x = grab_count(spec, "x")
    y = grab_count(spec, "y")

    auf = ufloat(a, np.sqrt(a))
    buf = ufloat(b, np.sqrt(b))
    cuf = ufloat(c, np.sqrt(c))
    muf = ufloat(m, np.sqrt(m))
    nuf = ufloat(n, np.sqrt(n))
    xuf = ufloat(x, np.sqrt(x))
    yuf = ufloat(y, np.sqrt(y))

    e_tos = auf/(auf+buf)
    e_hlt_tos_l0_tos = cuf/auf
    e_hlt_tos_l0_tis = muf/(muf+nuf)
    e_tis = xuf/(xuf+yuf)
    e_all = e_tos*e_hlt_tos_l0_tos + e_hlt_tos_l0_tis*e_tis

    dict_trigger = {
        "Spec": spec,
        "e TOS" : f"{e_tos*100.000:.3f}",
        "e TIS" : f"{e_tis*100.000:.3f}",
        "e HLT TOS | L0 TOS" : f"{e_hlt_tos_l0_tos*100.000:.3f}",
        "e HLT TOS | L0 TIS" : f"{e_hlt_tos_l0_tis*100.000:.3f}",
        "e All" : f"{e_all*100.000:.3f}",
        # "e TOS_REAL" : f"{e_tos_real*100.000:.3f}",
        # "e TIS_REAL" : f"{e_tis_real*100.000:.3f}",
        }

    print(dict_trigger)

    eff_df = pd.DataFrame(dict_trigger, index=[0])
    eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/sys_studies/txt_files/{spec}_TISTOS.txt")

def calc_trigger_data_TISTOS(spec, plot_flag, name):

    file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/post_d_no_trigger/{id_to_spec_dict[spec]}_hh.root"
    DecayTree_List = [f"DecayTreeTuple_{name}"]
    tchain = grab_files_and_chain(file_path, DecayTree_List)

    ws = ROOT.RooWorkspace(spec)

    ws.factory(f"B_DTF_M[5230, 5330]")

    fit_var = ws.var("B_DTF_M")
    fit_args = ROOT.RooArgSet(fit_var)
    data_set = ROOT.RooDataSet(f"{spec}_events", f"{spec}_events", tchain, fit_args)

    ws.factory(f"Gaussian::s_fit(B_DTF_M, mean[5280,5270, 5290], width_a[10, 0.1, 30.0])")
    # ws.factory(f"Gaussian::fit_b(B_DTF_M, mean, width_b[15.0, 0.1, 20.0])")
    # ws.factory(f"SUM::s_fit(a_frac[0.5, 0.1, 0.9]*fit_a, fit_b)")
    #
    # ws.factory(f"CBShape::s_fit(B_DTF_M,mean[5280,5270,5290],width[10,0.01,20], alpha[2,0.01,5.0], n[5,0,100])")
    ws.factory(f"Exponential:b_fit(B_DTF_M, c0[0, -1, 1])")
    ws.factory("SUM::a_fit(s_yield[500,0,1000]*s_fit, b_yield[500,0,1000]*b_fit)")

    fit_pdf = ws.pdf("a_fit")


    fit_result = fit_pdf.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())


    if plot_flag:
        frame = fit_var.frame(ROOT.RooFit.Title(f"{spec}"), ROOT.RooFit.Binning(50))
        data_set.plotOn(frame, ROOT.RooFit.Name("data"))
        fit_pdf.plotOn(frame, ROOT.RooFit.Name("pdf"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kBlue))
        fit_pdf.plotOn(frame, ROOT.RooFit.Name("sig"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Components(f"s_fit"))
        fit_pdf.plotOn(frame, ROOT.RooFit.Name("bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Components(f"b_fit"))

        p = residualPlot()
        p.pt.cd()
        xaxis = frame.GetXaxis()
        xaxis.SetTickLength(0)
        xaxis.SetNdivisions(0)
        xaxis.SetLabelSize(0)

        frame.Draw()
        p.pb.cd()

        hpull = frame.pullHist(f"data", f"pdf")
        pull = fit_var.frame()
        pull.addPlotable(hpull, "P")
        pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
        pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize() * p.blabelratio)
        pull.GetYaxis().SetLabelSize(
            ROOT.gStyle.GetLabelSize("Y") * (1 + p.padratio) / (2 * p.padratio)
        )
        pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
        pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
        pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset() / p.blabelratio)
        pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength() * p.blabelratio)
        pull.GetYaxis().SetTitle("Pull")
        pull.GetYaxis().SetNdivisions(202, False)
        pull.GetYaxis().SetRangeUser(-3, 3)
        pull.GetYaxis().CenterTitle()
        pull.Draw("AP")

        # pull.GetXaxis().SetTitle(f" m({data_class.rec_decay_string}) [MeV]")

        save_pdf(p, f"fits", f"{spec}_{name}", rpflag = 1)

    sy = ws.var("s_yield")
    syerr = sy.getPropagatedError(fit_result)
    # syerr = np.sqrt(sy.getValV())
    r = ufloat(sy.getValV(), syerr)
    return(r)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Check trigger efficiency in data vs MC (Control).")
    parser.add_argument('--Spec', choices=["Z_m_p","01_Z_m_p","norm7","norm8"], help = 'Data Spec')
    parser.add_argument('--Plot', action='store_true')

    args = parser.parse_args()

    spec = args.Spec
    plot_flag = args.Plot

    if id_to_spec_dict[spec] == spec:
        type = "data"
        auf = calc_trigger_data_TISTOS(spec, plot_flag, "a")
        buf = calc_trigger_data_TISTOS(spec, plot_flag, "b")
        cuf = calc_trigger_data_TISTOS(spec, plot_flag, "c")
        muf = calc_trigger_data_TISTOS(spec, plot_flag, "m")
        nuf = calc_trigger_data_TISTOS(spec, plot_flag, "n")
        xuf = calc_trigger_data_TISTOS(spec, plot_flag, "x")
        yuf = calc_trigger_data_TISTOS(spec, plot_flag, "y")

        e_tos = auf/(auf+buf)
        e_hlt_tos_l0_tos = cuf/auf
        e_hlt_tos_l0_tis = muf/(muf+nuf)
        e_tis = xuf/(xuf+yuf)
        e_all = e_tos*e_hlt_tos_l0_tos + e_hlt_tos_l0_tis*e_tis

        dict_trigger = {
            "Spec": spec,
            "e TOS" : f"{e_tos*100.000:.3f}",
            "e TIS" : f"{e_tis*100.000:.3f}",
            "e HLT TOS | L0 TOS" : f"{e_hlt_tos_l0_tos*100.000:.3f}",
            "e HLT TOS | L0 TIS" : f"{e_hlt_tos_l0_tis*100.000:.3f}",
            "e All" : f"{e_all*100.000:.3f}",
            }

        print(dict_trigger)

        eff_df = pd.DataFrame(dict_trigger, index=[0])
        eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/sys_studies/txt_files/{spec}_TISTOS.txt")

    else:
        type = "mc"
        id_to_spec_dict
        # file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/*/post_d_no_trigger/{id_to_spec_dict[spec]}.root"
        # DecayTree_List = ["DecayTreeTuple"]
        # tchain = grab_files_and_chain(file_path, DecayTree_List)
        calc_trigger_mc_TISTOS(id_to_spec_dict[spec])
