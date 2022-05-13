import sys
import os

basedir = os.getcwd().split('ntuple_building')[0]
sys.path.append(basedir)

from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

def apply_dwindow_cuts(spec, inputfile, outputfile, type, snap_flag, txt_flag):

        base_file = ROOT.TFile(inputfile, "READ")
        base_tree = base_file.Get(f"DecayTreeTuple")
        rdf_base = RDF(base_tree)

        rdf_base = rdf_numbers(rdf_base, "No D")

        sig_sig_line = get_dwindow_values(spec, "sig")
        sb2_line = get_dwindow_values(spec, "sb")

        print(sig_sig_line)

        if "signal_and_sb" in outputfile:

            rdf_ToT_dsig = rdf_base.apply_filter(sig_sig_line, f"d_sig_for_{spec}")
            rdf_ToT_dsb = rdf_base.apply_filter(sb2_line, f"d_sb_for_{spec}")

            if snap_flag:
                if os.path.exists(outputfile):
                    os.remove(outputfile)
                    print("First deleting old ", outputfile)
                else:
                    print("making ", outputfile)

                clist = rdf_base.rdf.GetColumnNames()
                if type == "mc":
                    clist_f = ["B_dtf_c_nPV","nPV", "B_dtf_nPV"]
                if type == "data":
                    clist_f = ["nPV", "B_dtf_nPV"]
                for name in clist:
                    if name not in clist_f:
                        clist_f.append(name)

                print(f"Starting snapshot for {outputfile}")

                rdf_ToT_dsig.rdf.Snapshot(f"DecayTreeTuple_SIG", outputfile, clist_f, opts)
                rdf_ToT_dsb.rdf.Snapshot(f"DecayTreeTuple_SB", outputfile, clist_f, opts)

                print(f"finished snapshot for {outputfile}")

            if txt_flag:

                    eff_dict_dwin = {
                        "Scheme ID" : id_to_scheme_dict[spec],
                        "Year": year,
                        "$\epsilon_{D} Sig$": f"{rdf_ToT_dsig.calc_can_eff()*100.000:.3f}",
                        "$\epsilon_{D} Sb$": f"{rdf_ToT_dsb.calc_can_eff()*100.000:.3f}",
                    }
                    eff_df = pd.DataFrame(eff_dict_dwin , index=[0])
                    eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/d_win/can_{spec}_{year}.txt")

                    eff_dict_dwin = {
                        "Scheme ID" : id_to_scheme_dict[spec],
                        "Year": year,
                        "$\epsilon_{D} Sig$": f"{rdf_ToT_dsig.calc_event_eff()*100.000:.3f}",
                        "$\epsilon_{D} Sb$": f"{rdf_ToT_dsb.calc_event_eff()*100.000:.3f}",
                    }
                    eff_df = pd.DataFrame(eff_dict_dwin , index=[0])
                    eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/d_win/event_{spec}_{year}.txt")


        if "d_no_trigger" in outputfile:

            trigger_TOS_ALL_L0 = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)"
            # trigger_TOS_H_L0 = "(B_L0HadronDecision_TOS == 1)"
            trigger_TIS_ALL_L0 = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"
            trigger_TOS_HLT1 = "(B_Hlt1TrackMVADecision_TOS == 1 || B_Hlt1TwoTrackMVADecision_TOS == 1)"
            # trigger_TIS_HLT1 = "(B_Hlt1TrackMVADecision_TIS == 1 || B_Hlt1TwoTrackMVADecision_TIS == 1)"
            trigger_TIS_HLT1 = "(B_Hlt1Global_TIS == 1)"
            trigger_TOS_HLT2 = "(B_Hlt2Topo2BodyDecision_TOS == 1 || B_Hlt2Topo3BodyDecision_TOS == 1 || B_Hlt2Topo4BodyDecision_TOS == 1)"
            # trigger_TIS_HLT2 = "(B_Hlt2Topo2BodyDecision_TIS == 1 || B_Hlt2Topo3BodyDecision_TIS == 1 || B_Hlt2Topo4BodyDecision_TIS == 1)"
            trigger_TIS_HLT2 = "(B_Hlt2Global_TIS == 1)"

            alltis_l0tos_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && {trigger_TOS_ALL_L0}"
            alltis_nl0tos_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && !{trigger_TOS_ALL_L0}"

            nalltis_l0tos_line =  f"!({trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2}) && {trigger_TOS_ALL_L0}"
            nalltis_nl0tos_line = f"!({trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2}) && !{trigger_TOS_ALL_L0}"

            c_line = f"({a_line}) && {trigger_TOS_HLT1} && {trigger_TOS_HLT2}"
            d_line = f"({b_line}) && {trigger_TOS_HLT1} && {trigger_TOS_HLT2}"

            m_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && {trigger_TOS_HLT1} && {trigger_TOS_HLT2}"
            n_line = f"{trigger_TIS_ALL_L0} && {trigger_TIS_HLT1} && {trigger_TIS_HLT2} && !({trigger_TOS_HLT1} && {trigger_TOS_HLT2})"

            x_line = f"{trigger_TOS_ALL_L0} && {trigger_TOS_HLT1} && {trigger_TOS_HLT2} && {trigger_TIS_ALL_L0}"
            y_line = f"{trigger_TOS_ALL_L0} && {trigger_TOS_HLT1} && {trigger_TOS_HLT2} && !{trigger_TIS_ALL_L0}"

            rdf_ToT_dsig = rdf_base.Filter(sig_sig_line, f"d_sig_for_{spec}")
            rdf_ToT_dsig = rdf_ToT_dsig.Filter("B_DTF_M <= 5330 && B_DTF_M >= 5230")

            a = rdf_ToT_dsig.Filter(alltis_l0tos_line)
            b = rdf_ToT_dsig.Filter(alltis_nl0tos_line)
            c = rdf_ToT_dsig.Filter(c_line)
            d = rdf_ToT_dsig.Filter(d_line)
            m = rdf_ToT_dsig.Filter(m_line)
            n = rdf_ToT_dsig.Filter(n_line)
            x = rdf_ToT_dsig.Filter(x_line)
            y = rdf_ToT_dsig.Filter(y_line)

            # ap = rdf_ToT_dsig.Filter(ap_line)
            # bp = rdf_ToT_dsig.Filter(bp_line)
            # print(f"base: {rdf_ToT_dsig.Count().GetValue()}")
            # print(f"a: {a.Count().GetValue()}")
            # print(f"b: {b.Count().GetValue()}")
            # print(f"ap: {ap.Count().GetValue()}")
            # print(f"bp: {bp.Count().GetValue()}")

            print(f"Starting snapshotfor {outputfile}")
            clist = ["B_DTF_M"]

            a.Snapshot(f"DecayTreeTuple_a", outputfile, clist, opts)
            b.Snapshot(f"DecayTreeTuple_b", outputfile, clist, opts)
            c.Snapshot(f"DecayTreeTuple_c", outputfile, clist, opts)
            d.Snapshot(f"DecayTreeTuple_d", outputfile, clist, opts)
            m.Snapshot(f"DecayTreeTuple_m", outputfile, clist, opts)
            n.Snapshot(f"DecayTreeTuple_n", outputfile, clist, opts)
            x.Snapshot(f"DecayTreeTuple_x", outputfile, clist, opts)
            y.Snapshot(f"DecayTreeTuple_y", outputfile, clist, opts)

            print(f"finished snapshot for {outputfile}")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Apply D Window Cuts Conditions")
    parser.add_argument("--Spec_List", choices = id_to_spec_dict.keys(),  nargs="+", help = 'Spec')
    parser.add_argument("--Output_Type",choices = ["signal_and_sb","d_no_trigger"])
    parser.add_argument('--Snap_Flag', action='store_true')
    parser.add_argument('--Txt_Flag', action='store_true')

    args = parser.parse_args()
    spec_list = args.Spec_List
    output_type = args.Output_Type
    snap_flag = args.Snap_Flag
    txt_flag = args.Txt_Flag

    for spec in spec_list:
        if id_to_spec_dict[spec] == spec:
            type = "data"
        else:
            type = "mc"
        for year in ["2016","2017","2018"]:
            if output_type == "d_no_trigger":
                inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/pre_d/{id_to_spec_dict[spec]}.root"
            if output_type == "signal_and_sb":
                outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_{output_type}/{id_to_spec_dict[spec]}.root"
                inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_trigger/{id_to_spec_dict[spec]}.root"
            apply_dwindow_cuts(id_to_spec_dict[spec], inputfile, outputfile, type, snap_flag, txt_flag)
