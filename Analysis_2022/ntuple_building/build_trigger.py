import sys
import os

basedir = os.getcwd().split('ntuple_building')[0]
sys.path.append(basedir)

from essential_functions import *

RDF = ROOT.ROOT.RDataFrame

trigger_TOS_ALL_L0 = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)"
trigger_TOS_H_L0 = "(B_L0HadronDecision_TOS == 1)"

trigger_TIS_ALL_L0 = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"

trigger_TOS_HLT1 = "(B_Hlt1TrackMVADecision_TOS == 1 || B_Hlt1TwoTrackMVADecision_TOS == 1)"
trigger_TIS_HLT1 = "(B_Hlt1TrackMVADecision_TIS == 1 || B_Hlt1TwoTrackMVADecision_TIS == 1)"

trigger_TOS_HLT2 = "(B_Hlt2Topo2BodyDecision_TOS == 1 || B_Hlt2Topo3BodyDecision_TOS == 1 || B_Hlt2Topo4BodyDecision_TOS == 1)"
trigger_TIS_HLT2 = "(B_Hlt2Topo2BodyDecision_TIS == 1 || B_Hlt2Topo3BodyDecision_TIS == 1 || B_Hlt2Topo4BodyDecision_TIS == 1)"

def apply_trigger(spec, type, year, snap_flag, txt_flag):

    inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/pre_d/{spec}.root"
    outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_trigger/{spec}.root"

    base_file = ROOT.TFile(inputfile, "READ")

    base_tree = base_file.Get(f"DecayTreeTuple")

    rdf_temp = RDF(base_tree)

    rdf_base = rdf_numbers(rdf_temp, "No Trigger")
    rdf_trigger_h = rdf_base.apply_filter(f"B_L0HadronDecision_TOS == 1","L0_TOS_H")
    rdf_trigger_m = rdf_base.apply_filter(f"!(B_L0HadronDecision_TOS == 1) && B_L0MuonDecision_TOS == 1","L0_TOS_M")
    rdf_trigger_p = rdf_base.apply_filter(f"!(B_L0HadronDecision_TOS == 1) && B_L0PhotonDecision_TOS == 1","L0_TOS_P")
    rdf_trigger_e = rdf_base.apply_filter(f"!(B_L0HadronDecision_TOS == 1) && B_L0ElectronDecision_TOS == 1","L0_TOS_E")

    # rdf_trigger_L0_tos_ALL = rdf_base.apply_filter(f"{trigger_TOS_ALL_L0} || {trigger_TIS_ALL_L0}","L0_TOS_ALL")
    # rdf_trigger_L0_tos_h = rdf_base.apply_filter(f"{trigger_TOS_H_L0} || {trigger_TIS_ALL_L0}","L0_TOS_H")
    #
    # rdf_trigger_HLT1_tos = rdf_trigger_L0_tos_ALL.apply_filter(f"{trigger_TOS_HLT1}","HLT1_TOS")
    # rdf_trigger_HLT2_tos = rdf_trigger_HLT1_tos.apply_filter(f"{trigger_TOS_HLT2}","HLT2_TOS")

    print('All stats:')
    allCutsReport = rdf_base.rdf.Report()
    allCutsReport.Print()

    if txt_flag:
            eff_dict_trigger_temp = {
                "Scheme ID" : id_to_scheme_dict[spec],
                "Year": year,
                "$\epsilon_{trigger} ALL$": f"{rdf_trigger_L0_tos_ALL.calc_can_eff()*100.000:.3f}",
                "$\epsilon_{trigger} H$": f"{rdf_trigger_L0_tos_h.calc_can_eff()*100.000:.3f}",
                "$\epsilon_{trigger} HLT1$": f"{rdf_trigger_HLT1_tos.calc_can_eff()*100.000:.3f}",
                "$\epsilon_{trigger} HLT2$": f"{rdf_trigger_HLT2_tos.calc_can_eff()*100.000:.3f}",
            }
            eff_df = pd.DataFrame(eff_dict_trigger_temp , index=[0])
            eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/trigger/can_{spec}_{year}.txt")

            eff_dict_trigger_temp = {
                "Scheme ID" : id_to_scheme_dict[spec],
                "Year": year,
                "$\epsilon_{trigger} ALL$": f"{rdf_trigger_L0_tos_ALL.calc_event_eff()*100.000:.3f}",
                "$\epsilon_{trigger} H$": f"{rdf_trigger_L0_tos_h.calc_event_eff()*100.000:.3f}",
                "$\epsilon_{trigger} HLT1$": f"{rdf_trigger_HLT1_tos.calc_event_eff()*100.000:.3f}",
                "$\epsilon_{trigger} HLT2$": f"{rdf_trigger_HLT2_tos.calc_event_eff()*100.000:.3f}",
            }
            eff_df = pd.DataFrame(eff_dict_trigger_temp , index=[0])
            eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/trigger/event_{spec}_{year}.txt")

    if snap_flag:
        clist = rdf_trigger_HLT2_tos.rdf.GetColumnNames()
        clist_f = ["B_dtf_c_nPV","nPV", "B_dtf_nPV"]
        for name in clist:
            if name not in clist_f:
                clist_f.append(name)
        print(f"Starting snapshot for {outputfile}")
        rdfsnap = rdf_trigger_HLT2_tos.rdf.Snapshot(f"DecayTreeTuple", outputfile, clist_f)
        print(f"finished snapshot for {outputfile}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Apply Trigger Conditions")
    parser.add_argument("--Spec_List", choices = id_to_spec_dict.keys(),  nargs="+", help = 'Spec')
    parser.add_argument('--Snap_Flag', action='store_true')
    parser.add_argument('--Txt_Flag', action='store_true')

    args = parser.parse_args()
    spec_list = args.Spec_List
    snap_flag = args.Snap_Flag
    txt_flag = args.Txt_Flag

    for spec in spec_list:
        if id_to_spec_dict[spec] == spec:
            type = "data"
        else:
            type = "mc"
        for year in ["2016","2017","2018"]:
            apply_trigger(id_to_spec_dict[spec], type, year, snap_flag, txt_flag)
            break
