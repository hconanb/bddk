import ROOT
import os, sys, shutil
import glob as glob
import pandas
import itertools

print("Max tree size", ROOT.TTree.GetMaxTreeSize())
ROOT.TTree.SetMaxTreeSize(2000000000000000000)
print("Updated tree size", ROOT.TTree.GetMaxTreeSize())

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/")
from essentials import *

trigger_TOS = "B_L0Global_TOS == 1"
trigger_TIS = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"
trigger_HLT1 = "B_Hlt1Global_TOS == 1"
trigger_HLT2 = "B_Hlt2Global_TOS == 1"


def get_truth_con(name, id, m_id):
        newstring = f"abs({name}_TRUEID) == {id} && abs({name}_MC_MOTHER_ID) == {m_id}"
        return newstring
def get_key_con(p1,p2,p3,p4,n):
    if n == 2:
        newstring = f"{p1}_MC_MOTHER_KEY == {p2}_MC_MOTHER_KEY"
    if n == 3:
        newstring = f"{p1}_MC_MOTHER_KEY == {p2}_MC_MOTHER_KEY && {p1}_MC_MOTHER_KEY == {p3}_MC_MOTHER_KEY && {p2}_MC_MOTHER_KEY == {p3}_MC_MOTHER_KEY"
    if n == 4:
        newstring = f"{p1}_MC_MOTHER_KEY == {p2}_MC_MOTHER_KEY && {p1}_MC_MOTHER_KEY == {p3}_MC_MOTHER_KEY && {p1}_MC_MOTHER_KEY == {p4}_MC_MOTHER_KEY && {p2}_MC_MOTHER_KEY == {p3}_MC_MOTHER_KEY && {p2}_MC_MOTHER_KEY == {p4}_MC_MOTHER_KEY && {p3}_MC_MOTHER_KEY == {p4}_MC_MOTHER_KEY"
    return newstring


truth_MC_B0 = "abs(B_TRUEID) == 511"
truth_MC_B0normkey = get_key_con("D1","D2","K","",3)
truth_MC_B0sigkey = get_key_con("D1","D2","KST","",3)

truth_MC_Dp = get_truth_con("D1", 411, 511)
truth_MC_DpH1 = get_truth_con("D1H1", 321, 411)
truth_MC_DpH2 = get_truth_con("D1H2", 211, 411)
truth_MC_DpH3 = get_truth_con("D1H3", 211, 411)
truth_MC_Dpkey = get_key_con("D1H1","D1H2","D1H3","",3)
truth_MC_Dp = f"{truth_MC_Dp} && {truth_MC_DpH1} && {truth_MC_DpH2} && {truth_MC_DpH3} && {truth_MC_Dpkey}"

truth_MC_D0 = get_truth_con("D2", 421, 511)
truth_MC_D0H1 = get_truth_con("D2H1", 321, 421)
truth_MC_D0H2 = get_truth_con("D2H2", 211, 421)
truth_MC_D0H3 = get_truth_con("D2H3", 211, 421)
truth_MC_D0H4 = get_truth_con("D2H4", 211, 421)
truth_MC_D02key = get_key_con("D2H1","D2H2","","",2)
truth_MC_D04key = get_key_con("D2H1","D2H2","D2H3","D2H4",4)
truth_MC_D02 = f"{truth_MC_D0} && {truth_MC_D0H1} && {truth_MC_D0H2}"
truth_MC_D04 = f"{truth_MC_D0} && {truth_MC_D0H1} && {truth_MC_D0H2} && {truth_MC_D0H3} && {truth_MC_D0H4} && {truth_MC_D04key}"

# truth_MC_KSTAR = get_truth_con("KST", 313, 511)
# truth_MC_K_KSTAR = get_truth_con("KST", 321, 313)
# truth_MC_pi_KSTAR = get_truth_con("KST", 211, 313)
truth_MC_K_B = get_truth_con("K", 321, 511)
# truth_MC_KstKey = get_key_con("","D2H2","","",2)


all_cuts_string = f"{truth_MC_B0} && {truth_MC_B0normkey} && {truth_MC_Dp} && {truth_MC_D04} && {truth_MC_K_B}"
# loose_filters = "abs(KST_M - 895) < 50 && D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && KSTH1_ProbNNk > 0.3"
# loose_filters_norm = "D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && K_ProbNNk > 0.3"

mc_spec_list = [
    # "01_z_11198006",
    # "02_z_11198400",
    # "02_p_11198005",
    # "04_z_11198401",
    # "04_p_11198410",
    # "04_zz_11198022",
    # "04_st_11198022",
    # "05_p_12197023",
    # "06_p_12197410",
    # "07_p_12197400",
    # "07_zz_12197024",
    # "07_st_12197024",
    # "08_p_12197401",
    # "08_zz_12197422",
    # "08_st_12197422",
    # "09_zz_11196019",
    # "10_zz_11196413",
    # "12_zz_11196414",
    # "13_s_13198040",
    # "14_s_13198200",
    # "15_s_13198400",
    # "16_s_13198600",
    # "norm7_norm7_12197008",
    "norm8_norm8_11198007",
]


def filter(name, type_flag, year, trigger_flag):

    if type_flag == "MC":
        spec = name.split("_")[1]
        if "_st_" in name:
            newname = name.replace("_st_", "_zz_")
            file_list = glob.glob(
                f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{newname}/*_{year}_*/*/ntuple.root"
            )
        else:
            file_list = glob.glob(
                f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{name}/*_{year}_*/*/ntuple.root"
            )
        tc = ROOT.TChain(f"{name}/DecayTreeTuple")
        tcmc = ROOT.TChain("MCDecayTreeTuple/MCDecayTree")
        for file_name in file_list:
            tc.Add(file_name)
            tcmc.Add(file_name)
            print(f"added {file_name}")
        rdf_mc_base = RDF(tc)
        rdf_MC_GEN = RDF(tcmc)
    #### Add MC Lines ####

    # trigger_cut = fl_dict[f"{trigger_flag}"]
    # if "norm" not in spec:
    #     loose_cut = loose_filters
    # if "norm" in spec:
    #     loose_cut = loose_filters_norm

    ### Dira Cuts

    print(all_cuts_string)
    rdf_final = rdf_mc_base.Filter(
        all_cuts_string
    )

    outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_test_root/test_spectrum_filtered_{year}.root"
    if os.path.exists(outputfile):
        os.remove(outputfile)
        print("First deleting old ", outputfile)
    else:
        print("making ", outputfile)

    clist = rdf_final.GetColumnNames()
    rdfsnap = rdf_final.Snapshot(f"DecayTreeTuple", outputfile, clist, opts)
    print(f"finished snapshot for {spec}")

for name in mc_spec_list:
    for year in ["2016","2017","2018"]:
        for trigger in ["T","nTaT"]:
            filter(name, "MC", year, trigger)
