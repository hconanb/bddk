import sys
from createXFD import *

sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/")
from essentials import *

print("Max tree size", ROOT.TTree.GetMaxTreeSize())
ROOT.TTree.SetMaxTreeSize(2000000000000000000)
print("Updated tree size", ROOT.TTree.GetMaxTreeSize())

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

trigger_TOS = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)"
trigger_TIS = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"
trigger_HLT1 = "(B_Hlt1TrackMVADecision_TOS == 1 || B_Hlt1TwoTrackMVADecision_TOS == 1)"  # "B_Hlt1Global_TOS == 1"
trigger_HLT2 = "(B_Hlt2Topo2BodyDecision_TOS == 1 || B_Hlt2Topo3BodyDecision_TOS == 1 || B_Hlt2Topo4BodyDecision_TOS == 1)"  # "B_Hlt2Global_TOS == 1"
loose_filters = "abs(KST_M - 895) < 50 && D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && KSTH1_ProbNNk > 0.3"
loose_filters_norm = "D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && K_ProbNNk > 0.3"

fl_dict = {
    "ToT": f"({trigger_TOS} || {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    "nTaT": f"(!{trigger_TOS} && {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    "T": f"({trigger_TOS}) && {trigger_HLT1} && {trigger_HLT2}",
}

mc_spec_list = [
    "01_z_11198006",
    "02_z_11198400",
    "02_p_11198005",
    "04_z_11198401",
    "04_p_11198410",
    "04_zz_11198022",
    "04_st_11198022",
    "05_p_12197023",
    "06_p_12197410",
    "07_p_12197400",
    "07_zz_12197024",
    "07_st_12197024",
    "08_p_12197401",
    "08_zz_12197422",
    "08_st_12197422",
    "09_zz_11196019",
    "10_zz_11196413",
    "12_zz_11196414",
    "13_s_13198040",
    "14_s_13198200",
    "15_s_13198400",
    "16_s_13198600",
    "norm7_norm7_12197008",
    "norm8_norm8_11198007",
]

data_spec_list = ["z"]#, "zz", "p", "m", "st", "s", "norm7", "norm8"]

p0_mc_list = [
    "01_z_11198006",
    "05_p_12197023",
    "09_zz_11196019",
    "13_s_13198040",
    "norm7_norm7_12197008",
    "norm8_norm8_11198007",
]
p1_1e_mc_list = ["02_p_11198005", "06_p_12197410", "10_zz_11196413", "14_s_13198200"]
p1_2e_mc_list = [
    "07_p_12197400",
    "07_zz_12197024",
    "07_st_12197024",
    "15_s_13198400",
]
p1_flip_mc_list = ["02_z_11198400"]
p2_mc_list = [
    "04_z_11198401",
    "04_p_11198410",
    "04_zz_11198022",
    "04_st_11198022",
    "08_p_12197401",
    "08_zz_12197422",
    "08_st_12197422",
    "12_zz_11196414",
    "16_s_13198600",
]

def get_key_con(plist, gd1="_", gd2="_", gd3="_"):
    n = len(plist)
    if n == 2:
        newstring = f"{plist[0]}_MC{gd1}MOTHER_KEY == {plist[1]}_MC{gd2}MOTHER_KEY"
    if n == 3:
        newstring = f"{plist[0]}_MC{gd1}MOTHER_KEY == {plist[1]}_MC{gd2}MOTHER_KEY && {plist[0]}_MC{gd1}MOTHER_KEY == {plist[2]}_MC{gd3}MOTHER_KEY && {plist[1]}_MC{gd2}MOTHER_KEY == {plist[2]}_MC{gd3}MOTHER_KEY"
    if n == 4:
        newstring = f"{plist[0]}_MC_MOTHER_KEY == {plist[1]}_MC_MOTHER_KEY && {plist[0]}_MC_MOTHER_KEY == {plist[2]}_MC_MOTHER_KEY && {plist[0]}_MC_MOTHER_KEY == {plist[3]}_MC_MOTHER_KEY && {plist[1]}_MC_MOTHER_KEY == {plist[2]}_MC_MOTHER_KEY && {plist[1]}_MC_MOTHER_KEY == {plist[3]}_MC_MOTHER_KEY && {plist[2]}_MC_MOTHER_KEY == {plist[3]}_MC_MOTHER_KEY"
    return newstring

def build_truth_strings(name):
    spec = name.split("_")[1]
    b0l = [
        "01_z_11198006",
        "02_z_11198400",
        "02_p_11198005",
        "04_z_11198401",
        "04_p_11198410",
        "04_zz_11198022",
        "04_st_11198022",
        "09_zz_11196019",
        "10_zz_11196413",
        "12_zz_11196414",
        "norm8_norm8_11198007",
    ]
    bpl = [
        "05_p_12197023",
        "06_p_12197410",
        "07_p_12197400",
        "07_zz_12197024",
        "07_st_12197024",
        "08_p_12197401",
        "08_zz_12197422",
        "08_st_12197422",
        "norm7_norm7_12197008",
    ]
    bsl = [
        "13_s_13198040",
        "14_s_13198200",
        "15_s_13198400",
        "16_s_13198600",
    ]
    if name in b0l:
        B_ID = B0_ID
    if name in bpl:
        B_ID = Bp_ID
    if name in bsl:
        B_ID = Bs_ID
    truth_MC = f"abs(B_TRUEID) == {B_ID}"
    if spec == "z":
        c1_ID = Dp_ID
        c2_ID = Dp_ID
    if spec == "zz" or spec == "st" or spec == "norm7":
        c1_ID = D0_ID
        c2_ID = D0_ID
    if spec == "p":
        c1_ID = D0_ID
        c2_ID = Dp_ID
    if spec == "m" or spec == "norm8":
        c1_ID = Dp_ID
        c2_ID = D0_ID
    if spec == "st":
        c1st_ID = Dpst_ID
    if spec == "s":
        c1_ID = Ds_ID
        c2_ID = Dp_ID
    # Conditions that only work if no excited particles
    if name in p0_mc_list:
        truth_MC = (
            f"{truth_MC} && abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {B_ID}"
        )
        truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c1_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {B_ID}"
        if spec != "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {B_ID}"
        if spec == "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {k_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = (
            f"{truth_MC} && abs(D2_TRUEID) == {c2_ID} && abs(D2_MC_MOTHER_ID) == {B_ID}"
        )
        truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}  && abs(D2H1_MC_MOTHER_ID) == {c2_ID} && abs(D2H1_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID} && abs(D2H2_MC_MOTHER_ID) == {c2_ID} && abs(D2H2_MC_GD_MOTHER_ID) == {B_ID}"
        if c1_ID == Dp_ID or c1_ID == Ds_ID:
            truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID} && abs(D1H3_MC_MOTHER_ID) == {c1_ID} && abs(D1H3_MC_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2','D1H3'])}"
        if c2_ID == Dp_ID:
            truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID} && abs(D2H3_MC_MOTHER_ID) == {c2_ID} && abs(D2H3_MC_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2','D2H3'])}"
        if c1_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2'])}"
        if c2_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2'])}"
        if spec == "norm8" or spec == "norm7":
            truth_MC = f"{truth_MC} && abs(D2H4_TRUEID) == {pi_ID} && abs(D2H4_MC_MOTHER_ID) == {c2_ID} && abs(D2H4_MC_GD_MOTHER_ID) == {B_ID}"
        if "norm" not in spec:
            truth_MC = f"{truth_MC} && {get_key_con(['D1','D2','KST'])}"
    if name in p1_1e_mc_list:
        if "02_p_" in name:
            c1_st_ID = Dpst_ID
        if "06_p_" in name or "10_zz_" in name:
            c1_st_ID = D0st_ID
        if "14_s_" in name:
            c1_st_ID = Dsst_ID
        truth_MC = f"{truth_MC} && abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {c1_st_ID} && abs(D1_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c1_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if spec != "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if spec == "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {k_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        truth_MC = (
            f"{truth_MC} && abs(D2_TRUEID) == {c2_ID} && abs(D2_MC_MOTHER_ID) == {B_ID}"
        )
        truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}  && abs(D2H1_MC_MOTHER_ID) == {c2_ID} && abs(D2H1_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID} && abs(D2H2_MC_MOTHER_ID) == {c2_ID} && abs(D2H2_MC_GD_MOTHER_ID) == {B_ID}"
        if c1_ID == Dp_ID or c1_ID == Ds_ID:
            truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID} && abs(D1H3_MC_MOTHER_ID) == {c1_ID} && abs(D1H3_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H3_MC_GD_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2','D1H3'])}"
        if c2_ID == Dp_ID:
            truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID} && abs(D2H3_MC_MOTHER_ID) == {c2_ID} && abs(D2H3_MC_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2','D2H3'])}"
        if c1_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2'])}"
        if c2_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2'])}"
        truth_MC = f"{truth_MC} && {get_key_con(['D1','D2','KST'], gd1 = '_GD_')}"
    if name in p1_2e_mc_list:
        c2_st_ID = Dpst_ID
        truth_MC = (
            f"{truth_MC} && abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {B_ID}"
        )
        truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c1_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {B_ID}"
        if spec != "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {B_ID}"
        if spec == "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {k_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D2_TRUEID) == {c2_ID} && abs(D2_MC_MOTHER_ID) == {c2_st_ID} && abs(D2_MC_GD_MOTHER_ID) == {B_ID} "
        truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}  && abs(D2H1_MC_MOTHER_ID) == {c2_ID} && abs(D2H1_MC_GD_MOTHER_ID) == {c2_st_ID} && abs(D2H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID} && abs(D2H2_MC_MOTHER_ID) == {c2_ID} && abs(D2H2_MC_GD_MOTHER_ID) == {c2_st_ID} && abs(D2H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if c1_ID == Dp_ID or c1_ID == Ds_ID:
            truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID} && abs(D1H3_MC_MOTHER_ID) == {c1_ID} && abs(D1H3_MC_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2','D1H3'])}"
        if c2_ID == Dp_ID:
            truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID} && abs(D2H3_MC_MOTHER_ID) == {c2_ID} && abs(D2H3_MC_GD_MOTHER_ID) == {c2_st_ID} && abs(D2H3_MC_GD_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2','D2H3'])}"
        if c1_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2'])}"
        if c2_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2'])}"
        truth_MC = f"{truth_MC} && {get_key_con(['D1','D2','KST'], gd2 = '_GD_')}"
    if name in p2_mc_list:
        if "04_" in name:
            c1_st_ID = Dpst_ID
            c2_st_ID = Dpst_ID
        if "08_" in name:
            c1_st_ID = D0st_ID
            c2_st_ID = Dpst_ID
        if "12_zz_" in name:
            c1_st_ID = D0st_ID
            c2_st_ID = D0st_ID
        if "16_s_" in name:
            c1_st_ID = Dsst_ID
            c2_st_ID = Dpst_ID
        truth_MC = f"{truth_MC} && abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {c1_st_ID} && abs(D1_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c1_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if spec != "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if spec == "s":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {k_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D2_TRUEID) == {c2_ID} && abs(D2_MC_MOTHER_ID) == {c2_st_ID} && abs(D2_MC_GD_MOTHER_ID) == {B_ID} "
        truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}  && abs(D2H1_MC_MOTHER_ID) == {c2_ID} && abs(D2H1_MC_GD_MOTHER_ID) == {c2_st_ID} && abs(D2H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID} && abs(D2H2_MC_MOTHER_ID) == {c2_ID} && abs(D2H2_MC_GD_MOTHER_ID) == {c2_st_ID} && abs(D2H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if c1_ID == Dp_ID or c1_ID == Ds_ID:
            truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID} && abs(D1H3_MC_MOTHER_ID) == {c1_ID} && abs(D1H3_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H3_MC_GD_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2','D1H3'])}"
        if c2_ID == Dp_ID:
            truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID} && abs(D2H3_MC_MOTHER_ID) == {c2_ID} && abs(D2H3_MC_GD_MOTHER_ID) == {c2_st_ID} && abs(D2H3_MC_GD_GD_MOTHER_ID) == {B_ID}"
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2','D2H3'])}"
        if c1_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2'])}"
        if c2_ID == D0_ID:
            truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2'])}"
        truth_MC = f"{truth_MC} && {get_key_con(['D1','D2','KST'], gd1 = '_GD_', gd2 = '_GD_')}"
    if name in p1_flip_mc_list:
        c1_st_ID = Dpst_ID
        c2_st_ID = Dpst_ID

        f1_d = f"abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {c1_st_ID} && abs(D1_MC_GD_MOTHER_ID) == {B_ID} && abs(D2_TRUEID) == {c2_ID} && abs(D2_MC_MOTHER_ID) == {B_ID}"
        f2_d = f"abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {B_ID} && abs(D2_TRUEID) == {c2_ID} && abs(D2_MC_MOTHER_ID) == {c2_st_ID} && abs(D2_MC_GD_MOTHER_ID) == {B_ID} "
        truth_MC = f"{truth_MC} && (({f1_d}) || ({f2_d}))"

        f1_d1h1 = f"abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c1_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        f1_d1h2 = f"abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        f1_d1h3 = f"abs(D1H3_TRUEID) == {pi_ID} && abs(D1H3_MC_MOTHER_ID) == {c1_ID} && abs(D1H3_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H3_MC_GD_GD_MOTHER_ID) == {B_ID}"
        f1_d2h1 = f"abs(D2H1_TRUEID) == {k_ID}  && abs(D2H1_MC_MOTHER_ID) == {c2_ID} && abs(D2H1_MC_GD_MOTHER_ID) == {B_ID}"
        f1_d2h2 = f"abs(D2H2_TRUEID) == {pi_ID} && abs(D2H2_MC_MOTHER_ID) == {c2_ID} && abs(D2H2_MC_GD_MOTHER_ID) == {B_ID}"
        f1_d2h3 = f"abs(D2H3_TRUEID) == {pi_ID} && abs(D2H3_MC_MOTHER_ID) == {c2_ID} && abs(D2H3_MC_GD_MOTHER_ID) == {B_ID}"

        f2_d1h1 = f"abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c2_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {B_ID}"
        f2_d1h2 = f"abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c2_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {B_ID}"
        f2_d1h3 = f"abs(D1H3_TRUEID) == {pi_ID} && abs(D1H3_MC_MOTHER_ID) == {c2_ID} && abs(D1H3_MC_GD_MOTHER_ID) == {B_ID}"
        f2_d2h1 = f"abs(D2H1_TRUEID) == {k_ID}  && abs(D2H1_MC_MOTHER_ID) == {c1_ID} && abs(D2H1_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D2H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        f2_d2h2 = f"abs(D2H2_TRUEID) == {pi_ID} && abs(D2H2_MC_MOTHER_ID) == {c1_ID} && abs(D2H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D2H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        f2_d2h3 = f"abs(D2H3_TRUEID) == {pi_ID} && abs(D2H3_MC_MOTHER_ID) == {c1_ID} && abs(D2H3_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D2H3_MC_GD_GD_MOTHER_ID) == {B_ID}"

        f1_d1h = f"{f1_d1h1} && {f1_d1h2} && {f1_d1h3} && {f1_d2h1} && {f1_d2h2} && {f1_d2h3}"
        f2_d2h = f"{f2_d1h1} && {f2_d1h2} && {f2_d1h3} && {f2_d2h1} && {f2_d2h2} && {f2_d2h3}"

        truth_MC = f"{truth_MC} && (({f1_d1h}) || ({f2_d2h}))"
        truth_MC = f"{truth_MC} && {get_key_con(['D1H1','D1H2','D1H3'])}"
        truth_MC = f"{truth_MC} && {get_key_con(['D2H1','D2H2','D2H3'])}"

        f1key = get_key_con(["D1", "D2", "KST"], gd1="_GD_")
        f2key = get_key_con(["D1", "D2", "KST"], gd2="_GD_")

        truth_MC = f"{truth_MC} && (({f1key}) || ({f2key}))"

    # #These Conditions should work for all spectrum
    if "norm" not in spec:
        truth_MC = f"{truth_MC} && abs(KST_TRUEID) == {kst0_ID} && abs(KST_MC_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(KSTH1_TRUEID) == {k_ID} && abs(KSTH1_MC_MOTHER_ID) == {kst0_ID} && abs(KSTH1_MC_GD_MOTHER_ID) == {B_ID} && abs(KSTH2_TRUEID) == {pi_ID} && abs(KSTH2_MC_MOTHER_ID) == {kst0_ID} && abs(KSTH2_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && {get_key_con(['KSTH1','KSTH2'])}"
    if "norm" in spec:
        truth_MC = (
            f"{truth_MC} && abs(K_TRUEID) == {k_ID} && abs(K_MC_MOTHER_ID) == {B_ID}"
        )
        truth_MC = f"{truth_MC} && {get_key_con(['D1','D2','K'])}"
    return truth_MC

def convertTuple(tup, spec):
    pname = ""
    for p in tup:
        pname += p
    return pname

def invmass(plist):
    """arguments:
    plist -- list of particles, e.g., ["K1", "K2", ...]
    """
    sflag = 0
    for p in plist:
        if sflag == 0:
            e = f"{p}_PE"
            x = f"{p}_PX"
            y = f"{p}_PY"
            z = f"{p}_PZ"
        if sflag == 1:
            e = f"{e}+{p}_PE"
            x = f"{x}+{p}_PX"
            y = f"{y}+{p}_PY"
            z = f"{z}+{p}_PZ"
        sflag = 1
    e2 = f"pow({e},2)"
    x2 = f"pow({x},2)"
    y2 = f"pow({y},2)"
    z2 = f"pow({z},2)"
    s = f"(sqrt({e2}-{x2}-{y2}-{z2}))"
    return s

def get_dwindow_values(dwindow_ws, tag):

    if tag == "z" or tag == "zz" or tag == "st" or tag == "norm7":

        d1_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
        d2_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
        d1_std = dwindow_ws.var(f"width_a_{tag}_D1").getValV()
        d2_std = dwindow_ws.var(f"width_a_{tag}_D1").getValV()

    if tag == "p" or tag == "m" or tag == "s" or tag == "norm8":

        d1_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
        d2_mstart = dwindow_ws.var(f"mean_{tag}_D2").getValV()

        d1_std = dwindow_ws.var(f"width_a_{tag}_D1").getValV()
        d2_std = dwindow_ws.var(f"width_a_{tag}_D2").getValV()

    if tag == "st":
        d3_mstart = dwindow_ws.var(f"mean_{tag}_D3").getValV()
        d3_std = dwindow_ws.var(f"width_a_{tag}_D3").getValV()
        d3window = d3_std * 2.5

    d1window = d1_std * 2
    d2window = d2_std * 2

    if tag == "st":
        mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window} && abs(DstmD_M - {d3_mstart}) < {d2window}) "
    else:
        mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window})"

    return mass_cut

def make_folders():
    base_old = "/mnt/c/Users/Harris/Desktop/rootfiles/raw_mc_2021"
    old_path_list = glob.glob(f"{base_old}/*/*/ntuple.root")
    for mc_name in mc_name_list:
        new_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{mc_name}/"
        if os.path.isdir(new_path) == False:
            os.makedirs(new_path)
        for old_path in old_path_list:
            if os.path.isfile(old_path):
                temp_file = ROOT.TFile.Open(old_path)
                if temp_file.GetListOfKeys().Contains(mc_name):
                    temp_file.Close()
                    temp_path = old_path.split("/")[-3]
                    if os.path.isdir(f"{base_old}/{temp_path}"):
                        shutil.move(f"{base_old}/{temp_path}", new_path)
                else:
                    temp_file.Close()

def rename_folders():
    folder_path_list = glob.glob(
        f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/*/*"
    )
    num_names_df = pandas.read_csv(
        "ganga_names.txt", sep=" ", engine="python", header=None
    )
    column_list = ["ganga_n", "new_name"]
    num_names_df.columns = column_list
    print(num_names_df)
    for name in folder_path_list:
        num = name.rsplit("/")[-1]
        base = name.rsplit("/", 1)[-2]
        newfoldername = num_names_df.loc[num_names_df["ganga_n"] == int(num)].iat[0, 1]
        print(f"move {name} to {base}/{newfoldername}")
        shutil.move(name, f"{base}/{newfoldername}")

def filter(name, type_flag, year, trigger_flag):

    if type_flag == "DATA":
        spec = name
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_data/1*/*/ntuple.root")
        tree_chain = ROOT.TChain(f"data_{name}_Data/DecayTreeTuple")
        file = ROOT.TFile(
            f"/mnt/c/Users/Harris/Desktop/rootfiles/data_run2/{spec}_spectrum.root"
        )
        tc = ROOT.TChain(f"data_{name}_Data/DecayTreeTuple")
        if "_st_" not in name:
            for file_name in file_list:
                tc.Add(file_name)
                tcmc.Add(file_name)
        rdf_data_base = RDF(tree)


    if type_flag == "MC":
        spec = name.split("_")[1]
        if "_st_" in name:
            newname = name.replace("_st_", "_zz_")
            file_list = glob.glob(
                f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{newname}/*_{year}_*/*/ntuple.root"
            )
            for file in file_list:
                CreateXFD(
                    file,
                    f"{name}/DecayTreeTuple",
                    file.replace("ntuple.root", "ntuple_st_fix.root"),
                    f"{name}/DecayTreeTuple",
                    mc_flag=1,
                )
            rec_file_list = glob.glob(
                f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{newname}/*_{year}_*/*/ntuple_st_fix.root"
            )
            mc_file_list = file_list
        else:
            file_list = glob.glob(
                f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{name}/*_{year}_*/*/ntuple.root"
            )
        tc = ROOT.TChain(f"{name}/DecayTreeTuple")
        tcmc = ROOT.TChain("MCDecayTreeTuple/MCDecayTree")
        if "_st_" not in name:
            for file_name in file_list:
                tc.Add(file_name)
                tcmc.Add(file_name)
                # print(f"added {file_name}")
        if "_st_" in name:
            for rec_file_name, mc_file_name in zip(rec_file_list, mc_file_list):
                tc.Add(rec_file_name)
                tcmc.Add(mc_file_name)
                # print(f"added {rec_file_name} and {mc_file_name}")
        rdf_mc_base = RDF(tc)
        rdf_MC_GEN = RDF(tcmc)

    # print(f"{name} reconstructed : {rdf_mc_base.Count().GetValue()}")
    # print(f"{name} gen : {rdf_MC_GEN.Count().GetValue()}")

    ### Add MC Lines ####
    TRUE_PE_REC = "(pow(D1_TRUEP_E + D2_TRUEP_E + KST_TRUEP_E,2))"
    TRUE_PX_REC = "(pow(D1_TRUEP_X + D2_TRUEP_X + KST_TRUEP_X,2))"
    TRUE_PY_REC = "(pow(D1_TRUEP_Y + D2_TRUEP_Y + KST_TRUEP_Y,2))"
    TRUE_PZ_REC = "(pow(D1_TRUEP_Z + D2_TRUEP_Z + KST_TRUEP_Z,2))"
    TRUE_BM_REC = (f"sqrt({TRUE_PE_REC} - ({TRUE_PX_REC} + {TRUE_PY_REC} + {TRUE_PZ_REC}))")

    if "norm" not in spec and "st" not in spec and type_flag == "MC":
        rdf_rec = (
            rdf_mc_base.Define("B_DTF_M", "B_dtf_M[0]")
            .Define("TRUE_BM_REC", TRUE_BM_REC)
            .Define("Res", f"B_DTF_M - TRUE_BM_REC")
        )
    if "st" in spec and type_flag == "MC":
        rdf_rec = (
            rdf_mc_base.Define("B_DTF_M", "B_dtf_M")
            .Define("TRUE_BM_REC", TRUE_BM_REC)
            .Define("Res", f"B_DTF_M - TRUE_BM_REC")
        )
    if "norm" in spec and type_flag == "MC":
        rdf_rec = rdf_mc_base.Define("B_DTF_M", "B_dtf_M[0]")

    if type_flag == "DATA":
        rdf_rec = rdf_data_base
    #### Add Submass Branches ####
    if "norm" not in spec:
        tl = ["D1", "D2", "KSTH1", "KSTH2"]
        nlist = [2, 3, 4]
        for n in nlist:
            # Final_n_List = []
            all_n_combinations = itertools.combinations(tl, n)
            tupflag = 0
            for tup in all_n_combinations:
                tname = convertTuple(tup, spec)
                rdf_temp = rdf_rec.Define(tname, invmass(tup))
                rdf_rec = rdf_temp
    ## Get D Window Cuts
    # dwindow_file = ROOT.TFile(f"d_window_root_files/d_{spec}_mass_fits.root", "READ")
    # dwindow_ws = dwindow_file.Get(f"d_{spec}_mass_fits")
    # dwindow_cut = get_dwindow_values(dwindow_ws, spec)
    # dwindow_file.Close()
    # print(dwindow_cut)
    #### Get BKG Cut
    if spec == "zz":
        bkg_cut = "(D1KSTH2 > 2050)"
        ###Add st cut
    else:
        bkg_cut = "(D1_M > 0)"
    #### Get Trigger and Loose Cuts
    trigger_cut = fl_dict[f"{trigger_flag}"]
    if "norm" not in spec:
        loose_cut = loose_filters
    if "norm" in spec:
        loose_cut = loose_filters_norm
    ### Dira Cuts
    if "st" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX > 0)"
    if "st" in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX_FIXED > 0)"
    ### Apply All Cuts
    all_cuts_string = (
        f"{bkg_cut} && {dira_cut} && {loose_cut} && {trigger_cut}"
        # f"{bkg_cut} && {dwindow_cut} && {dira_cut} && {loose_cut} && {trigger_cut}"

    )
    # all_cuts_string = bkg_cut
    if type_flag == "MC":
        all_cuts_string = f"{all_cuts_string} && {build_truth_strings(name)}"
    # print(all_cuts_string)
    rdf_final = rdf_rec.Filter(all_cuts_string)

    # print(f"{name}_{year}_{trigger_flag} reconstructed after cuts: {rdf_final.Count().GetValue()}")
    # print(f"{name}_{year}_{trigger_flag} generated: {rdf_MC_GEN.Count().GetValue()}")
    # print(f"eff: {rdf_final.Count().GetValue()/rdf_MC_GEN.Count().GetValue()}")

    if type_flag == "DATA":
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root_nod/{type_flag}/{name}_spectrum_filtered.root"
    else:
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type_flag}/{name}_{year}_{trigger_flag}_spectrum_filtered.root"

    if os.path.exists(outputfile):
        os.remove(outputfile)
        print("First deleting old ", outputfile)
    else:
        print("making ", outputfile)

    clist = rdf_final.GetColumnNames()

    rdfsnap = rdf_final.Snapshot(f"DecayTreeTuple", outputfile, clist, opts)
    if type_flag == "MC":
        mc_clist = rdf_MC_GEN.GetColumnNames()
        rdfsnap_MC = rdf_MC_GEN.Snapshot(
            f"MCDecayTreeTuple", outputfile, mc_clist, opts
        )
    print(f"finished snapshot for {spec}")

def zz_olapfit(name, year, trigger):

    file_base = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type_flag}/{spec_name}_{year}_{trigger_flag}_spectrum_filtered.root")

    tree_zz = file_base.Get(f"DecayTreeTuple")

    st_name = spec_name.replace("_zz_","_st_")

    tree_st = file_base.Get(f"DecayTreeTuple")

    rdf_zz = RDF(tree_zz)
    rdf_st = RDF(tree_st)

    np_zz = rdf_zz.AsNumpy()
    np_st = rdf_st.AsNumpy()

    zz_en_base = np.char.mod('%d', np_zz["eventNumber"])
    zz_nCan_base = np.char.mod('%d', np_zz["runNumber"])

    st_en_base = np.char.mod('%d', np_st["eventNumber"])
    st_nCan_base = np.char.mod('%d', np_st["runNumber"])

    zz_all = np.core.defchararray.add(np.core.defchararray.add(zz_en_base, "_"),zz_nCan_base)
    st_all = np.core.defchararray.add(np.core.defchararray.add(st_en_base, "_"),st_nCan_base)

    olap_list = np.intersect1d(zz_all,st_all)
    fline = "B_M > 0"
    for olap_entry in olap_list:
        eno = olap_entry.split("_")[0]
        ncan = olap_entry.split("_")[1]
        fline = f"{fline} && (eventNumber != {eno} || (eventNumber == {eno} && nCandidate != {ncan}))"
         # if old_evvent.eventNumer != eno or (old_event.eventNumber = eno and old_event.nCandidate != ncan):
    rdf_zz_next = rdf_zz.Filter(fline)

    print(len(zz_en_base))
    print(len(olap_list))
    print(rdf_zz_next.Count().GetValue())

    clist = rdf_zz_next.GetColumnNames()

    if os.path.exists(file_name.replace("ntuple.root", "ntuple_zz_fix.root")):
        os.remove(file_name.replace("ntuple.root", "ntuple_zz_fix.root"))
        print("First deleting old ", file_name.replace("ntuple.root", "ntuple_zz_fix.root"))
    else:
        print("making ", file_name.replace("ntuple.root", "ntuple_zz_fix.root"))
    rdfsnap = rdf_zz_next.Snapshot(f"{spec_name}/DecayTreeTuple", file_name.replace("ntuple.root", "ntuple_zz_fix.root"), clist, opts)
#
# for name in mc_spec_list:
#     for year in ["2016", "2017", "2018"]:
#         for trigger in ["T", "nTaT"]:
#             filter(name, "MC", year, trigger)
#             if name in ["04_zz_11198022","07_zz_12197024","08_zz_12197422"]:
#                 zz_olapfit(name, year, trigger)

for spec in data_spec_list:
    for year in ["2016", "2017", "2018"]:
        for trigger in ["T", "nTaT"]:
            filter(name, "Data", year, trigger)

for spec in data_spec_list:
    print(f"Building {spec} ntuple")
    filter(spec, "DATA", "none", "ToT")
