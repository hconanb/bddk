import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/")
from essentials import *
from createXFD import *

# print("Max tree size", ROOT.TTree.GetMaxTreeSize())
# ROOT.TTree.SetMaxTreeSize(2000000000000000000)
# print("Updated tree size", ROOT.TTree.GetMaxTreeSize())

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
data_spec_list = [
    # "Z_m_p",
    # "Zs_sm_p",
    "Z_mst_p",
    "Z_m_pst",
    # "Z_mst_pst",
    "Z_z_z",
    "P_z_p",
    "M_m_z",
    "P_z_pst",
    "M_mst_z",
    "norm7",
    # "norm8"
    ]
mc_spec_list = [
    # "01_Z_m_p_11198006",
    # "02_Z_m_p_11198400",
    # "02_P_z_p_11198005",
    "02_Z_mst_p_11198005",
    # "04_Z_m_p_11198401",
    # "04_P_z_p_11198410",
    "04_Z_mst_p_11198410",
    # "04_Z_z_z_11198022",
    # "04_P_z_pst_11198022",
    # "05_P_z_p_12197023",
    # "06_P_z_p_12197410",
    # "07_P_z_p_12197400",
    # "07_Z_z_z_12197024",
    # "07_P_z_pst_12197024",
    # "08_P_z_p_12197401",
    # "08_Z_z_z_12197422",
    # "08_P_z_pst_12197422",
    # "09_Z_z_z_11196019",
    # "10_Z_z_z_11196413",
    # "12_Z_z_z_11196414",
    "13_Zs_sm_p_13198040",
    "14_Zs_sm_p_13198200",
    "15_Zs_sm_p_13198400",
    "16_Zs_sm_p_13198600",
    # "norm7_norm7_12197008",
    # "norm8_norm8_11198007",
]
p0_mc_list = [
    "01_Z_m_p_11198006",
    "05_P_z_p_12197023",
    "09_Z_z_z_11196019",
    "13_Zs_sm_p_13198040",
    "norm7_norm7_12197008",
    "norm8_norm8_11198007",
]
p1_1e_mc_list = [
"02_P_z_p_11198005", "06_P_z_p_12197410", "10_Z_z_z_11196413", "14_Zs_sm_p_13198200","02_Z_mst_p_11198005","04_Z_mst_p_11198410"]
p1_2e_mc_list = [
    "07_P_z_p_12197400",
    "07_Z_z_z_12197024",
    "07_P_z_pst_12197024",
    "15_Zs_sm_p_13198400",
]
p1_flip_mc_list = [
"02_Z_m_p_11198400"]
p2_mc_list = [
    "04_Z_m_p_11198401",
    "04_P_z_p_11198410",
    "04_Z_z_z_11198022",
    "04_P_z_pst_11198022",
    "08_P_z_p_12197401",
    "08_Z_z_z_12197422",
    "08_P_z_pst_12197422",
    "12_Z_z_z_11196414",
    "16_Zs_sm_p_13198600",
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
    na = name.split("_")
    if "norm" not in name:
        spec = f"{na[1]}_{na[2]}_{na[3]}"
    if "norm" in name:
        spec = name.split("_")[1]
    print(spec)
    b0l = [
        "01_Z_m_p_11198006",
        "02_Z_m_p_11198400",
        "02_P_z_p_11198005",
        "04_Z_m_p_11198401",
        "04_P_z_p_11198410",
        "04_Z_z_z_11198022",
        "04_P_z_pst_11198022",
        "09_Z_z_z_11196019",
        "10_Z_z_z_11196413",
        "12_Z_z_z_11196414",
        "norm8_norm8_11198007",
        "02_Z_mst_p_11198005",
        "04_Z_mst_p_11198410",
    ]
    bpl = [
        "05_P_z_p_12197023",
        "06_P_z_p_12197410",
        "07_P_z_p_12197400",
        "07_Z_z_z_12197024",
        "07_P_z_pst_12197024",
        "08_P_z_p_12197401",
        "08_Z_z_z_12197422",
        "08_P_z_pst_12197422",
        "norm7_norm7_12197008",
    ]
    bsl = [
        "13_Zs_sm_p_13198040",
        "14_Zs_sm_p_13198200",
        "15_Zs_sm_p_13198400",
        "16_Zs_sm_p_13198600",
    ]
    if name in b0l:
        B_ID = B0_ID
    if name in bpl:
        B_ID = Bp_ID
    if name in bsl:
        B_ID = Bs_ID
    truth_MC = f"abs(B_TRUEID) == {B_ID}"
    if spec == "Z_m_p":
        c1_ID = Dp_ID
        c2_ID = Dp_ID
    if spec == "Z_z_z" or spec == "P_z_pst" or spec == "norm7":
        c1_ID = D0_ID
        c2_ID = D0_ID
    if spec == "P_z_p" or spec == "Z_mst_p":
        c1_ID = D0_ID
        c2_ID = Dp_ID
    if spec == "M_m_z" or spec == "norm8":
        c1_ID = Dp_ID
        c2_ID = D0_ID
    if spec == "P_z_pst" or spec == "Z_mst_p":
        c1_st_ID = Dpst_ID
    if spec == "Zs_sm_p":
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
        if "02_P_z_p_" in name:
            c1_st_ID = Dpst_ID
        if "06_P_z_p_" in name or "10_Z_z_z_" in name:
            c1_st_ID = D0st_ID
        if "14_Zs_sm_p_" in name:
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
        if "Zs" not in spec:
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {B_ID}"
        if "Zs" in spec:
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
        if "12_Z_z_z_" in name:
            c1_st_ID = D0st_ID
            c2_st_ID = D0st_ID
        if "16_Zs_" in name:
            c1_st_ID = Dsst_ID
            c2_st_ID = Dpst_ID
        truth_MC = f"{truth_MC} && abs(D1_TRUEID) == {c1_ID} && abs(D1_MC_MOTHER_ID) == {c1_st_ID} && abs(D1_MC_GD_MOTHER_ID) == {B_ID}"
        truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}  && abs(D1H1_MC_MOTHER_ID) == {c1_ID} && abs(D1H1_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H1_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if spec != "Zs_sm_p":
            truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID} && abs(D1H2_MC_MOTHER_ID) == {c1_ID} && abs(D1H2_MC_GD_MOTHER_ID) == {c1_st_ID} && abs(D1H2_MC_GD_GD_MOTHER_ID) == {B_ID}"
        if spec == "Zs_sm_p":
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

def filter(file_list, spec, year, trigger):
    print(spec)
    if spec not in mc_spec_list:
        type = "DATA"
        print(type)
        tree_name = f"data_{spec}_Data/DecayTreeTuple"
        tree_chain = ROOT.TChain(tree_name)
        new_file_list = []
        for file_name in file_list:
            tempo = ROOT.TFile(file_name)
            dirlist = [b.GetTitle() for b in tempo.GetListOfKeys()]
            if tree_name.split("/")[0] in dirlist:
                new_file_list.append(file_name)
            tempo.Close()
        for next_file_name in new_file_list:
            tree_chain.Add(next_file_name)
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/pre_d/{spec}.root"

    if spec in mc_spec_list:
        type = "MC"
        print(type)
        ###Temp spec fix for MC####
        if "Z_m_p" in spec and "Z_m_pst" not in spec:
            temp_spec = spec.replace("Z_m_p","z")
            tree_name = f"{temp_spec}/DecayTreeTuple"
        elif "Z_z_z" in spec:
            temp_spec = spec.replace("Z_z_z","zz")
            tree_name = f"{temp_spec}/DecayTreeTuple"
        elif "P_z_p" in spec and "P_z_pst" not in spec and "_11198005" not in spec and "_11198410" not in spec:
            temp_spec = spec.replace("P_z_p","p")
            tree_name = f"{temp_spec}/DecayTreeTuple"
        elif "P_z_pst" in spec:
            temp_spec = spec.replace("P_z_pst","st")
            tree_name = f"{temp_spec}/DecayTreeTuple"
        else:
            tree_name = f"{spec}/DecayTreeTuple"

        print(tree_name)
        tree_chain = ROOT.TChain(tree_name)
        tree_name_gen = "MCDecayTreeTuple/MCDecayTree"
        tree_chain_gen = ROOT.TChain(tree_name_gen)
        for file_name in file_list:
            tree_chain.Add(file_name)
            tree_chain_gen.Add(file_name)
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/pre_d/{spec}_{trigger}.root"
        new_file_list = file_list

    if os.path.exists(outputfile):
        os.remove(outputfile)
        print("First deleting old ", outputfile)
    else:
        print("making ", outputfile)

    if "mst" in spec or "pst" in spec:
        print("creating xfd")
        new_file_name = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/st_fixes/{spec}_fix.root"
        CreateXFD(
            new_file_list,
            tree_name,
            new_file_name,
            tree_name,
        )
        new_inputfile = ROOT.TFile(new_file_name)
        tree = new_inputfile.Get(tree_name)
        rdf_base = RDF(tree)
        rdf_rec = rdf_base
    else:
        rdf_base = RDF(tree_chain)
        rdf_rec = rdf_base.Define("B_DTF_M", "B_dtf_M[0]")


    # #### Add Submass Branches ####
    # # if "norm" not in spec:
    # #     tl = ["D1", "D2", "KSTH1", "KSTH2"]
    # #     nlist = [2, 3, 4]
    # #     for n in nlist:
    # #         # Final_n_List = []
    # #         all_n_combinations = itertools.combinations(tl, n)
    # #         tupflag = 0
    # #         for tup in all_n_combinations:
    # #             tname = convertTuple(tup, spec)
    # #             rdf_temp = rdf_base.Define(tname, invmass(tup))
    # #             rdf_base = rdf_temp

    # # TRUE_PE_REC = "(pow(D1_TRUEP_E + D2_TRUEP_E + KST_TRUEP_E,2))"
    # # TRUE_PX_REC = "(pow(D1_TRUEP_X + D2_TRUEP_X + KST_TRUEP_X,2))"
    # # TRUE_PY_REC = "(pow(D1_TRUEP_Y + D2_TRUEP_Y + KST_TRUEP_Y,2))"
    # # TRUE_PZ_REC = "(pow(D1_TRUEP_Z + D2_TRUEP_Z + KST_TRUEP_Z,2))"
    # # TRUE_BM_REC = (f"sqrt({TRUE_PE_REC} - ({TRUE_PX_REC} + {TRUE_PY_REC} + {TRUE_PZ_REC}))")

    bkg_cut = "(D1_M > 0)"
    ###### Get Trigger and Loose Cuts
    if type == "MC":
        trigger_cut = fl_dict[trigger]
    if type == "DATA":
        trigger_cut = fl_dict[trigger]

    if "norm" not in spec:
        loose_cut = loose_filters
    if "norm" in spec:
        loose_cut = loose_filters_norm
    ##### Dira Cuts ###########
    if "_mst" and "_pst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX > 0)"
    if "_mst" in spec and "_pst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX_FIXED > 0 && D2_DIRA_ORIVX > 0)"
    if "_pst" in spec and "_mst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX_FIXED > 0)"
    if "_pst" in spec and "_mst" in spec:
        dira_cut = "(D1_DIRA_ORIVX_FIXED > 0 && D2_DIRA_ORIVX_FIXED > 0)"
    # #### Apply All Cuts #######
    all_cuts_string = (
        f"{bkg_cut} && {dira_cut} && {loose_cut} && {trigger_cut}"
    )
    if type == "MC":
        all_cuts_string = f"{all_cuts_string} && {build_truth_strings(spec)}"

    # print(all_cuts_string)
    # print(trigger_cut)

    rdf_final = rdf_rec.Filter(all_cuts_string)
    # rdf_final_trigger_test = rdf_final.Filter(trigger_cut)

    # base = rdf_base.Count().GetValue()
    # cut_pre_trigger = rdf_final.Count().GetValue()
    # cut_after_trigger = rdf_final_trigger_test.Count().GetValue()
    #
    # print(f"{spec}_{year} reconstructed before cuts: {base }")
    # print(f"{spec}_{year} reconstructed after cuts: {cut_pre_trigger}")
    # print(f"{spec}_{year} reconstructed after cuts + trigger condition: {cut_after_trigger}")
    # print(f"{spec}_{year} cut eff (no trigger) is: {cut_pre_trigger/base}")
    # print(f"{spec}_{year} cut eff (with trigger) is: {cut_after_trigger/base}")

    clist = rdf_final.GetColumnNames()

    print(f"Starting snapshot for {outputfile}")

    rdfsnap = rdf_final.Snapshot(f"DecayTreeTuple_{spec}", outputfile, clist, opts)
    if type == "MC":
        rdf_MC_GEN = RDF(tree_chain_gen)
        mc_clist = rdf_MC_GEN.GetColumnNames()
        rdfsnap_MC = rdf_MC_GEN.Snapshot(
            f"MCDecayTreeTuple", outputfile, mc_clist, opts
        )
    print(f"finished snapshot for {outputfile}")

# for spec in data_spec_list:
#     for year in ["2016", "2017", "2018"]:
#         file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_data/{year}/*/*/ntuple.root")
#         trigger = "ToT"
#         filter(file_list, spec, year, trigger)

for spec in mc_spec_list:
    for year in ["2016", "2017", "2018"]:
        # if spec == "02_Z_mst_p_11198005":
        #     file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/02_P_Z_m_p_p_11198005/*/*/ntuple.root"
        #     st_flag = 1
        # else:
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{spec}/*{year}*/*/ntuple.root")
        if "P_z_pst" in spec:
            file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{spec}/*{year}*/*/ntuple.root".replace("_P_z_pst", "_Z_z_z"))
        if "Z_mst_p" in spec:
            file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{spec}/*{year}*/*/ntuple.root".replace("_Z_mst_p", "_P_z_p"))
        for trigger in ["T","nTaT"]:
            filter(file_list, spec, year, trigger)
