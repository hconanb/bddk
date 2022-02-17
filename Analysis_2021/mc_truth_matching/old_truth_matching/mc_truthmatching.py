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
