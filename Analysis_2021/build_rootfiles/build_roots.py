import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
from createXFD import *

RDF = ROOT.ROOT.RDataFrame
# opts = ROOT.ROOT.RDF.RSnapshotOptions()
# opts.fMode = "UPDATE"

trigger_TOS = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)"
trigger_TOS_2 = "(B_L0HadronDecision_TOS == 1)"
trigger_TOS_3 = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1)"

trigger_TIS = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"

trigger_HLT1 = "(B_Hlt1TrackMVADecision_TOS == 1 || B_Hlt1TwoTrackMVADecision_TOS == 1)"
trigger_HLT2 = "(B_Hlt2Topo2BodyDecision_TOS == 1 || B_Hlt2Topo3BodyDecision_TOS == 1 || B_Hlt2Topo4BodyDecision_TOS == 1)"  # "B_Hlt2Global_TOS == 1"

loose_filters_Kst_Window = "(abs(KST_M - 895) < 50)"
loose_filters_Kaon = "(D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && KSTH1_ProbNNk > 0.3)"
loose_filters_Kaon_norm = loose_filters_Kaon.replace("KSTH1","K")

loose_filters_ghost_base = "(D1H1_ProbNNghost < 0.3 && D1H2_ProbNNghost < 0.3 && D2H1_ProbNNghost < 0.3 && D2H2_ProbNNghost < 0.3 && KSTH1_ProbNNghost < 0.3 && KSTH2_ProbNNghost < 0.3)"
loose_filters_ghost_d1m = f"{loose_filters_ghost_base} && (D1H3_ProbNNghost < 0.3)"
loose_filters_ghost_d2p = f"{loose_filters_ghost_base} && (D2H3_ProbNNghost < 0.3)"
loose_filters_ghost_d12mp = f"{loose_filters_ghost_base} && (D1H3_ProbNNghost < 0.3) && (D2H3_ProbNNghost < 0.3)"
loose_filters_ghost_norm7 = "(D1H1_ProbNNghost < 0.3 && D1H2_ProbNNghost < 0.3 && D2H1_ProbNNghost < 0.3 && D2H2_ProbNNghost < 0.3 && D2H3_ProbNNghost < 0.3 && D2H4_ProbNNghost < 0.3 && K_ProbNNghost < 0.3 && KSTH2_ProbNNghost < 0.3)"
loose_filters_ghost_norm8 = "(D1H1_ProbNNghost < 0.3 && D1H2_ProbNNghost < 0.3 && D1H3_ProbNNghost < 0.3 && D2H1_ProbNNghost < 0.3 && D2H2_ProbNNghost < 0.3 && D2H3_ProbNNghost < 0.3 && D2H4_ProbNNghost < 0.3 && K_ProbNNghost < 0.3 && KSTH2_ProbNNghost < 0.3)"

PE_mcdtt = "(pow(D1_TRUEP_E + D2_TRUEP_E + KST_TRUEP_E,2))"
PX_mcdtt = "(pow(D1_TRUEP_X + D2_TRUEP_X + KST_TRUEP_X,2))"
PY_mcdtt = "(pow(D1_TRUEP_Y + D2_TRUEP_Y + KST_TRUEP_Y,2))"
PZ_mcdtt = "(pow(D1_TRUEP_Z + D2_TRUEP_Z + KST_TRUEP_Z,2))"
M2_mcdtt = f"{PE_mcdtt} - ({PX_mcdtt} + {PY_mcdtt} + {PZ_mcdtt})"
True_BM = f"sqrt({M2_mcdtt})"

fl_dict = {
    "ToT": f"({trigger_TOS} || {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    "ToT_2":f"({trigger_TOS_2} || {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    # "nTaT": f"(!{trigger_TOS} && {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    # "T": f"({trigger_TOS}) && {trigger_HLT1} && {trigger_HLT2}",
}

data_spec_list = [
    # "Z_m_p",
    # "Zs_sm_p",
    # # # "Z_mst_p",
    # # # "Z_m_pst",
    # # # "Z_mst_pst",
    # "Z_z_z",
    # "P_z_p",
    # "M_m_z",
    # "P_z_pst",
    # # "M_mst_z",
    "norm7",
    "norm8"
    ]

mc_spec_list = [
    "01_Z_m_p_11198006",
    # "02_Z_m_p_11198400",
    # "02_P_z_p_11198005",
    # "04_Z_m_p_11198401",
    # "04_P_z_p_11198410",
    # "04_Z_z_z_11198023",
    # "04_P_z_pst_11198023",
    "05_P_z_p_12197023",
    # "06_P_z_p_12197410",
    # "07_P_z_p_12197400",
    # "07_Z_z_z_12197045",
    # "07_P_z_pst_12197045",
    # "08_P_z_p_12197401",
    # "08_Z_z_z_12197423",
    # "08_P_z_pst_12197423",
    "09_Z_z_z_11196019",
    # "10_Z_z_z_11196413",
    # "12_Z_z_z_11196414",
    "norm7_norm7_12197008",
    "norm8_norm8_11198007",
]

mc_spec_dict = {
        "01_Z_m_p_11198006" : "1",
        "02_Z_m_p_11198400" : "2a,3a",
        "02_P_z_p_11198005" : "2b,3b",
        "04_Z_m_p_11198401" : "4a",
        "04_P_z_p_11198410" : "4b",
        "04_Z_z_z_11198023" : "4c",
        "04_P_z_pst_11198023" : "4d",
        "05_P_z_p_12197023" : "5",
        "06_P_z_p_12197410" : "6",
        "07_P_z_p_12197400" : "7a",
        "07_Z_z_z_12197045": "7b",
        "07_P_z_pst_12197045": "7c",
        "08_P_z_p_12197401" : "8a",
        "08_Z_z_z_12197423" : "8b",
        "08_P_z_pst_12197423" : "8c",
        "09_Z_z_z_11196019" : "9",
        "10_Z_z_z_11196413" : "10",
        "12_Z_z_z_11196414" : "12",
        "13_Zs_sm_p_13198040" : "13",
        "14_Zs_sm_p_13198200" : "14",
        "15_Zs_sm_p_13198400" : "15",
        "16_Zs_sm_p_13198600" : "16",
        "norm7_norm7_12197008": "norm7",
        "norm8_norm8_11198007" : "norm8"
}

class rdf_numbers():
    def __init__(self, rdf_in, name):
        self.rdf = rdf_in
        self.name = name

        if name != "GEN":
            self.bm_hist =  self.rdf.Histo1D((f"{self.name}_bm", f"{self.name}_bm", 100, 4800, 5600), "B_DTF_M")
            self.bm_hist.SetStats(1)
            self.bm_hist.GetXaxis().SetTitle(f" m(B) [MeV]")
            self.bm_hist.SetTitle(name)

        base_count = self.rdf.Count().GetValue()
        base_count_err = math.sqrt(base_count)
        self.base_count_ufloat = ufloat(base_count, base_count_err)

        npy_temp = self.rdf.AsNumpy(columns=["eventNumber","runNumber"])
        df_temp = pandas.DataFrame(npy_temp)
        unique_events = df_temp[~df_temp.duplicated()].value_counts()
        n_unique_events = unique_events.size
        n_unique_events_err = math.sqrt(n_unique_events)
        self.n_unique_events_ufloat = ufloat(n_unique_events, n_unique_events_err)

        nodups = df_temp[~df_temp.duplicated(keep=False)].value_counts()
        dups = df_temp[df_temp.duplicated(keep=False)].value_counts()

        both = nodups.append(dups)
        can_per_event = both.mean()
        can_per_event_err = both.std()/math.sqrt(both.size)
        self.can_per_event_ufloat = ufloat(can_per_event, can_per_event_err)
        # if name == "GEN":
        #     self.can_per_event_ufloat = ufloat(1,0)

    # def give_prev_stuff(self, old_stuff):
    #     self.old_stuff = old_stuff

    def calc_can_eff(self):
        self.fail_count_ufloat = self.old_stuff.base_count_ufloat - self.base_count_ufloat
        next_eff = (self.base_count_ufloat/(self.base_count_ufloat + self.fail_count_ufloat))
        return next_eff

    def calc_event_eff(self):
        base_e_count_ufloat = self.old_stuff.n_unique_events_ufloat
        pass_e_count_ufloat = self.n_unique_events_ufloat
        fail_e_count_ufloat = base_e_count_ufloat - pass_e_count_ufloat
        next_eff = (pass_e_count_ufloat/(pass_e_count_ufloat + fail_e_count_ufloat))
        return next_eff

    def apply_filter(self, filter, filter_name):
        next_rdf = self.rdf.Filter(filter, filter_name)
        next_rdf_stuff = rdf_numbers(next_rdf, filter_name)
        next_rdf_stuff.old_stuff = self
        return next_rdf_stuff
def build_truth_base(b_spec, c1_spec, c2_spec, norm_flag = 0, st_flag = 0):
    truth_MC = f"abs(B_TRUEID) == {b_spec}"
    truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}"
    truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID}"
    if c1_spec == Dp_ID:
        truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID}"
    truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}"
    truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID}"
    if c2_spec == Dp_ID:
        truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID}"
    if norm_flag == 1:
        truth_MC = f"{truth_MC} && abs(D2H4_TRUEID) == {pi_ID}"
        truth_MC = f"{truth_MC} && abs(K_TRUEID) == {k_ID}"
    if norm_flag == 0:
        truth_MC = f"{truth_MC} && abs(KSTH1_TRUEID) == {k_ID}"
        truth_MC = f"{truth_MC} && abs(KSTH2_TRUEID) == {pi_ID}"
    return truth_MC
def build_truth_mom(b_spec, c1_spec, c2_spec, norm_flag = 0, st_flag = 0):
    truth_MC = f"abs(D1_TRUEID) == {c1_spec}"
    truth_MC = f"{truth_MC} && abs(D2_TRUEID) == {c2_spec}"
    if norm_flag == 0:
        truth_MC = f"{truth_MC} && abs(KST_TRUEID) == {kst0_ID}"
    return truth_MC
def zz_olap(rdf_zz, spec, year, type):

    if type == "data":
        tree_name = f"DecayTreeTuple"
        tree_chain = ROOT.TChain(tree_name)
        new_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/pre_d/P_z_pst.root")
        for next_file_name in new_file_list:
            tree_chain.Add(next_file_name)
    if type == "mc":
        tree_name = f"DecayTreeTuple"
        tree_chain = ROOT.TChain(tree_name)
        new_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/*/{spec}.root".replace("Z_z_z_", "P_z_pst_"))
        for next_file_name in new_file_list:
            tree_chain.Add(next_file_name)

    rdf_st = RDF(tree_chain)

    np_zz = rdf_zz.AsNumpy(columns=["eventNumber", "runNumber"])
    np_st = rdf_st.AsNumpy(columns=["eventNumber", "runNumber"])

    # eolap = np.intersect1d(np_zz["eventNumber"], np_st["eventNumber"])

    zz_en_base = np.char.mod('%d', np_zz["eventNumber"])
    zz_nCan_base = np.char.mod('%d', np_zz["runNumber"])

    st_en_base = np.char.mod('%d', np_st["eventNumber"])
    st_nCan_base = np.char.mod('%d', np_st["runNumber"])

    zz_all = np.core.defchararray.add(np.core.defchararray.add(zz_en_base, "_"),zz_nCan_base)
    st_all = np.core.defchararray.add(np.core.defchararray.add(st_en_base, "_"),st_nCan_base)

    found = [i for i in zz_all if i in st_all]

    fline = "B_DTF_M > 0"

    for olap_entry in found:
        eno = olap_entry.split("_")[0]
        cno = olap_entry.split("_")[1]
        fline = f"{fline} && !(eventNumber == {eno} && runNumber == {cno})"

    rdf_zz_next = rdf_zz.Filter(f"{fline}")

    print(rdf_zz.Count().GetValue())
    print(rdf_zz_next.Count().GetValue())
    print(len(found))

    dict_lap = {
        "Spec": spec,
        "Year": year,
        "Type": type,
        "Candidates in Z_z_z" : rdf_zz.Count().GetValue(),
        "Candidates in P_z_pst" : rdf_st.Count().GetValue(),
        "Shared Candidates Removed from Z_z_z": len(found),
    }

    eff_df = pandas.DataFrame(dict_lap , index=[0])
    eff_df.to_csv(f"lap_txt_files/lap_{spec}_{type}_{year}_txt")

    # clist = rdf_zz_next.GetColumnNames()
    #
    # newfilename = zzz.replace("post_veto", "nolap").replace("postveto","nolap")
    #
    # if os.path.exists(newfilename):
    #     os.remove(newfilename)
    #     print("First deleting old ", newfilename)
    # else:
    #     print("making ", newfilename)
    #
    # rdfsnap = rdf_zz_next.Snapshot("DecayTreeTuple", newfilename, clist)

    return rdf_zz_next
def build_ufloat_a_eff(base_rdf, filter_line, sr_flag = 0):

    base_count = base_rdf.Count().GetValue()
    base_count_err = math.sqrt(base_count)
    base_count_ufloat = ufloat(base_count, base_count_err)

    next_rdf = base_rdf.Filter(filter_line)

    pass_count = next_rdf.Count().GetValue()
    pass_count_err = math.sqrt(pass_count)
    pass_count_ufloat = ufloat(pass_count, pass_count_err)

    fail_count_ufloat = base_count_ufloat - pass_count_ufloat

    next_eff = (pass_count_ufloat/(pass_count_ufloat + fail_count_ufloat))*100

    # print(filter_line, pass_count_ufloat, fail_count_ufloat, next_eff)
    if sr_flag == 0:
        return next_rdf, next_eff
    if sr_flag == 1:
        inverse_next_rdf = base_rdf.Filter(f"!({filter_line})")
        return next_rdf, next_eff, inverse_next_rdf
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
def get_gen_eff(spec, year):
    inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/build_rootfiles/eff/b20c_gen_eff.xlsx"
    df = pandas.read_excel(inbook, index_col = 0)
    nkey = spec.split("_1")[0]
    gen_info_Up = ufloat(df[f'gen_eff_{year}_Up'][nkey], df[f'err_gen_eff_{year}_Up'][nkey])
    gen_info_Down = ufloat(df[f'gen_eff_{year}_Down'][nkey], df[f'err_gen_eff_{year}_Down'][nkey])
    return ((gen_info_Up + gen_info_Down) / 2)
def filter(file_list, spec, year, type, eff_flag = 1, snap_flag = 1):

    ### Build Base Tree ####
    if type == "DATA":
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
        #For now 2021
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/pre_d/{spec}.root"

        ### Prep Final Destination ###
        if os.path.exists(outputfile):
            os.remove(outputfile)
            print("First deleting old ", outputfile)
        else:
            print("making ", outputfile)


        rdf_base = RDF(tree_chain)
        rdf_rec = rdf_base.Define("B_DTF_M", "B_dtf_M[0]")
        ##############################
    if type == "MC":

        tree_name = f"{spec}/DecayTreeTuple"

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

        # print(f"Starting snapshot for outputfile_wrong_tbase and outputfile_wrong_tmom for {year}")
        ###### Get D Window Cut ########
        na = spec.split("_")
        if "norm" not in spec:
            tspec = f"{na[1]}_{na[2]}_{na[3]}"
        if "norm" in spec:
            na = spec.split("_")
            tspec = f"{na[0]}"
        dst_flag = False
        if tspec == "Z_m_p":
            d1_flag = "mp"
            d2_flag = "mp"
        if tspec == "Z_z_z":
            d1_flag = "z"
            d2_flag = "z"
        if tspec == "P_z_p":
            d1_flag = "z"
            d2_flag = "mp"
        if tspec == "M_m_z":
            d1_flag = "mp"
            d2_flag = "z"
        if tspec == "P_z_pst":
            d1_flag = "z"
            d2_flag = "z"
            dst_flag = True
        if tspec == "Zs_sm_p":
            d1_flag = "sm"
            d2_flag = "mp"
        if tspec == "norm7":
            d1_flag = "z"
            # d2_flag = "d0k3pi"
            d2_flag = "z"
        if tspec == "norm8":
            d1_flag = "mp"
            # d2_flag = "d0k3pi"
            d2_flag = "z"

        # D_Window_Cut = get_dwindow_values(tspec, d1_flag, d2_flag, dst_flag)
        # print(D_Window_Cut)
        ########################################################################

        rdf_base = RDF(tree_chain)
        rdf_rec = rdf_base.Define("B_DTF_M", "B_dtf_M[0]")
        if "norm" in spec:
            rdf_rec = rdf_rec.Define("B_TrueMass", True_BM.replace("KST","K"))
        else:
            if "_pst_" in spec:
                rdf_rec = rdf_rec.Define("B_TrueMass", True_BM.replace("D2_","D2st_"))
            else:
                rdf_rec = rdf_rec.Define("B_TrueMass", True_BM)
    ####################################################################
    ### Add Submass Branches ####
    nlist = [2, 3, 4, 5]
    if "norm" not in spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "KSTH1", "KSTH2"]
        nlist = [2, 3, 4, 5]
    if "norm7" in spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "D2H3", "D2H4","K"]
        nlist = nlist + [6]
    if "norm8" in spec:
        tl = ["D1H1", "D1H2", "D1H3", "D2H1", "D2H2", "D2H3", "D2H4","K"]
        nlist = nlist + [6, 7]
    if "Z_m_p" in spec:
        tl = tl + ["D1H3", "D2H3"]
        nlist = nlist + [6, 7]
    if "P_z_p" in spec:
        tl = tl + ["D2H3"]
        nlist = nlist + [6]
    if "M_m_z" in spec:
        tl = tl + ["D1H3"]
        nlist = nlist + [6]
    for n in nlist:
        # Final_n_List = []
        all_n_combinations = itertools.combinations(tl, n)
        tupflag = 0
        for tup in all_n_combinations:
            tname = convertTuple(tup, spec)
            rdf_temp = rdf_rec.Define(tname, invmass(tup))
            rdf_rec = rdf_temp
    #############################

    #### Create Correct DIRA Variables for soft pion reconstruction ###
    #### add accessible dtf variable ###
    #### build base rdf ###D1_TRUEP_E

    if "mst" in spec or "pst" in spec:

        fd_d2 = "pow((B_ENDVERTEX_X - D2_ENDVERTEX_X),2) + pow((B_ENDVERTEX_Y - D2_ENDVERTEX_Y),2) + pow((B_ENDVERTEX_Z - D2_ENDVERTEX_Z),2)"
        p2d2 = "sqrt(pow(D2_PX,2) + pow(D2_PY,2) + pow(D2_PZ,2))"
        bd2dot = "(D2_PX * (D2_ENDVERTEX_X - B_ENDVERTEX_X))+ (D2_PY * (D2_ENDVERTEX_Y - B_ENDVERTEX_Y))+ (D2_PZ * (D2_ENDVERTEX_Z - B_ENDVERTEX_Z))"
        D2_DIRA_ORIVX_FIXED = "bd2dot/(fd_d2*p2d2)"
        rdf_rec = rdf_rec.Define("fd_d2", fd_d2)\
                         .Define("p2d2", p2d2)\
                         .Define("bd2dot", bd2dot)\
                         .Define("D2_DIRA_ORIVX_FIXED", D2_DIRA_ORIVX_FIXED)\
                         .Define("D2stmD_M", "D2st_M - D2_M")
    if "norm" not in spec:
        rdf_rec = rdf_rec.Define("B_VtxChi2_013_fix","if (KSTH1_ID > 0) return B_VtxChi2_013; else return B_VtxChi2_023;") \
                         .Define("B_VtxChi2_023_fix","if (KSTH1_ID > 0) return B_VtxChi2_023; else return B_VtxChi2_013;") \
                         .Define("KSTPI_IPCHI2", "log(B_ENDVERTEX_CHI2 - B_VtxChi2_013_fix)") \
                         .Define("DSTM_DM",'D1H1D1H2KSTH2 - D1H1D1H2')

    ### Get Loose Cuts ###
    # if "norm" not in spec:
    #     loose_cut = loose_filters
    # if "norm" in spec:
    #     loose_cut = loose_filters_norm
    ######################

    ##### Dira Cuts ################################
    if "_mst" not in spec and "_pst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX > 0)"
    if "_mst" in spec and "_pst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX_FIXED > 0 && D2_DIRA_ORIVX > 0)"
    if "_pst" in spec and "_mst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX_FIXED > 0)"
    if "_pst" in spec and "_mst" in spec:
        dira_cut = "(D1_DIRA_ORIVX_FIXED > 0 && D2_DIRA_ORIVX_FIXED > 0)"
    #################################################

    ###### Set Windows of DTF_M
    tot_trigger_cut = fl_dict["ToT"]
    if "norm" not in spec:
        b_dtf_window_cut = "B_DTF_M <= 5600 && B_DTF_M >= 4800"
    if "norm" in spec:
        b_dtf_window_cut = "B_DTF_M <= 5330 && B_DTF_M >= 5230"
    ###################

    #####Ghost Cuts ########
    if "Z_m_p" in spec:
        ghost_cut = loose_filters_ghost_d12mp
    if "Z_z_z" in spec:
        ghost_cut = loose_filters_ghost_base
    if "P_z_p" in spec:
        ghost_cut = loose_filters_ghost_d2p
    if "M_m_z" in spec:
        ghost_cut = loose_filters_ghost_d1m
    if "norm7" in spec:
        ghost_cut = loose_filters_ghost_norm7
    if "norm8" in spec:
        ghost_cut = loose_filters_ghost_norm8

    if spec == "Z_z_z" or spec == "P_z_p":
        veto_cut = "(DSTM_DM > 150)"
    else:
        veto_cut = "(B_DTF_M > -10000)"

    #### Apply All Cuts except d_window to Data #######
    if type == "DATA":
        if "norm" not in spec:
            no_veto_string = f"({tot_trigger_cut}) && ({dira_cut}) && ({b_dtf_window_cut}) && ({loose_filters_Kst_Window}) && ({loose_filters_Kaon})" #({ghost_cut})
            all_cuts_string = f"({veto_cut} && ({no_veto_string}))"
        else:
            all_cuts_string = f"(({tot_trigger_cut}) && {veto_cut} && ({dira_cut}) && ({b_dtf_window_cut}) && ({loose_filters_Kaon_norm}))"

        if spec == "Z_z_z":
            print(f"Starting snapshot for  for no veto")
            rdf_nv_final = rdf_rec.Filter(no_veto_string)
            rdf_nv_final = zz_olap(rdf_nv_final, spec, year, type)
            clist = rdf_nv_final.GetColumnNames()
            rdfsnapnv = rdf_nv_final.Snapshot(f"DecayTreeTuple", outputfile.replace("pre_d","nveto"), clist)
            print(f"finished snapshot for for no veto")


        rdf_final = rdf_rec.Filter(all_cuts_string)
        clist = rdf_final.GetColumnNames()

        print(f"Starting snapshot for {outputfile}")
        rdfsnap = rdf_final.Snapshot(f"DecayTreeTuple", outputfile, clist)
        print(f"finished snapshot for {outputfile}")

    ###################################################
    #### Apply all cuts including d_window to MC ######
    if type == "MC":
        ##### Get Generator Count ######
        tree_name_gen = "MCDecayTreeTuple/MCDecayTree"
        tree_chain_gen = ROOT.TChain(tree_name_gen)
        for file_name in file_list:
            tree_chain_gen.Add(file_name)

        rdf_MC_GEN = RDF(tree_chain_gen)
        rdf_GEN = rdf_numbers(rdf_MC_GEN, "GEN")
        rdf_REC = rdf_numbers(rdf_rec, "REC")

        ###Possible truth matching
        ### MC specs with b0 moms

        b0_mom_list = ["01_","02_","04_","09_", "10_", "12_", "norm8"]
        bp_mom_list = ["05_","06_","07_","08_", "norm7_"]
        norm_flag = 0
        if any(ele in spec for ele in b0_mom_list):
            b_mom = B0_ID
        if any(ele in spec for ele in bp_mom_list):
            b_mom = Bp_ID
        if "_m_p_" in spec:
            d1_id = Dp_ID
            d2_id = Dp_ID
        if "_z_p_" in spec:
            d1_id = D0_ID
            d2_id = Dp_ID
        if "_z_z_" in spec or "_z_pst_" in spec or "norm7" in spec:
            d1_id = D0_ID
            d2_id = D0_ID
        if "norm8" in spec:
            d1_id = Dp_ID
            d2_id = D0_ID
        if "norm" in spec:
            norm_flag = 1

        truth_base_string = build_truth_base(b_mom, d1_id, d2_id, norm_flag)
        truth_mom_string = build_truth_mom(b_mom, d1_id, d2_id, norm_flag)

        print(f"({truth_base_string}) && ({truth_mom_string})")

        if spec == "04_Z_z_z_11198023" or spec == "07_Z_z_z_12197045" or spec == "08_Z_z_z_12197423":
            print(spec, "remove olap")
            print(type)
            rdf_REC.rdf = zz_olap(rdf_REC.rdf, spec, year, "mc")
            print("done removing olap")
        rdf_REC.old_stuff = rdf_GEN

        ###For charmless tests
        d1_mstart, d1_std, d2_mstart, d2_std = get_dwindow_values(tspec, d1_flag, d2_flag, dst_flag, rflag = "print")

        # print(d1_mstart, d2_mstart, d1_std, d2_std)

        d1_sig_max = 2*d1_std
        d2_sig_max = 2*d2_std

        d1_sb_min = 3*d1_std
        d2_sb_min = 3*d2_std

        d1_sb_max = 5*d1_std
        d2_sb_max = 5*d2_std

        d1_sb_min_line = f"(abs({d1_mstart} - D1_M) > {d1_sb_min})"
        d2_sb_min_line = f"(abs({d2_mstart} - D2_M) > {d2_sb_min})"

        d1_sb_max_line = f"(abs({d1_mstart} - D1_M) < {d1_sb_max})"
        d2_sb_max_line = f"(abs({d2_mstart} - D2_M) < {d2_sb_max})"

        d1_sb_line = f"({d1_sb_min_line} && {d1_sb_max_line})"
        d2_sb_line = f"({d2_sb_min_line} && {d2_sb_max_line})"

        d1_sig_line = f"(abs({d1_mstart} - D1_M) < {d1_sig_max})"
        d2_sig_line = f"(abs({d2_mstart} - D2_M) < {d2_sig_max})"

        sig_sig_line = f"({d1_sig_line} && {d2_sig_line})"
        sig_sb_line = f"({d1_sig_line} && {d2_sb_line})"
        sb_sig_line = f"({d1_sb_line} && {d2_sig_line})"
        sb2_line = f"({sig_sb_line} || {sb_sig_line})"

        # rdf_REC_tmc = rdf_REC.apply_filter(f"({truth_base_string}) && ({truth_mom_string})", f"full truth")

        # rdf_trigger_tos_1 = rdf_REC.apply_filter(trigger_TOS , f"trigger_line_TOS_ALL_tmc")
        # rdf_trigger_tos_2 = rdf_REC.apply_filter(trigger_TOS_2 , f"trigger_line_TOS_Had_tmc")

        rdf_bwindow = rdf_REC.apply_filter(b_dtf_window_cut, f"b_window")
        #
        rdf_dira = rdf_bwindow.apply_filter(dira_cut, "DIRA CUT")

        rdf_dwindow = rdf_dira.apply_filter(sig_sig_line, "sig_dwasjdf")

        if "norm" not in spec:
            rdf_kstwindow = rdf_dwindow.apply_filter(loose_filters_Kst_Window, "Kst Window")
        if "norm" in spec:
            rdf_kstwindow = rdf_dwindow
        #
        # # rdf_ghost = rdf_kstwindow.apply_filter(ghost_cut, "ghost Prob NN cut")
        #
        if "norm" not in spec:
            rdf_kaon = rdf_kstwindow.apply_filter(loose_filters_Kaon, "Kaon ProbNN cut")
        if "norm" in spec:
            rdf_kaon = rdf_kstwindow.apply_filter(loose_filters_Kaon_norm, "Kaon ProbNN cut")
        #
        rdf_tt = rdf_kaon.apply_filter(f"({trigger_TIS}) && ({trigger_HLT1}) &&({trigger_HLT2})")
        #
        # if "norm" not in spec:
        #     mc_all_cuts_string = f"(({tot_trigger_cut}) && {veto_cut} && ({dira_cut}) && ({b_dtf_window_cut}) && ({loose_filters_Kaon} && ({loose_filters_Kst_Window})))"
        # if "norm" in spec:
        #     mc_all_cuts_string = f"(({tot_trigger_cut}) && {veto_cut} && ({dira_cut}) && ({b_dtf_window_cut}) && ({loose_filters_Kaon_norm}))"

        rdf_trigger_tos_1 = rdf_tt.apply_filter(trigger_TOS, f"trigger_line_TOS_ALL_tmc")
        rdf_trigger_tos_2 = rdf_tt.apply_filter(trigger_TOS_2, f"trigger_line_TOS_Had_tmc")
        rdf_trigger_tos_3 = rdf_tt.apply_filter(trigger_TOS_3, f"trigger_line_TOS_HadMu_tmc")

        print(rdf_trigger_tos_1.calc_event_eff(), rdf_trigger_tos_2.calc_event_eff(), rdf_trigger_tos_3.calc_event_eff())

        # rdf_ToT_dsig = rdf_REC.apply_filter(f"({mc_all_cuts_string}) && ({sig_sig_line})", f"d_for_{spec}")
        # rdf_ToT_dsb = rdf_REC.apply_filter(f"({mc_all_cuts_string}) && ({sb2_line})", f"sb_for_{spec}")
        # rdf_ToT_r23 = rdf_REC.apply_filter(f"({mc_all_cuts_string}) && !({sig_sig_line}) && !({sb2_line})", f"r23_for_{spec}")
        # # rdf_ToT_dsigdsb = rdf_off.apply_filter(f"({D_Window_Cut}) || ({sb2_line})",f"d_and_sb_for_{spec}")
        #
        # print(rdf_off.calc_event_eff())

        print('All stats:')
        allCutsReport = rdf_base.Report()
        allCutsReport.Print()


        # eff_dict_ctest_temp = {
        #     "Scheme ID" : mc_spec_dict[spec],
        #     "Year": year,
        #     "D_WIN": f"{rdf_ToT_dsig.calc_can_eff()*100.000:.3f}",
        #     "not in D_WIN": f"{rdf_ToT_dsb.calc_can_eff()*100.000:.3f}",
        #     "D_WIN or not in DWIN": f"{rdf_ToT_dsigdsb.calc_can_eff()*100.000:.3f}",
        # }
        # eff_dict_etest_temp = {
        #     "Scheme ID" : mc_spec_dict[spec],
        #     "Year": year,
        #     "D_WIN": f"{rdf_ToT_dsig.calc_event_eff()*100.000:.3f}",
        #     "not in D_WIN": f"{rdf_ToT_dsb.calc_event_eff()*100.000:.3f}",
        #     "D_WIN or not in DWIN": f"{rdf_ToT_dsigdsb.calc_event_eff()*100.000:.3f}",
        # }
        # eff_dict_trigger = {
        #     "Scheme ID" : mc_spec_dict[spec],
        #     "Year": year,
        #     "$\epsilon_{T}$": f"{rdf_T.calc_can_eff()*100.000:.3f}",
        #     "$\epsilon_{nTaT}$":f"{rdf_nTaT.calc_can_eff()*100.000:.3f}"
        # }
        # if eff_flag == 1:
        #     gen_info = get_gen_eff(spec, year)
        #     eff_dict_gen_temp = {
        #         "Scheme ID" : mc_spec_dict[spec],
        #         "Year": year,
        #         "Number Accepted": f"{rdf_GEN.base_count_ufloat.n}",
        #         "$\epsilon_{geometrical}$": f"{gen_info*100.000:.3f}",
        #         # "$\epsilon_{stripping}$": f"{rdf_REC.base_count_ufloat/rdf_GEN.base_count_ufloat*100.000:.3f}",
        #     }
        #     eff_dict_event_temp = {
        #         "Scheme ID" : mc_spec_dict[spec],
        #         "Year": year,
        #         "$\epsilon_{stripping}$": f"{rdf_REC.calc_event_eff()*100.000:.3f}",
        #         "$\epsilon_{offline}$": f"{rdf_off.calc_event_eff()*100.000:.3f}",
        #         "$\epsilon_{D Windo w}$": f"{rdf_ToT_dsig.calc_event_eff()*100.000:.3f}",
        #     }
        #     eff_dict_can_temp = {
        #         "Scheme ID" : mc_spec_dict[spec],
        #         "Year": year,
        #         "$\epsilon_{stripping}$": f"{rdf_REC.calc_can_eff()*100.000:.3f}",
        #         "$\epsilon_{offline}$": f"{rdf_off.calc_can_eff()*100.000:.3f}",
        #         "$\epsilon_{D Window}$": f"{rdf_ToT_dsig.calc_can_eff()*100.000:.3f}",
        #     }
        #     # now = datetime.datetime.now()
        #     folder_path = '/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_efficiencies/build_root_effs'
        #     # if not os.path.exists(folder_path):
        #     #     os.makedirs(folder_path)
        #
        #     eff_df = pandas.DataFrame(eff_dict_gen_temp , index=[0])
        #     eff_df.to_csv(f"{folder_path}/{spec}_{year}_gen_eff_txt")
        #
        #     eff_df = pandas.DataFrame(eff_dict_can_temp , index=[0])
        #     eff_df.to_csv(f"{folder_path}/{spec}_{year}_can_eff_txt")
        #
        #     eff_df = pandas.DataFrame(eff_dict_event_temp , index=[0])
        #     eff_df.to_csv(f"{folder_path}/{spec}_{year}_event_eff_txt")

            # eff_df = pandas.DataFrame(eff_dict_ctest_temp , index=[0])
            # eff_df.to_csv(f"{folder_path}/charmt/{spec}_{year}_can_dwin_t_txt")
            #
            # eff_df = pandas.DataFrame(eff_dict_etest_temp , index=[0])
            # eff_df.to_csv(f"{folder_path}/charmt/{spec}_{year}_event_dwin_t_txt")

            # trig_df = pandas.DataFrame(eff_dict_trigger , index=[0])
            # trig_df.to_csv(f"eff_txt_files/{now.month}_{now.day}/{spec}_{year}_trigger_txt")

        if snap_flag == 1:

            outputfile_dsig = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/sig_d/{spec}.root"
            outputfile_dsb = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/sb_d/{spec}.root"
            outputfile_r23 = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/r23/{spec}.root"

            of_dsig = outputfile_dsig.split(f"{spec}.root")[0]
            of_dsb = outputfile_dsb.split(f"{spec}.root")[0]
            of_r23 = outputfile_r23.split(f"{spec}.root")[0]

            if not os.path.isdir(of_dsig):
                os.makedirs(of_dsig)
            if not os.path.isdir(of_dsb):
                os.makedirs(of_dsb)
            if not os.path.isdir(of_r23):
                os.makedirs(of_r23)

            if os.path.exists(outputfile_dsig):
                os.remove(outputfile_dsig)
                print("First deleting old ", outputfile_dsig)
            else:
                print("making ", outputfile_dsig)
            if os.path.exists(outputfile_dsb):
                os.remove(outputfile_dsb)
                print("First deleting old ", outputfile_dsb)
            else:
                print("making ", outputfile_dsb)
            if os.path.exists(outputfile_r23):
                os.remove(outputfile_r23)
                print("First deleting old ", outputfile_r23)
            else:
                print("making ", outputfile_r23)

            clist = rdf_off.rdf.GetColumnNames()
            print(f"Starting snapshot for {outputfile_dsig} and {outputfile_dsb} and {outputfile_r23}")

            rdfsnap_dsig = rdf_ToT_dsig.rdf.Snapshot(f"DecayTreeTuple", outputfile_dsig, clist)
            rdfsnap_dsb = rdf_ToT_dsb.rdf.Snapshot(f"DecayTreeTuple", outputfile_dsb, clist)
            rdfsnap_r23 = rdf_ToT_r23.rdf.Snapshot(f"DecayTreeTuple", outputfile_r23, clist)

            print(f"finished snapshot for {outputfile_dsig} and {outputfile_dsb}")
    ###################################################

#
# for spec in data_spec_list:
#     for year in ["2016", "2017", "2018"]:
#         file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_data/{year}/*/*/ntuple.root")
#         type = "DATA"
#         filter(file_list, spec, year, type)

for spec in mc_spec_list:
    for year in ["2016","2017","2018"]:
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_mc/{spec}/{year}/*/*/ntuple.root")
        if "P_z_pst" in spec:
            string = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_mc/{spec}/{year}/*/*/ntuple.root"
            string = string.replace("_P_z_pst","_Z_z_z")
            file_list = glob.glob(string)
        type = "MC"
        filter(file_list, spec, year, type, eff_flag = 1, snap_flag = 0)
        break
