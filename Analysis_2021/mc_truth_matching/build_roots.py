import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

trigger_TOS = "(B_L0HadronDecision_TOS == 1 || B_L0MuonDecision_TOS == 1 ||  B_L0ElectronDecision_TOS == 1 ||  B_L0PhotonDecision_TOS == 1)"
trigger_TIS = "(B_L0HadronDecision_TIS == 1 || B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"
trigger_HLT1 = "(B_Hlt1TrackMVADecision_TOS == 1 || B_Hlt1TwoTrackMVADecision_TOS == 1)"
trigger_HLT2 = "(B_Hlt2Topo2BodyDecision_TOS == 1 || B_Hlt2Topo3BodyDecision_TOS == 1 || B_Hlt2Topo4BodyDecision_TOS == 1)"  # "B_Hlt2Global_TOS == 1"
loose_filters = "abs(KST_M - 895) < 50 && D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && KSTH1_ProbNNk > 0.3"
loose_filters_norm = "D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && K_ProbNNk > 0.3"

fl_dict = {
    "ToT": f"({trigger_TOS} || {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    "nTaT": f"(!{trigger_TOS} && {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
    "T": f"({trigger_TOS}) && {trigger_HLT1} && {trigger_HLT2}",
}

mc_spec_dict = {
    "01_Z_m_p_11198006" : "1",
    # "02_Z_m_p_11198400" : "2a,3a",
    # # "02_P_z_p_11198005" : "2b,3b",
    # # # # "02_Z_mst_p_11198005",
    # "04_Z_m_p_11198401" : "4a",
    # # "04_P_z_p_11198410",
    # # "04_Z_mst_p_11198410",
    # "04_Z_z_z_11198022",
    # # "04_P_z_pst_11198022",
    # "05_P_z_p_12197023" : "5",
    # # "06_P_z_p_12197410",
    # # "07_P_z_p_12197400",
    # # "07_Z_z_z_12197024",
    # # "07_P_z_pst_12197024",
    # # "08_P_z_p_12197401",
    # # "08_Z_z_z_12197422",
    # # "08_P_z_pst_12197422",
    # "09_Z_z_z_11196019" : "9",
    # "10_Z_z_z_11196413" : "10",
    # "12_Z_z_z_11196414" : "12",
    # # "13_Zs_sm_p_13198040",
    # # "14_Zs_sm_p_13198200",
    # # "15_Zs_sm_p_13198400",
    # # "16_Zs_sm_p_13198600",
    # # "norm7_norm7_12197008",
    # "norm8_norm8_11198007" : "norm8"
}

class rdf_numbers():
    def __init__(self, rdf_in, name):
        self.rdf = rdf_in
        self.name = name

        if name != "gen":
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
        if name == "gen":
            self.can_per_event_ufloat = ufloat(1,0)

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

def mult_can_checks(rdf):

    npy_temp = rdf.AsNumpy(columns=["eventNumber","runNumber"])
    df_temp = pandas.DataFrame(npy_temp)
    unique_events = df_temp[~df_temp.duplicated()].value_counts()

    n_unique_events = unique_events.size
    n_unique_events_err = math.sqrt(n_unique_events)
    n_unique_events_ufloat = ufloat(n_unique_events, n_unique_events_err)

    nodups = df_temp[~df_temp.duplicated(keep=False)].value_counts()
    dups = df_temp[df_temp.duplicated(keep=False)].value_counts()
    both = nodups.append(dups)
    can_per_event = both.mean()
    can_per_event_err = both.std()
    cen_per_event_ufloat = ufloat(can_per_event, can_per_event_err)


    return n_unique_events_ufloat, can_per_event_ufloat

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

def get_key_con(plist, gd1="_", gd2="_", gd3="_"):
    n = len(plist)
    if n == 2:
        newstring = f"{plist[0]}_MC{gd1}MOTHER_KEY == {plist[1]}_MC{gd2}MOTHER_KEY"
    if n == 3:
        newstring = f"{plist[0]}_MC{gd1}MOTHER_KEY == {plist[1]}_MC{gd2}MOTHER_KEY && {plist[0]}_MC{gd1}MOTHER_KEY == {plist[2]}_MC{gd3}MOTHER_KEY && {plist[1]}_MC{gd2}MOTHER_KEY == {plist[2]}_MC{gd3}MOTHER_KEY"
    if n == 4:
        newstring = f"{plist[0]}_MC_MOTHER_KEY == {plist[1]}_MC_MOTHER_KEY && {plist[0]}_MC_MOTHER_KEY == {plist[2]}_MC_MOTHER_KEY && {plist[0]}_MC_MOTHER_KEY == {plist[3]}_MC_MOTHER_KEY && {plist[1]}_MC_MOTHER_KEY == {plist[2]}_MC_MOTHER_KEY && {plist[1]}_MC_MOTHER_KEY == {plist[3]}_MC_MOTHER_KEY && {plist[2]}_MC_MOTHER_KEY == {plist[3]}_MC_MOTHER_KEY"
    return newstring

def build_truth_base(b_spec, c1_spec, c2_spec, s_spec, norm_flag = 0, st_flag = 0):
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

def build_truth_mom(b_spec, c1_spec, c2_spec, s_spec, norm_flag = 0, st_flag = 0):
    truth_MC = f"abs(D1_TRUEID) == {c1_spec}"
    truth_MC = f"{truth_MC} && abs(D2_TRUEID) == {c2_spec}"
    if norm_flag == 0:
        truth_MC = f"{truth_MC} && abs(KST_TRUEID) == {kst0_ID}"
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

def filter(file_list, spec, year, type):

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

        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_MC/{year}/{spec}.root"
        outputfile_wrong_tbase = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_MC/{year}/{spec}_tbase.root"

        if os.path.exists(outputfile):
            os.remove(outputfile)
            print("First deleting old ", outputfile)
        else:
            print("making ", outputfile)

        ###### Get D Window Cut ########
        na = spec.split("_")
        if "norm" not in spec:
            tspec = f"{na[1]}_{na[2]}_{na[3]}"
        if "norm" in spec:
            na = spec.split("_")
            tspec = f"{na[0]}"
        dwindow_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/build_rootfiles/d_window_root_files/d_{tspec}_mass_fits.root", "READ")
        dwindow_ws = dwindow_file.Get(f"d_{tspec}_mass_fits")
        D_Window_Cut = get_dwindow_values(dwindow_ws, tspec, rflag = "apply")
        ########################################################################
        ###### Get Truth Cut for MC ###########3
        if tspec == "Z_m_p":
            truth_string_id = build_truth_base(B0_ID, Dp_ID, Dp_ID, kst0_ID)
            truth_string_mom = build_truth_mom(B0_ID, Dp_ID, Dp_ID, kst0_ID)
        #
        # if tspec == "P_z_p":
        #     truth_string_id = build_truth_base(Bp_ID, D0_ID, Dp_ID, kst0_ID)
        #     truth_string_mom = build_truth_mom(Bp_ID, D0_ID, Dp_ID, kst0_ID)
        #
        # if tspec == "Z_z_z":
        #     truth_string_id = build_truth_base(B0_ID, D0_ID, D0_ID, kst0_ID)
        #     truth_string_mom = build_truth_mom(B0_ID, D0_ID, D0_ID, kst0_ID)
        #
        # if tspec == "norm8":
        #     truth_string_id = build_truth_base(B0_ID, Dp_ID, D0_ID, kst0_ID, norm_flag = 1)
        #     truth_string_mom = build_truth_mom(B0_ID, Dp_ID, D0_ID, kst0_ID, norm_flag = 1)

    #### Create Correct DIRA Variables for soft pion reconstruction ###
    #### add accessible dtf variable ###
    #### build base rdf ###
    if "mst" in spec or "pst" in spec:
        print("creating xfd")
        new_file_name = outputfile.replace(f"{spec}.root",f"{spec}_st_fix.root")
        CreateXFD(
            new_file_list,
            tree_name,
            new_file_name,
            tree_name,
        )
        new_inputfile = ROOT.TFile(new_file_name)
        tree = new_inputfile.Get(tree_name)
        rdf_base = RDF(tree)
    else:
        rdf_base = RDF(tree_chain)

    ####################################################################
    rdf_rec = rdf_base.Define("B_DTF_M", "B_dtf_M[0]") \
                      .Define("B_VtxChi2_013_fix","if (KSTH1_ID > 0) return B_VtxChi2_013; else return B_VtxChi2_023;") \
                      .Define("B_VtxChi2_023_fix","if (KSTH1_ID > 0) return B_VtxChi2_023; else return B_VtxChi2_013;") \
                      .Define("KSTPI_IPCHI2", "log(B_ENDVERTEX_CHI2 - B_VtxChi2_013_fix)")


    ### Add Submass Branches ####
    if "norm" not in spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "KSTH1", "KSTH2"]
        nlist = [2, 3, 4, 5]
        if "Z_m_p" in spec:
            tl = tl + ["D1H3", "D2H3"]
            nlist = nlist + [6, 7]
        if "P_z_p" in spec:
            tl = tl + ["D2H3"]
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

    ### Get Loose Cuts ###
    if "norm" not in spec:
        loose_cut = loose_filters
    if "norm" in spec:
        loose_cut = loose_filters_norm
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

    tot_trigger_cut = fl_dict["ToT"]
    b_dtf_window_cut = "B_DTF_M <= 5600 && B_DTF_M >= 4800"

    ###################
    all_cuts_string = (
        f"{loose_cut} && {tot_trigger_cut} && {dira_cut} && {b_dtf_window_cut}"
    )
    ###################################################
    #### Apply all cuts including d_window to MC ######
    if type == "MC":

        ##### Get Generator Count ######
        tree_name_gen = "MCDecayTreeTuple/MCDecayTree"
        tree_chain_gen = ROOT.TChain(tree_name_gen)
        for file_name in file_list:
            tree_chain_gen.Add(file_name)

        rdf_MC_GEN = RDF(tree_chain_gen)

        # gen_stuff = rdf_numbers(rdf_MC_GEN, "gen")
        # gen_stuff.base_count_ufloat = gen_stuff.n_unique_events_ufloat
        # gen_stuff.old_stuff = gen_stuff
        #
        # rec_stuff = rdf_numbers(rdf_rec, "rec")
        # rec_stuff.old_stuff = gen_stuff

        rdf_rec_1 = rdf_rec.Filter(all_cuts_string, f"all_before_d_{spec}")
        rdf_final_ToT = rdf_rec_1.Filter(D_Window_Cut, f"d_for_{spec}")
        rdf_final_T = rdf_final_ToT.Filter(fl_dict["T"], "T Trigger Cut")
        rdf_final_nTaT = rdf_final_ToT.Filter(fl_dict["nTaT"], "nTaT Trigger Cut")

        bm_rec = rdf_final_ToT.Histo1D((f"2016 D+D-K*0 MC: No Truth Matching", f"rec_bm", 100, 5200, 5400), "B_DTF_M")
        bm_rec.SetTitle("2016 D+D-K*0 MC: No Truth Matching")
        bm_rec.SetStats(1)
        bm_rec.GetXaxis().SetTitle(f" m_{{DTF}}(B) [MeV]")

        rdf_bad_tmc = rdf_final_ToT.Filter(f"!({truth_string_id})")

        bm_tmc = rdf_bad_tmc.Histo1D((f"2016 D+D-K*0 MC: Failed Truth Matching", f"tmc_bm", 100, 5200, 5400), "B_DTF_M")
        bm_tmc.SetTitle("2016 D+D-K*0 MC: Events Failing Truth Matching on B and D's")
        bm_tmc.SetStats(1)
        bm_tmc.GetXaxis().SetTitle(f" m_{{DTF}}(B) [MeV]")

        c1 = ROOT.TCanvas("c1","c1")
        c1.Divide(1,2)
        c1.cd(1)
        ROOT.gStyle.SetOptStat("ne")
        bm_rec.Draw()
        c1.cd(2)
        bm_tmc.Draw()
        save_png(c1, f"bm_test", f"bm_{spec}_{year}_nodira", rpflag=0)

        # for stuff in [rec_stuff, res_pass_rest, rec_tmc1_pass, rec_tmc1_fail, rec_tmc2_pass, rec_tmc2_fail]:
        #     print(f"for {stuff.name} n Candidates: {stuff.base_count_ufloat}")
        #     print(f"for {stuff.name} n unique events: {stuff.n_unique_events_ufloat}")
        #     print(f"for {stuff.name} event eff: {stuff.calc_event_eff()}")
        #     print(f"for {stuff.name} can eff: {stuff.calc_can_eff()}")
        #     print(f"for {stuff.name} ncan/event: {stuff.can_per_event_ufloat}\n")

        print('All stats:')
        allCutsReport = rdf_base.Report()
        allCutsReport.Print()

        # rdfsnap = rec_tmc1_fail.rdf.Snapshot(f"DecayTreeTuple_{spec}", "temp_fail.root")

        ##### Old Test 9/8/2021
        # tm_id_p = rec_stuff.apply_filter(truth_string_id,"tmc_id_pass")
        # tm_id_f = rec_stuff.apply_filter(f"!({truth_string_id})","tmc_id_fail")
        #
        # tm_mom_p = rec_stuff.apply_filter(truth_string_mom,"tmc_mom_pass")
        # tm_mom_f = rec_stuff.apply_filter(f"!({truth_string_mom})","tmc_mom_fail")
        #
        # tm_mom_pp = tm_id_p.apply_filter(truth_string_mom,"tmc_id_pass__tmc_mom_pass")
        # tm_mom_pf = tm_id_p.apply_filter(f"!({truth_string_mom})","tmc_id_pass__tmc_mom_fail")
        #
        # tm_mom_fp = tm_id_f.apply_filter(truth_string_mom,"tmc_id_fail__tmc_mom_pass")
        # tm_mom_ff = tm_id_f.apply_filter(f"!({truth_string_mom})","tmc_id_fail__tmc_mom_fail")
        #
        # final_stuff_list = []
        # for i in [gen_stuff, rec_stuff, tm_id_p, tm_id_f, tm_mom_p, tm_mom_f, tm_mom_pp, tm_mom_pf, tm_mom_fp, tm_mom_ff]:
        #     print(f"for {i.name} n unique events: {i.n_unique_events_ufloat}")
        #     print(f"for {i.name} event eff: {i.calc_event_eff()}")
        #     print(f"for {i.name} can eff: {i.calc_can_eff()}")
        #     print(f"for {i.name} ncan/event: {i.can_per_event_ufloat}\n")
        #
        # tm_id_p_r = tm_id_p.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_id_p.name}_passrest")
        # tm_id_f_r = tm_id_f.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_id_f.name}_passrest")
        #
        # tm_mom_p_r = tm_mom_p.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_mom_p.name}_passrest")
        # tm_mom_f_r = tm_mom_f.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_mom_f.name}_passrest")
        #
        # tm_mom_pp_r = tm_mom_pp.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_mom_pp.name}_passrest")
        # tm_mom_pf_r = tm_mom_pf.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_mom_pf.name}_passrest")
        #
        # tm_mom_fp_r = tm_mom_fp.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_mom_fp.name}_passrest")
        # tm_mom_ff_r = tm_mom_ff.apply_filter(f"{D_Window_Cut} && {dira_cut} && {loose_cut}",f"{tm_mom_ff.name}_passrest")
        #
        # for final_stuff in [tm_id_p_r, tm_id_f_r, tm_mom_p_r, tm_mom_f_r, tm_mom_pp_r, tm_mom_pf_r, tm_mom_fp_r, tm_mom_ff_r]:
        #     print(f"for {final_stuff.name} n unique events: {final_stuff.n_unique_events_ufloat}")
        #     print(f"for {final_stuff.name} event eff: {final_stuff.calc_event_eff()}")
        #     print(f"for {final_stuff.name} can eff: {final_stuff.calc_can_eff()}")
        #     print(f"for {final_stuff.name} ncan/event: {final_stuff.can_per_event_ufloat}\n")

        # ######################################
        # ##### rec -> truth_1 -> truth_2 -> loose
        #
        # rdf_next, etruth_base, rdf_tbase_fail = build_ufloat_a_eff(rdf_rec, truth_string_id, sr_flag = 1)
        # rdf_next, etruth_mom, rdf_tmom_fail = build_ufloat_a_eff(rdf_next, truth_string_mom, sr_flag = 1)
        #
        # ####################################
        #
        # rdf_rec_wrong_tbase = rdf_tbase_fail.Filter(f"{trigger_line_fwrong} && {D_Window_Cut} && {dira_cut} && {loose_cut}")
        # rdf_rec_wrong_tmom = rdf_tmom_fail.Filter(f"{trigger_line_fwrong} && {D_Window_Cut} && {dira_cut} && {loose_cut}")
        #
        # clist_wrong_tbase = rdf_rec_wrong_tbase.GetColumnNames()
        # clist_wrong_tmom = rdf_rec_wrong_tmom.GetColumnNames()
        #
        # # rdfsnap = rdf_rec_wrong_tbase.Snapshot(f"DecayTreeTuple_{spec}", outputfile_wrong_tbase, clist_wrong_tbase, opts)
        # # rdfsnap = rdf_rec_wrong_tmom.Snapshot(f"DecayTreeTuple_{spec}", outputfile_wrong_tmom, clist_wrong_tmom, opts)
        # # print(f"Done snapshot for outputfile_wrong_tbase and outputfile_wrong_tmom for {year}")
        #
        # #####################################
        #
        # rdf_next, eloose = build_ufloat_a_eff(rdf_next, loose_cut)
        # rdf_next, edira = build_ufloat_a_eff(rdf_next, dira_cut)
        # rdf_next, edwin = build_ufloat_a_eff(rdf_next, D_Window_Cut)
        # print(D_Window_Cut)
        #
        ### Apply Trigger Condtions ################
        # rdf_final_ToT, eToT = build_ufloat_a_eff(rdf_next, fl_dict["ToT"])
        # rdf_final_T, eT = build_ufloat_a_eff(rdf_final_ToT, fl_dict["T"])
        # rdf_final_nTaT, enTaT = build_ufloat_a_eff(rdf_final_ToT, fl_dict["nTaT"])
        # ############################################
        # ####### Final Eff's with everything #############
        # final_count_ToT = rdf_final_ToT.Count().GetValue()
        # final_count_ToT_err = math.sqrt(final_count_ToT)
        # final_count_ToT_ufloat = ufloat(final_count_ToT, final_count_ToT_err)
        # efinal_ToT = (final_count_ToT_ufloat/gen_count_ufloat)*100
        #
        # final_count_T = rdf_final_T.Count().GetValue()
        # final_count_T_err = math.sqrt(final_count_T)
        # final_count_T_ufloat = ufloat(final_count_T, final_count_T_err)
        # efinal_T = (final_count_T_ufloat/gen_count_ufloat)*100
        #
        # final_count_nTaT = rdf_final_nTaT.Count().GetValue()
        # final_count_nTaT_err = math.sqrt(final_count_nTaT)
        # final_count_nTaT_ufloat = ufloat(final_count_nTaT, final_count_nTaT_err)
        # efinal_nTaT = (final_count_nTaT_ufloat/gen_count_ufloat)*100
        # ###############################################
        #
        # ##### Get Generator Eff ######
        # geo_eff_book = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_gen_eff.xlsx"
        # geo_eff_df = pandas.read_excel(geo_eff_book)
        #
        # geo_slice = geo_eff_df.loc[geo_eff_df['eff_script_id'] == spec]
        # geo_info_year_Up = ufloat(geo_slice[f'gen_eff_{year}_Up'].values, geo_slice[f'err_gen_eff_{year}_Up'].values)
        # geo_info_year_Down = ufloat(geo_slice[f'gen_eff_{year}_Down'].values, geo_slice[f'err_gen_eff_{year}_Down'].values)
        #
        # egeo = ((geo_info_year_Up + geo_info_year_Down)/2)*100
        #
        # ######## Dict for final print
        # mcdf_dict_temp = {'Scheme ID' : mc_spec_dict[spec],
        #                   'Year': year,
        #                   'Number Accepted' : gen_count_ufloat.n,
        #                   '$\epsilon_{geometrical}$': f'${egeo:.3fL}$',
        #                   '$\epsilon_{reconstruction}$': f'${ereco:.3fL}$',
        #                   '$\epsilon_{truth_base}$': f'${etruth_base:.3fL}$',
        #                   '$\epsilon_{truth_mom}$': f'${etruth_mom:.3fL}$',
        #                   '$\epsilon_{loose}$': f'${edira:.3fL}$',
        #                   '$\epsilon_{dira}$': f'${eloose:.3fL}$',
        #                   '$\epsilon_{dwin}$': f'${edwin:.3fL}$',
        #                   '$\epsilon_{ToT}$': f'${eToT:.3fL}$',
        #                   '$\epsilon_{T}$': f'${eT:.3fL}$',
        #                   '$\epsilon_{nTaT}$': f'${enTaT:.3fL}$',
        #                   '$\epsilon_{final_ToT}$': f'${efinal_ToT:.3fL}$',
        #                   '$\epsilon_{final_T}$': f'${efinal_T:.3fL}$',
        #                   '$\epsilon_{final_nTaT}$': f'${efinal_nTaT:.3fL}$',
        #                   }
        #
        # ######## Dict for final spread
        # mcdf_dict_spread = {'Scheme ID' : mc_spec_dict[spec],
        #                   'Year': year,
        #                   'epsilon_geometrical': egeo,
        #                   'epsilon_final_ToT': efinal_ToT,
        #                   'epsilon_final_T': efinal_T,
        #                   'epsilon_final_nTaT': efinal_nTaT,
        #                   }
        #
        # pprint.pprint(mcdf_dict_temp,  sort_dicts=False)
        #
        # ############# Get event eff
        # reco_event_eff, ncan_reco = mult_can_checks(rdf_rec)
        # base_tbase_event_eff, ncan_base_tbase = mult_can_checks(rdf_tbase_fail)
        # base_tmom_event_eff, ncan_base_tmom = mult_can_checks(rdf_tmom_fail)
        # final_tot_event_eff, ncan_final_tot = mult_can_checks(rdf_final_ToT)
        # final_tbase_event_eff, ncan_final_tbase = mult_can_checks(rdf_rec_wrong_tbase)
        # final_tmom_event_eff, ncan_final_tmom = mult_can_checks(rdf_rec_wrong_tmom)
        #
        # print(f"reco event eff: {reco_event_eff} and ncan/event: {ncan_reco}")
        # print(f"fail first tmc event eff, no cuts: {base_tbase_event_eff} and ncan/event: {ncan_base_tbase}")
        # print(f"fail second tmc event eff, no cuts: {base_tmom_event_eff} and ncan/event: {ncan_base_tmom}")
        #
        # print(f"final good event eff: {final_tot_event_eff} and ncan/event: {ncan_final_tot}")
        # print(f"final fail tmc1 event eff: {final_tbase_event_eff} and ncan/event: {ncan_final_tbase}")
        # print(f"final fail tmc2 event eff: {final_tmom_event_eff} and ncan/event: {ncan_final_tmom}")

        # mcdf = pandas.DataFrame(mcdf_dict_temp, index=[0])
        #
        # outefftxt = f"eff_txt_files/{year}_{spec}"
        # outefftxt_test = f"eff_txt_files/{year}_{spec}_test"
        # mcdf.to_csv(outefftxt)
        #
        # mcdf = pandas.DataFrame(mcdf_dict_spread, index=[0])
        # mcdf.to_csv(outefftxt_test)
        #

    ###################################################

# for spec in data_spec_list:
#     for year in ["2016", "2017", "2018"]:
#         file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_data/{year}/*/*/ntuple.root")
#         type = "DATA"
#         filter(file_list, spec, year, type)

for key in mc_spec_dict:
    spec = key
    for year in ["2016", "2017", "2018"]:
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_MC/{spec}/{year}/*/*/ntuple.root")
        if "P_z_pst" in spec:
            string = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_MC/{spec}/{year}/*/*/ntuple.root"
            string = string.replace("_P_z_pst","_Z_z_z")
            file_list = glob.glob(string)
        # if "Z_mst_p" in spec:
        #     file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_MC/{spec}/{year}/*/ntuple.root".replace("_Z_mst_p", "_P_z_p"))
        type = "MC"
        # print(file_list)
        filter(file_list, spec, year, type)
