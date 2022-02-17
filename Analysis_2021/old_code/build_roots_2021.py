import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
from createXFD import *

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

data_spec_list = [
    "Z_m_p",
    # "Zs_sm_p",
    # "Z_mst_p",
    # "Z_m_pst",
    # "Z_mst_pst",
    # "Z_z_z",
    # "P_z_p",
    # "M_m_z",
    # "P_z_pst",
    # "M_mst_z",
    # "norm7",
    # "norm8"
    ]

mc_spec_list = [
    # "01_Z_m_p_11198006",
    # "02_Z_m_p_11198400",
    # # "02_P_z_p_11198005",
    # # "02_Z_mst_p_11198005",
    # "04_Z_m_p_11198401",
    # "04_P_z_p_11198410",
    # "04_Z_mst_p_11198410",
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
    # "13_Zs_sm_p_13198040",
    # "14_Zs_sm_p_13198200",
    # "15_Zs_sm_p_13198400",
    # "16_Zs_sm_p_13198600",
    # "norm7_norm7_12197008",
    # "norm8_norm8_11198007",
]

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

def filter(file_list, spec, year, trigger, type):

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
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/pre_d/{spec}.root"

    # if spec in mc_spec_list:
    #     type = "MC"
    #     ###Temp spec fix for MC####
    #     if "Z_m_p" in spec and "Z_m_pst" not in spec:
    #         temp_spec = spec.replace("Z_m_p","z")
    #         tree_name = f"{temp_spec}/DecayTreeTuple"
    #     elif "Z_z_z" in spec:
    #         temp_spec = spec.replace("Z_z_z","zz")
    #         tree_name = f"{temp_spec}/DecayTreeTuple"
    #     elif "P_z_p" in spec and "P_z_pst" not in spec and "_11198005" not in spec and "_11198410" not in spec:
    #         temp_spec = spec.replace("P_z_p","p")
    #         tree_name = f"{temp_spec}/DecayTreeTuple"
    #     elif "P_z_pst" in spec:
    #         temp_spec = spec.replace("P_z_pst","st")
    #         tree_name = f"{temp_spec}/DecayTreeTuple"
    #     elif "Zs_sm_p" in spec:
    #         temp_spec = spec.replace("Zs_sm_p","s")
    #         tree_name = f"{temp_spec}/DecayTreeTuple"
    #     else:
    #         tree_name = f"{spec}/DecayTreeTuple"
    #
    #     print(tree_name)
    #     tree_chain = ROOT.TChain(tree_name)
    #     new_file_list = []
    #     for file_name in file_list:
    #         tempo = ROOT.TFile(file_name)
    #         dirlist = [b.GetTitle() for b in tempo.GetListOfKeys()]
    #         if tree_name.split("/")[0] in dirlist:
    #             new_file_list.append(file_name)
    #         tempo.Close()
    #     for next_file_name in new_file_list:
    #         tree_chain.Add(next_file_name)
    #     tree_name_gen = "MCDecayTreeTuple/MCDecayTree"
    #     tree_chain_gen = ROOT.TChain(tree_name_gen)
    #     for file_name in file_list:
    #         tree_chain_gen.Add(file_name)
    #     outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/pre_d/{spec}_{trigger}.root"

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

    # if spec == "Z_z_z":
    # else:
    #     bkg_cut = "(D1_M > 0)"
    ###### Get Trigger and Loose Cuts
    # if type == "MC":
    #     trigger_cut = fl_dict[trigger]
    if type == "DATA":
        trigger_cut = fl_dict[trigger]

    if "norm" not in spec:
        loose_cut = loose_filters
    if "norm" in spec:
        loose_cut = loose_filters_norm
    ##### Dira Cuts ###########

    if "_mst" not in spec and "_pst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX > 0)"

    if "_mst" in spec and "_pst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX_FIXED > 0 && D2_DIRA_ORIVX > 0)"

    if "_pst" in spec and "_mst" not in spec:
        dira_cut = "(D1_DIRA_ORIVX > 0 && D2_DIRA_ORIVX_FIXED > 0)"

    if "_pst" in spec and "_mst" in spec:
        dira_cut = "(D1_DIRA_ORIVX_FIXED > 0 && D2_DIRA_ORIVX_FIXED > 0)"
    # #### Apply All Cuts #######
    all_cuts_string = (
        f"{dira_cut} && {loose_cut} && {trigger_cut}"
        # f"{bkg_cut} && {loose_cut} && {trigger_cut}"
    )

    # if type == "MC":
    #     all_cuts_string = f"{all_cuts_string} && {build_truth_strings(spec)}"


    rdf_final = rdf_rec.Filter(all_cuts_string)

    # base = rdf_base.Count().GetValue()
    # cut_pre_trigger = rdf_final.Count().GetValue()
    # cut_after_trigger = rdf_final_trigger_test.Count().GetValue()

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


for spec in data_spec_list:
    for year in ["2016", "2017", "2018"]:
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_data/{year}/*/*/ntuple.root")
        trigger = "ToT"
        type = "DATA"
        filter(file_list, spec, year, trigger, type)

# for spec in mc_spec_list:
#     for year in ["2016", "2017", "2018"]:
#         # if spec == "02_Z_mst_p_11198005":
#         #     file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/02_P_Z_m_p_p_11198005/*/*/ntuple.root"
#         #     st_flag = 1
#         # else:
#         file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{spec}/*{year}*/*/ntuple.root")
#         if "P_z_pst" in spec:
#             file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{spec}/*{year}*/*/ntuple.root".replace("_P_z_pst", "_Z_z_z"))
#         if "Z_mst_p" in spec:
#             file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc/{spec}/*{year}*/*/ntuple.root".replace("_Z_mst_p", "_P_z_p"))
#         for trigger in ["T","nTaT"]:
            # type = "MC"
            # filter(file_list, spec, year, trigger, type)
