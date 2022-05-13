import sys
import os

basedir = os.getcwd().split('ntuple_building')[0]
sys.path.append(basedir)

from essential_functions import *

RDF = ROOT.ROOT.RDataFrame

loose_filters_Kst_Window = "(abs(KST_M - 895) < 50)"
loose_filters_Kaon = "(D1H1_ProbNNk > 0.3 && D2H1_ProbNNk > 0.3 && KSTH1_ProbNNk > 0.3)"
loose_filters_Kaon_norm = loose_filters_Kaon.replace("KSTH1","K")

PE_mcdtt = "(pow(D1_TRUEP_E + D2_TRUEP_E + KST_TRUEP_E,2))"
PX_mcdtt = "(pow(D1_TRUEP_X + D2_TRUEP_X + KST_TRUEP_X,2))"
PY_mcdtt = "(pow(D1_TRUEP_Y + D2_TRUEP_Y + KST_TRUEP_Y,2))"
PZ_mcdtt = "(pow(D1_TRUEP_Z + D2_TRUEP_Z + KST_TRUEP_Z,2))"
M2_mcdtt = f"{PE_mcdtt} - ({PX_mcdtt} + {PY_mcdtt} + {PZ_mcdtt})"
True_BM = f"sqrt({M2_mcdtt})"

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
    tree_name = f"DecayTreeTuple"
    tree_chain = ROOT.TChain(tree_name)

    if type == "data":
        new_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/pre_d/P_z_pst.root")
        for next_file_name in new_file_list:
            tree_chain.Add(next_file_name)
    if type == "mc":
        new_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/pre_d/{spec}.root".replace("Z_z_z_", "P_z_pst_"))
        for next_file_name in new_file_list:
            tree_chain.Add(next_file_name)

    rdf_st = RDF(tree_chain)

    np_zz = rdf_zz.AsNumpy(columns=["eventNumber", "runNumber"])
    np_st = rdf_st.AsNumpy(columns=["eventNumber", "runNumber"])

    zz_en_base = np.char.mod('%d', np_zz["eventNumber"])
    zz_nCan_base = np.char.mod('%d', np_zz["runNumber"])

    st_en_base = np.char.mod('%d', np_st["eventNumber"])
    st_nCan_base = np.char.mod('%d', np_st["runNumber"])

    zz_all = np.core.defchararray.add(np.core.defchararray.add(zz_en_base, "_"),zz_nCan_base)
    st_all = np.core.defchararray.add(np.core.defchararray.add(st_en_base, "_"),st_nCan_base)

    found = [i for i in zz_all if i in st_all]

    fline = "B_DTF_M > 0"

    dict_lap = {
        "Spec": spec,
        "Year": year,
        "Type": type,
        "Candidates in Z_z_z" : rdf_zz.Count().GetValue(),
        "Candidates in P_z_pst" : rdf_st.Count().GetValue(),
        "Shared Candidates Removed from Z_z_z": len(found),
    }

    eff_df = pd.DataFrame(dict_lap , index=[0])
    eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/zz_olap/lap_{spec}_{year}.txt")

    for olap_entry in found:
        eno = olap_entry.split("_")[0]
        cno = olap_entry.split("_")[1]
        fline = f"{fline} && !(eventNumber == {eno} && runNumber == {cno})"

    return fline, found
def get_gen_eff(spec, year):
    inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/spreadsheets/b20c_gen_eff.xlsx"
    df = pd.read_excel(inbook, index_col = 0)
    nkey = spec.split("_1")[0]
    gen_info_Up = ufloat(df[f'gen_eff_{year}_Up'][nkey], df[f'err_gen_eff_{year}_Up'][nkey])
    gen_info_Down = ufloat(df[f'gen_eff_{year}_Down'][nkey], df[f'err_gen_eff_{year}_Down'][nkey])
    return ((gen_info_Up + gen_info_Down) / 2)
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
def filter(small_spec, year, type, eff_flag, snap_flag):

    spec = id_to_spec_dict[small_spec]
    if type == "mc":
        string = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_mc/{spec}/{year}/*/*/ntuple.root"
    if type == "data":
        string = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_data/{year}/*/*/ntuple.root"
    if "P_z_pst" in spec:
        string = string.replace("P_z_pst","Z_z_z")

    file_list = glob.glob(string)

    outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/pre_d/{spec}.root"

    if snap_flag:
        if os.path.exists(outputfile):
            os.remove(outputfile)
            print("First deleting old ", outputfile)
        else:
            print("making ", outputfile)
    ### Build Base Tree ####

    if type == "data":
        tree_name = f"data_{spec}_Data/DecayTreeTuple"
    if type == "mc":
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

    rdf_base = RDF(tree_chain)
    rdf_rec = rdf_base.Define("B_DTF_M", "B_dtf_M[0]")

    ### Add colums for MC only
    if type == "mc":
        if "norm" in spec:
            rdf_rec = rdf_rec.Define("B_TrueMass", True_BM.replace("KST","K"))
        else:
            if "_pst_" in spec:
                rdf_rec = rdf_rec.Define("B_TrueMass", True_BM.replace("D2_","D2st_"))
            else:
                rdf_rec = rdf_rec.Define("B_TrueMass", True_BM)
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
    if "norm" not in spec:
        b_dtf_window_cut = "B_DTF_M <= 5600 && B_DTF_M >= 4800"
    if "norm" in spec:
        b_dtf_window_cut = "B_DTF_M <= 5330 && B_DTF_M >= 5230"
    ###################

    #### Apply All Cuts#######
    if type == "data":
        if "norm" not in spec:
            all_cuts_string = f"(({dira_cut}) && ({b_dtf_window_cut}) && ({loose_filters_Kaon}))"
        else:
            all_cuts_string = f"(({dira_cut}) && ({b_dtf_window_cut}) && ({loose_filters_Kaon_norm}))"
        if spec == "Z_z_z":

            print(f"Starting snapshot for for no veto")

            rdf_nv_base = rdf_rec.Filter(all_cuts_string)

            fline, found = zz_olap(rdf_nv_base, spec, year, type)

            rdf_final = rdf_nv_base.Filter(f"{fline}")

            print(rdf_nv_base.Count().GetValue())
            print(rdf_final.Count().GetValue())
            print(len(found))

        else:
            rdf_final = rdf_rec.Filter(all_cuts_string)

    if type == "mc":
        ##### Get Generator Count ######
        tree_name_gen = "MCDecayTreeTuple/MCDecayTree"
        tree_chain_gen = ROOT.TChain(tree_name_gen)
        for file_name in file_list:
            tree_chain_gen.Add(file_name)

        print("Got Gen")

        rdf_MC_GEN = RDF(tree_chain_gen)
        rdf_GEN = rdf_numbers(rdf_MC_GEN, "GEN")
        rdf_REC = rdf_numbers(rdf_rec, "REC")

        if spec == "04_Z_z_z_11198023" or spec == "07_Z_z_z_12197045" or spec == "08_Z_z_z_12197423":

            fline, found = zz_olap(rdf_REC.rdf, spec, year, type)

            print(rdf_REC.rdf.Count().GetValue())

            rdf_REC = rdf_REC.apply_filter(f"{fline}","olap")

            print(rdf_REC.rdf.Count().GetValue())

        rdf_REC.old_stuff = rdf_GEN

        rdf_bwindow = rdf_REC.apply_filter(b_dtf_window_cut, f"b_window")

        rdf_dira = rdf_bwindow.apply_filter(dira_cut, "DIRA CUT")

        if "norm" not in spec:
            rdf_kstwindow = rdf_dira.apply_filter(loose_filters_Kst_Window, "Kst Window")
            rdf_final = rdf_kstwindow.apply_filter(loose_filters_Kaon, "Kaon ProbNN cut")
        if "norm" in spec:
            rdf_kstwindow = rdf_dira
            rdf_final = rdf_kstwindow.apply_filter(loose_filters_Kaon_norm, "Kaon ProbNN cut")

        print('All stats:')
        allCutsReport = rdf_base.Report()
        allCutsReport.Print()

    if txt_flag and type == "mc":
            gen_info = get_gen_eff(spec, year)
            eff_dict_gen_temp = {
                "Scheme ID" : id_to_scheme_dict[spec],
                "Year": year,
                "Number Accepted": f"{rdf_GEN.base_count_ufloat.n}",
                "$\epsilon_{geometrical}$": f"{gen_info*100.000:.3f}",
                }

            eff_dict_event_temp = {
                "Scheme ID" : id_to_scheme_dict[spec],
                "Year": year,
                "$\epsilon_{stripping}$": f"{rdf_REC.calc_event_eff()*100.000:.3f}",
                "$\epsilon_{offline}$": f"{rdf_final.calc_event_eff()*100.000:.3f}",
            }
            eff_dict_can_temp = {
                "Scheme ID" : id_to_scheme_dict[spec],
                "Year": year,
                "$\epsilon_{stripping}$": f"{rdf_REC.calc_can_eff()*100.000:.3f}",
                "$\epsilon_{offline}$": f"{rdf_final.calc_can_eff()*100.000:.3f}",
            }

            eff_df = pd.DataFrame(eff_dict_gen_temp , index=[0])
            eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/gen/{spec}_{year}.txt")

            eff_df = pd.DataFrame(eff_dict_can_temp , index=[0])
            eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/can/{spec}_{year}.txt")

            eff_df = pd.DataFrame(eff_dict_event_temp , index=[0])
            eff_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/event/{spec}_{year}.txt")


    if snap_flag:
        if type == "data":
            clist = rdf_final.GetColumnNames()
            clist_f = []
            for name in clist:
                if "B_dtf_" not in str(name):
                    clist_f.append(name)
            print(f"Starting snapshot for {outputfile}")
            rdfsnap = rdf_final.Snapshot(f"DecayTreeTuple", outputfile, clist_f)
            print(f"finished snapshot for {outputfile}")
        if type == "mc":
            clist = rdf_final.rdf.GetColumnNames()
            clist_f = ["B_dtf_c_nPV","nPV", "B_dtf_nPV"]
            for name in clist:
                if name not in clist_f:
                    clist_f.append(name)
            print(f"Starting snapshot for {outputfile}")
            rdfsnap = rdf_final.rdf.Snapshot(f"DecayTreeTuple", outputfile, clist_f)
            print(f"finished snapshot for {outputfile}")
    ###################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Apply PreSelections")
    parser.add_argument("--Spec_List", choices = id_to_spec_dict.keys(),  nargs="+", help = 'Spec')
    parser.add_argument('--Snap_Flag', action='store_true')
    parser.add_argument('--Txt_Flag', action='store_true')

    args = parser.parse_args()
    spec_list = args.Spec_List
    snap_flag = args.Snap_Flag
    txt_flag = args.Txt_Flag

    for small_spec in spec_list:
        if id_to_spec_dict[small_spec] == small_spec:
            type = "data"
        else:
            type = "mc"
        for year in ["2016","2017","2018"]:
            filter(small_spec, year, type, snap_flag, txt_flag)
