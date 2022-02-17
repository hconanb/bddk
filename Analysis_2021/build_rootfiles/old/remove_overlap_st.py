import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
from createXFD import *

def zz_olapfit(zzz, pzpst, year, type):

    zz_file = ROOT.TFile(zzz)
    st_file = ROOT.TFile(pzpst)

    tree_zz = zz_file.Get(f"DecayTreeTuple")


    tree_st = st_file.Get(f"DecayTreeTuple")

    rdf_zz = RDF(tree_zz)
    rdf_st = RDF(tree_st)

    np_zz = rdf_zz.AsNumpy(columns=["eventNumber", "runNumber"])
    np_st = rdf_st.AsNumpy(columns=["eventNumber", "runNumber"])

    eolap = np.intersect1d(np_zz["eventNumber"], np_st["eventNumber"])

    print(len(np_zz["eventNumber"]),len(np_st["eventNumber"]),len(eolap))


    # print(len(np_zz["eventNumber"]),len(np_st["eventNumber"]))
    #
    #

    #
    zz_en_base = np.char.mod('%d', np_zz["eventNumber"])
    zz_nCan_base = np.char.mod('%d', np_zz["runNumber"])

    st_en_base = np.char.mod('%d', np_st["eventNumber"])
    st_nCan_base = np.char.mod('%d', np_st["runNumber"])
    zz_all = np.core.defchararray.add(np.core.defchararray.add(zz_en_base, "_"),zz_nCan_base)
    st_all = np.core.defchararray.add(np.core.defchararray.add(st_en_base, "_"),st_nCan_base)

    found = [i for i in zz_all if i in st_all]
    #
    # print(found)
    #
    fline = "B_DTF_M > 1000"
    #
    for olap_entry in found:
        eno = olap_entry.split("_")[0]
    #     runnum = olap_entry.split("_")[1]
        fline = f"{fline} && !(eventNumber == {eno})"
    #     # if eno == "1042012908" or eno == 1042012908:
    #     #     break
         # if old_evvent.eventNumer != eno or (old_event.eventNumber = eno and old_event.nCandidate != ncan):
    # print(fline)

    rdf_zz_next = rdf_zz.Filter(f"{fline}")

    print(rdf_zz.Count().GetValue())
    print(rdf_zz_next.Count().GetValue())
    print(len(found))
    dict_lap = {
        "Spec": type,
        "Year": year,
        "Candidates in Z_z_z" : rdf_zz.Count().GetValue(),
        "Candidates in P_z_pst": rdf_st.Count().GetValue(),
        "Number of Duplicate Candidates": len(found),
    }

    eff_df = pandas.DataFrame(dict_lap , index=[0])
    eff_df.to_csv(f"lap_txt_files/lap_{type}_{year}_txt")

    clist = rdf_zz_next.GetColumnNames()

    newfilename = zzz.replace("post_veto", "nolap").replace("postveto","nolap")

    if os.path.exists(newfilename):
        os.remove(newfilename)
        print("First deleting old ", newfilename)
    else:
        print("making ", newfilename)
    rdfsnap = rdf_zz_next.Snapshot("DecayTreeTuple", newfilename, clist)

# fl_zzz = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_veto/Z_z_z_postveto.root")
# fl_pzpst  = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/P_z_pst_postdcuts.root")
#
year_list = ["2016","2017","2018"]

# for zzz, pzpst, year in zip(fl_zzz, fl_pzpst, year_list):
#     type = "Data"
#     zz_olapfit(zzz, pzpst, year, type)

for spec in ["04_P_z_pst_11198023","07_P_z_pst_12197045","08_P_z_pst_12197423"]:

    temp = spec.replace("P_z_pst","Z_z_z")
    fl_zzz_mc = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/*/post_veto/{temp}.root")
    fl_pzpst_mc  = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/*/post_d/{spec}.root")

    for zz, pzpst, year in zip(fl_zzz_mc, fl_pzpst_mc, year_list):
        zz_olapfit(zz, pzpst, year, spec)
