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

mc_spec_list = [
    "04_zz_11198022",
    "04_st_11198022",
]


def filter(name_zz, name_st, type_flag, year, trigger_flag):
    # if type_flag == "DATA":
    #     spec = name
    #     file = ROOT.TFile(
    #         f"/mnt/c/Users/Harris/Desktop/rootfiles/data_run2/{spec}_spectrum.root"
    #     )
    #     tree = file.Get("DecayTreeTuple")
    #     rdf_data_base = RDF(tree)

    if type_flag == "MC":

        file_name_zz = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_base_mc_root/MC/{name_zz}_{year}_{trigger_flag}"
        file_name_st = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC/{name_st}_{year}_{trigger_flag}"

        file_zz = ROOT.TFile(file_name_zz)
        tree_zz = file_zz.Get(f"{name_zz}/DecayTreeTuple")

        file_st = ROOT.TFile(file_name_st)
        tree_st = file_st.Get(f"{name_st}/DecayTreeTuple")

        rdf_zz = RDF(tree_zz)
        rdf_st = RDF(tree_st)

        np_zz = rdf_zz.AsNumpy(columns=["eventNumber","nCandidate"])
        np_st = rdf_st.AsNumpy(columns=["eventNumber","nCandidate"])

        zz_en_base = np.char.mod('%d', np_zz["eventNumber"])
        zz_nCan_base = np.char.mod('%d', np_zz["nCandidate"])

        st_en_base = np.char.mod('%d', np_st["eventNumber"])
        st_nCan_base = np.char.mod('%d', np_st["nCandidate"])

        zz_all = np.core.defchararray.add(zz_en_base, zz_nCan_base)
        st_all = np.core.defchararray.add(st_en_base, st_nCan_base)

        print(np.intersect1d(zz_all,st_all))



        # print(np_zz['eventNumber'][zz_ind], np_st['eventNumber'][st_ind])


filter(mc_spec_list[0], mc_spec_list[1], "MC", "2016", "T")
