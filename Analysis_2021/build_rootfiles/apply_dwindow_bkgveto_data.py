import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

data_spec_list = [
#     # "Z_m_p",
#     # "Zs_sm_p",
#     "Z_z_z",
#     "P_z_p",
#     # "M_m_z",
#     # "P_z_pst",
    "norm7",
    "norm8"
    ]
mc_spec_list = [
    # "01_Z_m_p_11198006",
    # "02_Z_m_p_11198400",
    # "02_P_z_p_11198005",
    # # # # "04_Z_m_p_11198401",
    # "04_P_z_p_11198410",
    # "04_Z_z_z_11198023",
    # # # "04_P_z_pst_11198023",
    # "05_P_z_p_12197023",
    # "06_P_z_p_12197410",
    # "07_P_z_p_12197400",
    # "07_Z_z_z_12197045",
    # # # # "07_P_z_pst_12197045",
    # "08_P_z_p_12197401",
    # "08_Z_z_z_12197423",
    # # # # "08_P_z_pst_12197423",
    # "09_Z_z_z_11196019",
    # "10_Z_z_z_11196413",
    # "12_Z_z_z_11196414",
    # # # "13_Zs_sm_p_13198040",
    # # # "14_Zs_sm_p_13198200",
    # # # "15_Zs_sm_p_13198400",
    # # # "16_Zs_sm_p_13198600",
    # "norm7_norm7_12197008",
    # "norm8_norm8_11198007",
]

def apply_dwindow_cuts(spec, inputfile, outputfile):

        base_file = ROOT.TFile(inputfile, "READ")
        base_tree = base_file.Get(f"DecayTreeTuple")
        rdf_base = RDF(base_tree)
        dst_flag = False
        if spec == "Z_m_p":
            d1_flag = "mp"
            d2_flag = "mp"
        if spec == "Z_z_z":
            d1_flag = "z"
            d2_flag = "z"
        if spec == "P_z_p":
            d1_flag = "z"
            d2_flag = "mp"
        if spec == "M_m_z":
            d1_flag = "mp"
            d2_flag = "z"
        if spec == "P_z_pst":
            d1_flag = "z"
            d2_flag = "z"
            dst_flag = True
        if spec == "Zs_sm_p":
            d1_flag = "sm"
            d2_flag = "mp"
        if spec == "norm7":
            d1_flag = "z"
            d2_flag = "d0k3pi"
        if spec == "norm8":
            d1_flag = "mp"
            d2_flag = "d0k3pi"

        dwindow_cut = get_dwindow_values(spec, d1_flag, d2_flag, dst_flag)

        print(dwindow_cut)

        rdf_post_dcuts = rdf_base.Filter(dwindow_cut)

        clist = rdf_post_dcuts.GetColumnNames()

        if os.path.exists(outputfile):
            os.remove(outputfile)
            print("First deleting old ", outputfile)
        else:
            print("making ", outputfile)

        rdfsnap = rdf_post_dcuts.Snapshot(f"DecayTreeTuple", outputfile, clist)

def apply_bkgveto_cuts(spec, inputfile, outputfile, flag = "None"):

    base_file = ROOT.TFile(inputfile, "READ")
    if spec in mc_spec_list:
        base_tree = base_file.Get(f"DecayTreeTuple")
    else:
        base_tree = base_file.Get(f"DecayTreeTuple")
    rdf_base = RDF(base_tree)

    if "Z_z_z" in spec:
        rdf_base = rdf_base.Define("DSTM_DM",'D1H1D1H2KSTH2 - D1H1D1H2')
        rdf_post_cuts = rdf_base.Filter("DSTM_DM > 150")

    if "P_z_p" in spec:
        rdf_base = rdf_base.Define("DSTM_DM",'D1H1D1H2KSTH2 - D1H1D1H2')
        rdf_post_cuts = rdf_base.Filter("DSTM_DM > 150")

    clist = rdf_post_cuts.GetColumnNames()

    print(rdf_post_cuts.Count().GetValue()/rdf_base.Count().GetValue())

    if os.path.exists(outputfile):
        os.remove(outputfile)
        print("First deleting old ", outputfile)
    else:
        print("making ", outputfile)

    rdfsnap = rdf_post_cuts.Snapshot(f"DecayTreeTuple", outputfile, clist)


for spec in data_spec_list:
    for year in ["2016","2017","2018"]:
        inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/pre_d/{spec}.root"
        outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/{year}/post_d/{spec}.root"
        apply_dwindow_cuts(spec, inputfile, outputfile)

# for spec in ["Z_z_z","P_z_p"]:
#     for year in ["2016","2017","2018"]:
#         inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/{year}/post_d/{spec}_postdcuts.root"
#         outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/{year}/post_veto/{spec}.root"
#         apply_bkgveto_cuts(spec, inputfile, outputfile)

# for spec in mc_spec_list:
#     for year in ["2016","2017","2018"]:
#         inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/post_d/{spec}.root"
#         outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/{year}/post_veto/{spec}.root"
#         apply_bkgveto_cuts(spec, inputfile, outputfile)
