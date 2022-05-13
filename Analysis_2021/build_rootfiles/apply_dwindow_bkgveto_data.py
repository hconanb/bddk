import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

data_spec_list = [
    # "Z_m_p",
    # "Z_z_z",
    # "P_z_p",
    # "M_m_z",
    # "P_z_pst",
    # "norm7",
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
            # d2_flag = "d0k3pi"
            d2_flag = "z"
        if spec == "norm8":
            d1_flag = "mp"
            # d2_flag = "d0k3pi"
            d2_flag = "z"


        d1_mstart, d1_std, d2_mstart, d2_std = get_dwindow_values(spec, d1_flag, d2_flag, dst_flag, rflag = "print")

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
        print(sig_sig_line)
        rdf_ToT_dsig = rdf_base.Filter(sig_sig_line, f"d_for_{spec}")
        rdf_ToT_dsb = rdf_base.Filter(sb2_line, f"sb_for_{spec}")

        # dwindow_cut = get_dwindow_values(spec, d1_flag, d2_flag, dst_flag)

        outputfile_dsig = outputfile
        outputfile_dsb = outputfile.replace("post_d","sb_d")

        of_dsig = outputfile_dsig.split(f"{spec}.root")[0]
        of_dsb = outputfile_dsb.split(f"{spec}.root")[0]

        if not os.path.isdir(of_dsig):
            os.makedirs(of_dsig)
        if not os.path.isdir(of_dsb):
            os.makedirs(of_dsb)

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

        clist = rdf_base.GetColumnNames()
        print(f"Starting snapshot for {outputfile_dsig} and {outputfile_dsb}")

        rdfsnap_dsig = rdf_ToT_dsig.Snapshot(f"DecayTreeTuple", outputfile_dsig, clist)
        rdfsnap_dsb = rdf_ToT_dsb.Snapshot(f"DecayTreeTuple", outputfile_dsb, clist)
        print(f"finished snapshot for {outputfile_dsig} and {outputfile_dsb}")

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
