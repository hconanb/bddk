import sys
import os

basedir = os.getcwd().split('ntuple_building')[0]
sys.path.append(basedir)

from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

def apply_bkgveto_cuts(spec, inputfile, outputfile):

    base_file = ROOT.TFile(inputfile, "READ")

    sig_tree = base_file.Get("DecayTreeTuple_SIG")
    sb_tree = base_file.Get("DecayTreeTuple_SB")

    rdf_base_signal = RDF(sig_tree)
    rdf_base_sb = RDF(sb_tree)

    if "Z_m_p" in spec:
        print("ASDF")
        rdf_post_cuts_SIG = rdf_base_signal.Filter("D1H1D2H1 > 1030 && D2H1KSTH1 > 1030")
        rdf_post_cuts_SB = rdf_base_sb.Filter("D1H1D2H1 > 1030 && D2H1KSTH1 > 1030")
    elif "Z_z_z" in spec or "P_z_p" in spec and "P_z_pst" not in spec:
        rdf_post_cuts_SIG = rdf_base_signal.Filter("DSTM_DM > 150")
        rdf_post_cuts_SB = rdf_base_sb.Filter("DSTM_DM > 150")
    else:
        rdf_post_cuts_SIG = rdf_base_signal
        rdf_post_cuts_SB = rdf_base_sb

    clist = rdf_base_signal.GetColumnNames()
    if type == "mc":
        clist_f = ["B_dtf_c_nPV","nPV", "B_dtf_nPV"]
    if type == "data":
        clist_f = ["nPV", "B_dtf_nPV"]
    for name in clist:
        if name not in clist_f:
            clist_f.append(name)

    print(rdf_post_cuts_SIG.Count().GetValue()/rdf_base_signal.Count().GetValue())

    if os.path.exists(outputfile):
        os.remove(outputfile)
        print("First deleting old ", outputfile)
    else:
        print("making ", outputfile)

    print(f"Starting snapshot for {outputfile}")
    rdf_post_cuts_SIG.Snapshot(f"DecayTreeTuple_SIG", outputfile, clist_f, opts)
    rdf_post_cuts_SB.Snapshot(f"DecayTreeTuple_SB", outputfile, clist_f, opts)
    print(f"finished snapshot for {outputfile}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Apply Bkg Vetos")
    parser.add_argument("--Spec_List", choices = id_to_spec_dict.keys(),  nargs="+", help = 'Spec')
    args = parser.parse_args()
    spec_list = args.Spec_List
    # snap_flag = args.Snap_Flag
    # txt_flag = args.Txt_Flag
    for spec in spec_list:
        if id_to_spec_dict[spec] == spec:
            type = "data"
        else:
            type = "mc"
        for year in ["2016","2017","2018"]:
            inputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_signal_and_sb/{id_to_spec_dict[spec]}.root"
            outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_veto/{id_to_spec_dict[spec]}.root"
            apply_bkgveto_cuts(id_to_spec_dict[spec], inputfile, outputfile)
