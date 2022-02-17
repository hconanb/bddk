import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

mc_spec_list = [
    "01_Z_m_p_11198006",
    # "02_Z_m_p_11198400",
    # "02_P_z_p_11198005",
    # "02_Z_mst_p_11198005",
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


def filter(spec):

    filein = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_MC/*/{spec}.root"
    file_list = glob.glob(filein)

    tree_name_tright = f"DecayTreeTuple_{spec}"
    tree_name_tbase = f"DecayTreeTuple_{spec}_tbase"
    tree_name_tmom = f"DecayTreeTuple_{spec}_tmom"

    tree_chain_tright = ROOT.TChain(tree_name_tright)
    tree_chain_tbase = ROOT.TChain(tree_name_tbase)
    tree_chain_tmom = ROOT.TChain(tree_name_tmom)


    for tree_chain, name in zip([tree_chain_tbase,tree_chain_tmom],["tmc1","tmc2"]):
        rdf = RDF(tree_chain)
        npy_temp = rdf.AsNumpy(columns=["eventNumber","runNumber"])

        df_temp = pandas.DataFrame(npy_temp)
        unique_events = df_temp[~df_temp.duplicated()].value_counts()
        n_unique_events = unique_events.size
        n_unique_events_err = math.sqrt(n_unique_events)
        n_unique_events_ufloat = ufloat(n_unique_events, n_unique_events_err)


        break


    # rdf_tbase = RDF(tree_chain_tbase)
    # rdf_tmom = RDF(tree_chain_tmom)
    #
    #
    # nbins = 100
    # max = 5600
    # min = 4800
    #
    # hist_tbase = rdf_tbase.Histo1D((f"BM_tbase", f"BM_tbase", nbins, min, max), f"B_DTF_M")
    # hist_tmom = rdf_tmom.Histo1D((f"BM_tmom", f"BM_tmom", nbins, min, max), f"B_DTF_M")
    #
    # bincan = ((max-min)/nbins)
    # c1 = ROOT.TCanvas("c1","c1")
    #
    # hist_tbase.GetXaxis().SetTitle("'B' Mass Wrong B or Track ID [MeV]")
    # hist_tbase.GetYaxis().SetTitle(f"Candidates / ({bincan} MeV)")
    #
    # hist_tmom.GetXaxis().SetTitle("'B' Mass Wrong D Mom [MeV]")
    # hist_tmom.GetYaxis().SetTitle(f"Candidates / ({bincan} MeV)")
    #
    # ROOT.gStyle.SetOptStat(111110)
    #
    # hist_tbase.UseCurrentStyle()
    # hist_tmom.UseCurrentStyle()
    #
    # c1.Divide(1,2)
    # c1.cd(1)
    # hist_tbase.Draw()
    # c1.cd(2)
    # hist_tmom.Draw()
    #
    # save_png(c1, f"bm_t_{spec}", f"bm_t_{spec}", rpflag=0)

for spec in mc_spec_list:
    filter(spec)
