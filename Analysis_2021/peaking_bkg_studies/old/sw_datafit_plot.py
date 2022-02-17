import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
run_name = "sw_test"

specs = ["Z_m_p"] #"Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]

rdf_list = []

ln_bin = 40

def get_mc_tree(mc_spec):
    mc_file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_MC/*/{mc_spec}*")
    mc_tree_chain = ROOT.TChain(f"DecayTreeTuple_{mc_spec}")
    for mcfile in mc_file_list:
        mc_tree_chain.Add(mcfile)
    return mc_tree_chain

def zmp_ip():
    spec = "Z_m_p"
    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/{spec}*.root")
    data_tchain = ROOT.TChain("DecayTreeTuple")
    for file_name in file_list:
        data_tchain.Add(file_name)

    sw_file = ROOT.TFile(f"root_files/sw_sw_test_{spec}.root")
    sw_tree = sw_file.Get(f"testtree")
    data_tchain.AddFriend(sw_tree, "sw_tree")

    rdf_data_all = RDF(data_tchain)
    rdf_data_all = rdf_data_all.Define("nextnny_1_sw", f"sw_tree.nny_1_sw")
    rdf_data_all = rdf_data_all.Define("nextnny_2_sw", f"sw_tree.nny_2_sw")
    rdf_data_all = rdf_data_all.Define("nextnny_3_sw", f"sw_tree.nny_3_sw")

    rdf_data_all_p0 = rdf_data_all.Filter("abs(B_DTF_M - 5280) < 50")
    rdf_data_all_p1 = rdf_data_all.Filter("abs(B_DTF_M - 5130) < 50")
    rdf_data_all_p2 = rdf_data_all.Filter("abs(B_DTF_M - 4950) < 50")

    mc_tree_01 = get_mc_tree("01_Z_m_p_11198006")
    mc_rdf_01 = RDF(mc_tree_01)

    mc_tree_02 = get_mc_tree("02_Z_m_p_11198400")
    mc_rdf_02 = RDF(mc_tree_02)

    mc_tree_04 = get_mc_tree("04_Z_m_p_11198401")
    mc_rdf_04 = RDF(mc_tree_04)



    hist_piip2_data_p0 = rdf_data_all_p0.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_data_p1 = rdf_data_all_p1.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_data_p2 = rdf_data_all_p2.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", ln_bin, -5, 5), 'KSTPI_IPCHI2')

    hist_piip2_data_sw_1 = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data_sw_1", f"KSTPI_IPCHI2_sw_1", ln_bin, -5, 5), 'KSTPI_IPCHI2', "nextnny_1_sw")
    hist_piip2_data_sw_2 = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data_sw_2", f"KSTPI_IPCHI2_sw_2", ln_bin, -5, 5), 'KSTPI_IPCHI2', "nextnny_2_sw")
    hist_piip2_data_sw_3 = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data_sw_3", f"KSTPI_IPCHI2_sw_3", ln_bin, -5, 5), 'KSTPI_IPCHI2', "nextnny_3_sw")

    hist_piip2_mc1 = mc_rdf_01.Histo1D((f"KSTPI_IPCHI2_mc1", f"KSTPI_IPCHI2_mc1", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc2 = mc_rdf_02.Histo1D((f"KSTPI_IPCHI2_mc2", f"KSTPI_IPCHI2_mc2", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc4 = mc_rdf_04.Histo1D((f"KSTPI_IPCHI2_mc4", f"KSTPI_IPCHI2_mc4", ln_bin, -5, 5), 'KSTPI_IPCHI2')

    c1 = ROOT.TCanvas("c1","c1")
    hist_piip2_data_p0.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    hist_piip2_data_p0.SetLineColor(ROOT.kRed)
    hist_piip2_data_sw_1.Scale(hist_piip2_data_p0.Integral()/hist_piip2_data_sw_1.Integral())
    hist_piip2_mc1.Scale(hist_piip2_data_p0.Integral()/hist_piip2_mc1.Integral())
    hlist = [hist_piip2_data_p0.GetPtr(), hist_piip2_data_sw_1.GetPtr(), hist_piip2_mc1.GetPtr()]
    hs = DrawStack(ROOT.gPad, hlist, drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"sw_test_01_{spec}", rpflag = 0)

    c1 = ROOT.TCanvas("c1","c1")
    hist_piip2_data_p1.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    # hist_piip2_data_sw_2.Scale(hist_piip2_data_p1.Integral()/hist_piip2_data_sw_2.Integral())
    hist_piip2_mc2.Scale(hist_piip2_data_p1.Integral()/hist_piip2_mc2.Integral())
    hlist2 = [hist_piip2_data_p1.GetPtr()  , hist_piip2_data_sw_2.GetPtr(), hist_piip2_mc2.GetPtr()]
    hs = DrawStack(ROOT.gPad, hlist2, drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"sw_test_02_{spec}", rpflag = 0)

    c1 = ROOT.TCanvas("c1","c1")
    hist_piip2_data_p2.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    hist_piip2_data_sw_3.Scale(hist_piip2_data_p2.Integral()/hist_piip2_data_sw_3.Integral())
    hist_piip2_mc4.Scale(hist_piip2_data_p2.Integral()/hist_piip2_mc4.Integral())
    hlist3 = [hist_piip2_data_p2.GetPtr(), hist_piip2_data_sw_3.GetPtr(), hist_piip2_mc4.GetPtr()]
    hs = DrawStack(ROOT.gPad, hlist3, drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"sw_test_04_{spec}", rpflag = 0)

def zzz_ip():
    spec = "Z_z_z"
    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/Z_z_z*.root")
    data_tchain = ROOT.TChain("DecayTreeTuple")
    for file_name in file_list:
        data_tchain.Add(file_name)

    sw_file = ROOT.TFile(f"root_files/sw_sw_test_{spec}.root")
    sw_tree = sw_file.Get(f"testtree")
    data_tchain.AddFriend(sw_tree, "sw_tree")

    rdf_data_all = RDF(data_tchain)

    rdf_data_all = rdf_data_all.Define("nextnny_4_sw", f"sw_tree.nny_4_sw")
    rdf_data_all = rdf_data_all.Define("nextnny_5_sw", f"sw_tree.nny_5_sw")
    rdf_data_all = rdf_data_all.Define("nextnny_6_sw", f"sw_tree.nny_6_sw")

    rdf_data_all_p0 = rdf_data_all.Filter("abs(B_DTF_M - 5280) < 50")
    rdf_data_all_p1 = rdf_data_all.Filter("abs(B_DTF_M - 5130) < 50")
    rdf_data_all_p2 = rdf_data_all.Filter("abs(B_DTF_M - 4950) < 50")
# "04_Z_z_z_11198022",
# "07_Z_z_z_12197024",
# "08_Z_z_z_12197422",
# "09_Z_z_z_11196019",
# "10_Z_z_z_11196413",
# "12_Z_z_z_11196414",

    mc_tree_04 = get_mc_tree("04_Z_z_z_11198022")
    mc_tree_07 = get_mc_tree("07_Z_z_z_12197024")
    mc_tree_08 = get_mc_tree("08_Z_z_z_12197422")
    mc_tree_09 = get_mc_tree("09_Z_z_z_11196019")
    mc_tree_10 = get_mc_tree("10_Z_z_z_11196413")
    mc_tree_12 = get_mc_tree("12_Z_z_z_11196414")

    mc_rdf_04 = RDF(mc_tree_04)
    mc_rdf_07 = RDF(mc_tree_07)
    mc_rdf_08 = RDF(mc_tree_08)
    mc_rdf_09 = RDF(mc_tree_09)
    mc_rdf_10 = RDF(mc_tree_10)
    mc_rdf_12 = RDF(mc_tree_12)

    hist_piip2_data_p0 = rdf_data_all_p0.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_data_p1 = rdf_data_all_p1.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_data_p2 = rdf_data_all_p2.Histo1D((f"KSTPI_IPCHI2_data", f"KSTPI_IPCHI2_data", ln_bin, -5, 5), 'KSTPI_IPCHI2')

    hist_piip2_data_sw_4 = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data_sw_1", f"KSTPI_IPCHI2_sw_1", ln_bin, -5, 5), 'KSTPI_IPCHI2', "nextnny_4_sw")
    hist_piip2_data_sw_5 = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data_sw_2", f"KSTPI_IPCHI2_sw_2", ln_bin, -5, 5), 'KSTPI_IPCHI2', "nextnny_5_sw")
    hist_piip2_data_sw_6 = rdf_data_all.Histo1D((f"KSTPI_IPCHI2_data_sw_3", f"KSTPI_IPCHI2_sw_3", ln_bin, -5, 5), 'KSTPI_IPCHI2', "nextnny_6_sw")

    hist_piip2_mc4 = mc_rdf_04.Histo1D((f"KSTPI_IPCHI2_mc4", f"KSTPI_IPCHI2_mc4", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc7 = mc_rdf_07.Histo1D((f"KSTPI_IPCHI2_mc7", f"KSTPI_IPCHI2_mc7", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc8 = mc_rdf_08.Histo1D((f"KSTPI_IPCHI2_mc8", f"KSTPI_IPCHI2_mc8", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc9 = mc_rdf_09.Histo1D((f"KSTPI_IPCHI2_mc9", f"KSTPI_IPCHI2_mc9", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc10 = mc_rdf_10.Histo1D((f"KSTPI_IPCHI2_mc10", f"KSTPI_IPCHI2_mc10", ln_bin, -5, 5), 'KSTPI_IPCHI2')
    hist_piip2_mc12 = mc_rdf_12.Histo1D((f"KSTPI_IPCHI2_mc12", f"KSTPI_IPCHI2_mc12", ln_bin, -5, 5), 'KSTPI_IPCHI2')

    c1 = ROOT.TCanvas("c1","c1")
    hist_piip2_data_p0.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    hist_piip2_data_p0.SetLineColor(ROOT.kRed)
    hist_piip2_data_sw_4.Scale(hist_piip2_data_p0.Integral()/hist_piip2_data_sw_4.Integral())
    hist_piip2_mc9.Scale(hist_piip2_data_p0.Integral()/hist_piip2_mc9.Integral())
    hlist = [hist_piip2_data_p0.GetPtr(), hist_piip2_data_sw_4.GetPtr(), hist_piip2_mc9.GetPtr()]
    hs = DrawStack(ROOT.gPad, hlist, drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"sw_test_04_{spec}", rpflag = 0)

    c1 = ROOT.TCanvas("c1","c1")
    hist_piip2_data_p1.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    hist_piip2_data_sw_5.Scale(hist_piip2_data_p1.Integral()/hist_piip2_data_sw_5.Integral())
    hist_piip2_mc7.Scale(hist_piip2_data_p1.Integral()/hist_piip2_mc7.Integral())
    hist_piip2_mc10.Scale(hist_piip2_data_p1.Integral()/hist_piip2_mc10.Integral())
    hlist2 = [hist_piip2_data_p1.GetPtr()  , hist_piip2_data_sw_5.GetPtr(), hist_piip2_mc7.GetPtr(), hist_piip2_mc10.GetPtr()]
    hs = DrawStack(ROOT.gPad, hlist2, drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"sw_test_05_{spec}", rpflag = 0)
    #
    c1 = ROOT.TCanvas("c1","c1")
    hist_piip2_data_p1.GetXaxis().SetTitle(f"log(B - DDK chi2)")
    hist_piip2_data_sw_6.Scale(hist_piip2_data_p1.Integral()/hist_piip2_data_sw_6.Integral())
    hist_piip2_mc4.Scale(hist_piip2_data_p1.Integral()/hist_piip2_mc4.Integral())
    hist_piip2_mc8.Scale(hist_piip2_data_p1.Integral()/hist_piip2_mc8.Integral())
    hist_piip2_mc12.Scale(hist_piip2_data_p1.Integral()/hist_piip2_mc12.Integral())
    hlist2 = [hist_piip2_data_p1.GetPtr() , hist_piip2_data_sw_6.GetPtr(), hist_piip2_mc4.GetPtr(), hist_piip2_mc8.GetPtr(),  hist_piip2_mc12.GetPtr()]
    hs = DrawStack(ROOT.gPad, hlist2, drawopts="nostack plc pmc")
    save_png(c1, "his_test", f"sw_test_06_{spec}", rpflag = 0)

# zmp_ip()
zzz_ip()
