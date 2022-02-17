# ["MC_01_Z_m_p", ["01_Z_m_p_11198006"], ["G","DG","BG"], 5280, 10, 50],
# ["MC_02_Z_m_p", ["02_Z_m_p_11198400"], ["BG","DG","BGEP"], 5130, 10, 50],
# # # ["03_Z_m_p", ["02_Z_m_p_11198400"], ["BG"], 5130, 10, 50],
# ["MC_04_Z_m_p", ["04_Z_m_p_11198401"], ["BG","BGEP","DG","G","GEP"], 4980, 20, 60],

# from glob import glob
# import ROOT

from essentials import *

def get_files(spec, type_flag, loc_flag = "pre_d"):
    if type_flag == "Data":
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_{type_flag}/*/{loc_flag}/{spec}.root")
        print(file_list)
    if type_flag == "MC" and loc_flag == "post_d":
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_{type_flag}/*/{loc_flag}/{spec}.root")
    if type_flag == "MC" and loc_flag == "pre_d":
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_{type_flag}/*/{loc_flag}/{spec}.root")
        print (file_list)
    if type_flag == "SW_data":
        file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/sw_files/sw_{spec}.root")
    return file_list

def get_tchain(file_list, tchain_name):
    tchain = ROOT.TChain(tchain_name)
    for file_name in file_list:
        tchain.Add(file_name)
    return tchain

class mc_spectrum_class:
    def __init__(self, name, flag):
        self.spec = name
        self.files = get_files(self.spec, "MC", flag)
        self.tchain = get_tchain(self.files, f"DecayTreeTuple_{name}")

class data_spectrum_class:
    def __init__(self, name, d1_string, d2_string, d1_mass, d2_mass, d_strat):
        self.spec = name

        self.d1_string = d1_string
        self.d2_string = d2_string

        self.d1_mass = d1_mass
        self.d2_mass = d2_mass

        if self.spec == "Z_m_pst" or self.spec == "P_z_pst":
            self.d3_string = "(D^{*+} #rightarrow D^{0} #pi+)"
            self.d3_mass = 150
        if self.spec == "Z_mst_p":
            self.d3_string = "(D^{*-} #rightarrow #bar{D^{0}} #pi-)"
            self.d3_mass = 150

        self.d_strat = d_strat
        self.files = get_files(self.spec, "Data", "pre_d")
        self.tchain = get_tchain(self.files, "DecayTreeTuple_norm8")


loc_flag = "pre_d"

z_c = data_spectrum_class("norm8", "D^{-}", "D^{+}", dpmass, dpmass, "e_dg_b", loc_flag)
Z_01 = mc_spectrum_class("norm8_norm8_11198007", loc_flag)

# ["MC_01_Z_m_p", ["01_Z_m_p_11198006"], ["G","DG","BG"], 5280, 10, 50],

# sw_data_file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/sw_files/sw_Z_m_p.root")
# sw_data_tchain = ROOT.TChain("SW_tree")
# for file_name in sw_data_file_list:
#     sw_data_tchain.Add(file_name)


# sw_file = get_files("Z_m_p", "SW_data", "pre_d")
# sw_tchain = get_tchain(sw_file, "SW_tree")
# z_c.tchain.AddFriend(sw_tchain)
# z_c.tchain.Scan("SW_tree.nny_1_sw")
z_data_rdf = RDF(z_c.tchain)
# z_data_rdf = z_data_rdf.Define("sig_sw", f"SW_tree.nny_1_sw")
#

z_01_mc_rdf = RDF(Z_01.tchain)


hist_data = z_data_rdf.Histo1D((f"data_B_DTF_M", f"data_B_DTF_M", 100, 1840, 1900), "D1_M")
hist_MC = z_01_mc_rdf.Histo1D((f"MC_B_DTF_M", f"MC_B_DTF_M", 100, 1840, 1900), "D1_M")


c1 = ROOT.TCanvas("c1","c1")

hist_MC.Scale(0.5*hist_data.Integral()/hist_MC.Integral())

temp_data = hist_data.DrawCopy("")
temp_mc = hist_MC.DrawCopy("Same")

# hist_data.SetMarkerColor(ROOT.kRed)
# hist_MC.SetMarkerColor(ROOT.kBlue)
temp_data.SetLineColor(ROOT.kRed)
temp_mc.SetMarkerColor(ROOT.kBlue)

# temp_data.update()
# temp_mc.update()

legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
legend.AddEntry(temp_data,f"Data","lp")
legend.AddEntry(temp_mc,f"MC","lp")
legend.Draw()
c1.Update()

save_png(c1, f"test", f"test", rpflag = 0)

hist_data = z_data_rdf.Histo1D((f"data_B_DTF_M", f"data_B_DTF_M", 100, 1840, 1900), "D2_M")
hist_MC = z_01_mc_rdf.Histo1D((f"MC_B_DTF_M", f"MC_B_DTF_M", 100, 1840, 1900), "D2_M")


c1 = ROOT.TCanvas("c1","c1")

hist_MC.Scale(0.5*hist_data.Integral()/hist_MC.Integral())

temp_data = hist_data.DrawCopy("")
temp_mc = hist_MC.DrawCopy("Same")

# hist_data.SetMarkerColor(ROOT.kRed)
# hist_MC.SetMarkerColor(ROOT.kBlue)
temp_data.SetLineColor(ROOT.kRed)
temp_mc.SetMarkerColor(ROOT.kBlue)

# temp_data.update()
# temp_mc.update()

legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
legend.AddEntry(temp_data,f"Data","lp")
legend.AddEntry(temp_mc,f"MC","lp")
legend.Draw()
c1.Update()

save_png(c1, f"test2", f"test2", rpflag = 0)
