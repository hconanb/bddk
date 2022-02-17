import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run')
from essentials import *

# mc_truth_match = "(abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313))"

inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_mc_base_test.xlsx"
outbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_eff_test.xlsx"

def build_rec_eff(self, tt_flag, year):

    eff_array = []
    bootrec_eff_array = []

    file_list = self.file_list

    for og_file_name in file_list:

        #Replace is for reused MC
        file_name = og_file_name.replace("03_Z_","02_Z_").replace("03_m_","02_p_").replace("04_m_","04_p_").replace("11_zz_","10_zz_")
        # mcfile = ROOT.TFile(f"{root_basepath}MC/{file_name}_{year}_{tt_flag}_spectrum_filtered.root")
        mcfile = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC_2021/{year}/post_d/{file_name}_{tt_flag}_postdcuts.root")
        mcfile_for_gen = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC_2021/{year}/pre_d/{file_name}_{tt_flag}.root")

        rec_tree = mcfile.Get("DecayTreeTuple")
        gen_tree = mcfile_for_gen.Get("MCDecayTreeTuple")

        rdf_rec = RDF(rec_tree)
        rdf_gen = RDF(gen_tree)

        #BOOTSTRAP for REDECAY
        if "norm" not in file_name:
            flag = file_name.split("_")[1]
            id = file_name.split("_")[2]
        else:
            flag  = file_name.split("_")[0]

        print(file_name, flag)

        rec_count = rdf_rec.Count().GetValue()
        gen_count = rdf_gen.Count().GetValue()

        print(rec_count, gen_count)

        if (os.path.isfile("TempFile.root")):
            os.remove("TempFile.root")
            print("Deleted TempFile")

        rdf_rec.Snapshot("TempTree", "TempFile.root", {"RD_org_eventNumber", "RD_org_runNumber"});
        tempfile = ROOT.TFile("TempFile.root")
        bootrec_tree = tempfile.Get("TempTree")

        bootrec_count = np.zeros(100)
        start_flag = 0

        for ev in bootrec_tree:
            # print (ev.RD_org_eventNumber, ev.RD_org_runNumber)
            if start_flag == 0:
                e_value = ev.RD_org_eventNumber
                r_value = ev.RD_org_runNumber
                count = 0
                start_flag = 1
            if e_value != ev.RD_org_eventNumber or r_value != ev.RD_org_runNumber:
                e_value = ev.RD_org_eventNumber
                r_value = ev.RD_org_runNumber
                bootrec_count += count*np.random.poisson(1,100)
                count = 0
            count += 1

        bootrec_count += count*np.random.poisson(1,100)
        print(bootrec_count.mean(), bootrec_count.std())

        bootrec_all = ufloat(bootrec_count.mean(), bootrec_count.std())
        rec_all = ufloat(rec_count, math.sqrt(rec_count))
        gen_all = ufloat(gen_count, math.sqrt(gen_count))
        bootrec_eff = bootrec_all/gen_all
        eff = rec_all/gen_all
        print(bootrec_eff, eff)

        eff_array.append(eff)
        bootrec_eff_array.append(bootrec_eff)
        print("Added ", og_file_name, " as ", file_name, year)
        tempfile.Close()

    self.df[f'rec_eff_{year}'] =  unumpy.nominal_values(eff_array)
    self.df[f'err_rec_eff_{year}'] = unumpy.std_devs(eff_array)

    self.df[f'bootrec_eff_{year}'] =  unumpy.nominal_values(bootrec_eff_array)
    self.df[f'err_bootrec_eff_{year}'] = unumpy.std_devs(bootrec_eff_array)

    return eff_array, bootrec_eff_array
def build_tot_eff(self, gen_info, rec_info, bootrec_info, year):
    tea = gen_info*rec_info
    boottea = gen_info*bootrec_info

    self.df[f'total_eff_{year}'] =  unumpy.nominal_values(tea)
    self.df[f'err_total_eff_{year}'] = unumpy.std_devs(tea)

    self.df[f'total_booteff_{year}'] =  unumpy.nominal_values(boottea)
    self.df[f'err_total_booteff_{year}'] = unumpy.std_devs(boottea)
class data_analysis():
    def __init__(self, data, tt_flag):
        self.df = pandas.DataFrame(data).dropna()

        self.gen_info_2016 = unumpy.uarray(self.df['gen_eff_2016'].values, self.df['err_gen_eff_2016'].values)
        self.gen_info_2017 = unumpy.uarray(self.df['gen_eff_2017'].values, self.df['err_gen_eff_2017'].values)
        self.gen_info_2018 = unumpy.uarray(self.df['gen_eff_2018'].values, self.df['err_gen_eff_2018'].values)

        self.front_ID = self.df['Front_ID']
        self.mc_event_ID = self.df['mc_event_ID'].astype(int).astype(str)
        self.file_list = self.front_ID + "_" + self.mc_event_ID
        print(self.file_list)

        self.rec_info_2016, self.bootrec_info_2016 = build_rec_eff(self, tt_flag, "2016")
        self.rec_info_2017, self.bootrec_info_2017 = build_rec_eff(self, tt_flag, "2017")
        self.rec_info_2018, self.bootrec_info_2018 = build_rec_eff(self, tt_flag, "2018")
        #
        # self.tea_2016, self.boottea_2016 = build_tot_eff(self, self.gen_info_2016, self.rec_info_2016, "2016")
        # self.tea_2017, self.boottea_2017 = build_tot_eff(self, self.gen_info_2017, self.rec_info_2017, "2017")
        # self.tea_2018, self.boottea_2018 = build_tot_eff(self, self.gen_info_2018, self.rec_info_2018, "2018")

        build_tot_eff(self, self.gen_info_2016, self.rec_info_2016,  self.bootrec_info_2016, "2016")
        build_tot_eff(self, self.gen_info_2017, self.rec_info_2017,  self.bootrec_info_2017, "2017")
        build_tot_eff(self, self.gen_info_2018, self.rec_info_2018,  self.bootrec_info_2018, "2018")

        print(self.df)

df_c1 = pandas.read_excel(inbook, sheet_name="MC_base")
df_c2 = df_c1.copy()
df_c3 = df_c1.copy()

# total_eff_tot = data_analysis(df_c1, "ToT")
total_eff_t = data_analysis(df_c2, "T")
total_eff_ntat = data_analysis(df_c3, "nTaT")
with pandas.ExcelWriter(outbook) as writer:
    # total_eff_tot.to_excel(writer, sheet_name='MC_eff_ToT')
    # print(f"Wrote: {outbook}, with sheet: MC_eff_TOSorTIS")
    total_eff_t.df.to_excel(writer, sheet_name='MC_eff_T')
    print(f"Wrote: {outbook}, with sheet: MC_eff_T")
    total_eff_ntat.df.to_excel(writer, sheet_name='MC_eff_nTaT')
    print(f"Wrote: {outbook}, with sheet: MC_eff_nTaT")
