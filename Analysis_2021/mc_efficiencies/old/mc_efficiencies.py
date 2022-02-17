import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021')
from essentials import *

inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_mc_base.xlsx"

def build_rec_eff(self, tt_flag, year, dira_cut = 0):

    eff_array = []
    bootrec_eff_array = []

    file_list = self.file_list

    for og_file_name in file_list:

        #Replace is for reused MC
        file_name = og_file_name.replace("03_Z_m_p","02_Z_m_p").replace("03_M_m_z","02_P_z_p").replace("04_M_m_z","04_P_z_p").replace("11_zz_","10_zz_")
        mcfile = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC_2021/{year}/post_d/{file_name}_{tt_flag}_postdcuts.root")
        mcfile_for_gen = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC_2021/{year}/pre_d/{file_name}_{tt_flag}.root")

        rec_tree = mcfile.Get("DecayTreeTuple")
        gen_tree = mcfile_for_gen.Get("MCDecayTreeTuple")

        rdf_rec = RDF(rec_tree)
        rdf_gen = RDF(gen_tree)

        ###Temp line for dira test
        rdf_rec = rdf_rec.Filter(f"(D1_DIRA_ORIVX*D1_FDCHI2_ORIVX > {dira_cut}) && (D2_DIRA_ORIVX*D2_FDCHI2_ORIVX > {dira_cut})")

        #BOOTSTRAP for REDECAY
        if "norm" not in file_name:
            flag = file_name.split("_")[1]
            id = file_name.split("_")[2]
        else:
            flag  = file_name.split("_")[0]

        rec_count = rdf_rec.Count().GetValue()
        gen_count = rdf_gen.Count().GetValue()

        if (os.path.isfile("TempFile.root")):
            os.remove("TempFile.root")
            print("Deleted TempFile")

        rdf_rec.Snapshot("TempTree", "TempFile.root", {"RD_org_eventNumber", "RD_org_runNumber"});
        tempfile = ROOT.TFile("TempFile.root")
        bootrec_tree = tempfile.Get("TempTree")

        bootrec_count = np.zeros(100)
        start_flag = 0

        for ev in bootrec_tree:
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

        bootrec_all = ufloat(bootrec_count.mean(), bootrec_count.std())
        rec_all = ufloat(rec_count, math.sqrt(rec_count))
        gen_all = ufloat(gen_count, math.sqrt(gen_count))
        bootrec_eff = bootrec_all/gen_all
        eff = rec_all/gen_all

        eff_array.append(eff)
        bootrec_eff_array.append(bootrec_eff)
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
    def __init__(self, data, tt_flag, dira_cut):
        self.df = pandas.DataFrame(data).dropna()
        # self.df = self.df[~self.df.Front_ID.str.contains("_Zs_")]
        self.df = self.df[self.df.Front_ID.str.contains("_Z_m_p")]

        self.gen_info_2016_Up = unumpy.uarray(self.df['gen_eff_2016_Up'].values, self.df['err_gen_eff_2016_Up'].values)
        self.gen_info_2017_Up = unumpy.uarray(self.df['gen_eff_2017_Up'].values, self.df['err_gen_eff_2017_Up'].values)
        self.gen_info_2018_Up = unumpy.uarray(self.df['gen_eff_2018_Up'].values, self.df['err_gen_eff_2018_Up'].values)

        self.gen_info_2016_Down = unumpy.uarray(self.df['gen_eff_2016_Down'].values, self.df['err_gen_eff_2016_Down'].values)
        self.gen_info_2017_Down = unumpy.uarray(self.df['gen_eff_2017_Down'].values, self.df['err_gen_eff_2017_Down'].values)
        self.gen_info_2018_Down = unumpy.uarray(self.df['gen_eff_2018_Down'].values, self.df['err_gen_eff_2018_Down'].values)

        self.gen_info_2016 = (self.gen_info_2016_Up + self.gen_info_2016_Down)/2
        self.gen_info_2017 = (self.gen_info_2017_Up + self.gen_info_2017_Down)/2
        self.gen_info_2018 = (self.gen_info_2018_Up + self.gen_info_2018_Down)/2

        print(self.gen_info_2016, self.gen_info_2017, self.gen_info_2018)

        self.front_ID = self.df['Front_ID']
        self.mc_event_ID = self.df['mc_event_ID'].astype(int).astype(str)
        self.file_list = self.front_ID + "_" + self.mc_event_ID

        self.rec_info_2016, self.bootrec_info_2016 = build_rec_eff(self, tt_flag, "2016", dira_cut)
        self.rec_info_2017, self.bootrec_info_2017 = build_rec_eff(self, tt_flag, "2017", dira_cut)
        self.rec_info_2018, self.bootrec_info_2018 = build_rec_eff(self, tt_flag, "2018", dira_cut)

        build_tot_eff(self, self.gen_info_2016, self.rec_info_2016,  self.bootrec_info_2016, "2016")
        build_tot_eff(self, self.gen_info_2017, self.rec_info_2017,  self.bootrec_info_2017, "2017")
        build_tot_eff(self, self.gen_info_2018, self.rec_info_2018,  self.bootrec_info_2018, "2018")

df_c1 = pandas.read_excel(inbook, sheet_name="MC_base")
df_c2 = df_c1.copy()
df_c3 = df_c1.copy()

for dira_cut in [-9 , -4 , -1, 0 , 1, 4, 9]

    total_eff_t = data_analysis(df_c2, "T", dira_cut)
    total_eff_ntat = data_analysis(df_c3, "nTaT", dira_cut)

    outbook = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_eff_{dira_cut}.xlsx"

    with pandas.ExcelWriter(outbook) as writer:
        total_eff_t.df.to_excel(writer, sheet_name='MC_eff_T')
        print(f"Wrote: {outbook}, with sheet: MC_eff_T")
        total_eff_ntat.df.to_excel(writer, sheet_name='MC_eff_nTaT')
        print(f"Wrote: {outbook}, with sheet: MC_eff_nTaT")
