import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
RDF = ROOT.ROOT.RDataFrame

mc_truth_match = "(abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313))"

inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_mc_base.xlsx"
outbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_eff.xlsx"

dwindow_file = ROOT.TFile(f"dmass/d_mass_fits.root","READ")
dwindow_ws = dwindow_file.Get("d_mass_fits")

def build_rec_eff(self, file_list, tt_flag):

    eff_array = []
    for file_name in file_list:
        print(file_name)
        mcfile = ROOT.TFile(mc_basepath + file_name)

        if "norm" not in file_name:
            flag = file_name.split("_")[1]
            id = file_name.split("_")[2]
            print(id)
        else:
            flag  = file_name.split("_")[0]

        print("getting reconstruction eff for :", file_name)

        rec_tree = mcfile.Get("DecayTreeTuple")
        gen_tree = mcfile.Get("MCDecayTreeTuple")

        rdf_rec = RDF(rec_tree)

        m_cut = get_dwindow_values(dwindow_ws, flag)
        print(flag, m_cut)

        rdf_rec_current_cut = rdf_rec.Filter(f"{m_cut} && {trigger_cuts[tt_flag]}")
        rdf_gen = RDF(gen_tree)

        rec_count = rdf_rec_current_cut.Count().GetValue()
        gen_count = rdf_gen.Count().GetValue()

        rec_all = ufloat(rec_count, math.sqrt(rec_count))
        gen_all = ufloat(gen_count, math.sqrt(gen_count))

        eff = rec_all/gen_all
        print(eff)
        eff_array.append(eff)

    self.df['rec_eff'] =  unumpy.nominal_values(eff_array)
    self.df['err_rec_eff'] = unumpy.std_devs(eff_array)
    return eff_array
def build_tot_eff(self):
    tea = self.gen_info*self.filter_info*self.rec_info
    self.df['total_eff'] =  unumpy.nominal_values(tea)
    self.df['err_total_eff'] = unumpy.std_devs(tea)
    return tea
class data_analysis():
    def __init__(self, data, tt_flag):
        self.df = pandas.DataFrame(data)
        self.gen_info = unumpy.uarray(self.df['gen_eff'].values, self.df['err_gen_eff'].values)
        self.filter_info = unumpy.uarray(self.df['filter_eff'].values, self.df['err_filter_eff'].values)
        self.file_names = self.df['root_file_name']
        self.rec_info = build_rec_eff(self, self.file_names, tt_flag)
        self.tea = build_tot_eff(self)

df_c1 = pandas.read_excel(inbook, sheet_name="MC_base")
df_c2 = df_c1.copy()
df_c3 = df_c1.copy()

total_eff_tot = data_analysis(df_c1, "ToT")
# total_eff_t = data_analysis(df_c2, "T")
# total_eff_ntat = data_analysis(df_c3, "nTaT")

with pandas.ExcelWriter(outbook) as writer:
    df_c1.to_excel(writer, sheet_name='MC_eff_ToT')
    print(f"Wrote: {outbook}, with sheet: MC_eff_TOSorTIS")
    df_c2.to_excel(writer, sheet_name='MC_eff_T')
    print(f"Wrote: {outbook}, with sheet: MC_eff_T")
    df_c3.to_excel(writer, sheet_name='MC_eff_nTaT')
    print(f"Wrote: {outbook}, with sheet: MC_eff_nTaT")


    #             rec_count += 1
    #         else:
    #             bad_count +=1
    #     #print(str(rec_count) + " vs " + str(bad_count))
    #     old_eff = rec_count/gen_events
    #     print("returning old")
    #     return old_eff
    #
    # for ev in rec_tree:
    #     if e_value != ev.RD_org_eventNumber or r_value != ev.RD_org_runNumber:
    #         e_value = ev.RD_org_eventNumber
    #         r_value = ev.RD_org_runNumber
    #         bootrec_count += count*np.random.poisson(1)
    #             #print (str(e_value) + " " + str(r_value) +"\n")
    #             #print (str(count))
    #         count = 0
    #     if abs(ev.D1_M - c1mass) < dwindow and abs(ev.D2_M - c2mass) < dwindow:
    #         count += 1
    #     #else:
    #         #print("BAD Hombre\n")

    # bootrec_count += count*np.random.poisson(1)
