import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
RDF = ROOT.ROOT.RDataFrame

mc_truth_match = "(abs(B_TRUEID) == 511 && abs(KST_TRUEID) == 313))"
# mc_file_list = ["9_zz_11196000.root"]#, "2a_z_11198400", "4a_z_11198401"]

inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_eff.xlsx"
bfbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_mc_base.xlsx"
outbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_effcorr.xlsx"


def get_n_from_ws(self, tt_flag):
    ws_base_file = ROOT.TFile(f"fits/fit_{tt_flag}.root")
    ws = ws_base_file.Get(f"fit_{tt_flag}")
    for n_name in self.df["harris_id"].values:
        n_v = ws.obj(n_name).valv()
        n_e = ws.obj(n_name).err()
        print(n_v)
        print(n_e)
def build_bf_values(self, bf_df):
    bf_value_array = []
    for cbar, c, s in zip(self.df['cbar'], self.df['c'], self.df['s']):

        n_cbar = bf_df.loc[bf_df['decay'] == cbar, 'bf_value'].values[0]
        u_cbar = bf_df.loc[bf_df['decay'] == cbar, 'err_bf_value'].values[0]

        print (f"{cbar} bf is: {n_cbar}")

        n_c = bf_df.loc[bf_df['cc_decay'] == c, 'bf_value'].values[0]
        u_c = bf_df.loc[bf_df['cc_decay'] == c, 'err_bf_value'].values[0]

        print (f"{c} bf is: {n_c}")

        cbar_value = ufloat(n_cbar, u_cbar)
        c_value = ufloat(n_c, u_c)

        e1_n_c = 1
        e1_u_c = 0
        e2_n_c = 1
        e2_u_c = 0

        if cbar == "D*-->D0bar_pi-" or cbar =="D*0bar":
            e1_n_c = bf_df.loc[bf_df['decay'] == 'D0bar', 'bf_value'].values[0]
            e1_u_c = bf_df.loc[bf_df['decay'] == 'D0bar', 'err_bf_value'].values[0]
            print ("got e1")
        if cbar == "D*-->D-pi0":
            e1_n_c = bf_df.loc[bf_df['decay'] == 'D-', 'bf_value'].values[0]
            e1_u_c = bf_df.loc[bf_df['decay'] == 'D-', 'err_bf_value'].values[0]
            print ("got e1")
        if c == "D*+->D0_pi+" or c == "D*0":
            e2_n_c = bf_df.loc[bf_df['cc_decay'] == 'D0', 'bf_value'].values[0]
            e2_u_c = bf_df.loc[bf_df['cc_decay'] == 'D0', 'err_bf_value'].values[0]
            print ("got e2")
        if c == "D*+->D+_pi0":
            e2_n_c = bf_df.loc[bf_df['cc_decay'] == 'D+', 'bf_value'].values[0]
            e2_u_c = bf_df.loc[bf_df['cc_decay'] == 'D+', 'err_bf_value'].values[0]
            print ("got e2")

        e1 = ufloat(e1_n_c, e1_u_c)
        e2 = ufloat(e2_n_c, e2_u_c)

        print(f"e1 is: {e1}")
        print(f"e2 is: {e2}")

        if s == "K*0":
            k_value = ufloat(2/3, 0)
            b_value = ufloat(1,0)
        if s == "K":
            k_value = ufloat(1, 0)
            if cbar == "D-":
                b_value = ufloat(1.07E-03, 1.10E-04)
            if cbar == "D0bar":
                b_value = ufloat(1.31E-03, 7.00E-05)

        value = c_value*cbar_value*e1*e2*k_value*b_value
        print (f"{cbar} + {c} is: {value}")
        bf_value_array.append(value)
    print(bf_value_array)
    self.df['bf_values'] =  unumpy.nominal_values(bf_value_array)
    self.df['err_bf_values'] = unumpy.std_devs(bf_value_array)
    return bf_value_array
def build_ecy(self):
    ecy = self.data_yield/(self.tea*self.bfa)
    self.df['eff_corrected_yield'] = unumpy.nominal_values(ecy)
    self.df['err_eff_corrected_yield'] = unumpy.std_devs(ecy)
    return ecy
def sanity_check(self):

    df_norm = self.df.loc[self.df["harris_id"].str.contains("norm")]
    norm7_row = df_norm.loc[df_norm['harris_id'] == 'n_norm7'].squeeze()
    norm8_row = df_norm.loc[df_norm['harris_id'] == 'n_norm8'].squeeze()
    #
    print(norm7_row)

    eff_b0 = ufloat(norm8_row['total_eff'], norm8_row['err_total_eff'])
    eff_bp = ufloat(norm7_row['total_eff'], norm7_row['err_total_eff'])

    nb0 = ufloat(norm8_row['yield'], norm8_row['err_yield'])
    nbp = ufloat(norm7_row['yield'], norm7_row['err_yield'])

    bfb0 = ufloat(1.07E-03, 1.10E-04)
    bfbp = ufloat(1.31E-03, 7.00E-05)
    bfd0 = ufloat(0.0395, 0.00031)
    bfdp = ufloat(0.0938, 0.0016)

    t1 = (eff_b0/eff_bp)*(nbp/nb0)
    t2 = (bfbp*bfd0)/(bfb0*bfdp)

    print(t1)
    print(t2)
class data_analysis():
    def __init__(self, data, bf_df, tt_flag):
        self.df = pandas.DataFrame(data)
        self.data_yield = get_n_from_ws(self, tt_flag)
        print(self.df["harris_id"].values)
        # self.data_yield = unumpy.uarray(self.df['yield'].values, self.df['err_yield'].values)
        # self.tea = unumpy.uarray(self.df['total_eff'].values, self.df['err_total_eff'].values)
        # print(self.data_yield)
        # print(self.tea)
        # self.bfa = build_bf_values(self, bf_df)
        # self.ecy = build_ecy(self)
        # self.sanity = sanity_check(self)

df_tot = pandas.read_excel(inbook, sheet_name="MC_eff_ToT")
df_t = pandas.read_excel(inbook, sheet_name="MC_eff_T")
df_ntat = pandas.read_excel(inbook, sheet_name="MC_eff_nTaT")

# df_y =  pandas.read_csv(normyieldtxt, sep=" ", header=None)
# df_y.columns = ["harris_id", "tt_flag", "yield", "err_yield"]
#
# df_y_tot = df_y.loc[df_y["tt_flag"] == "tot"]
# df_y_t = df_y.loc[df_y["tt_flag"] == "t"]
# df_y_ntat = df_y.loc[df_y["tt_flag"] == "ntat"]

bf_df = pandas.read_excel(bfbook, sheet_name="BF_Values_PDG")

# df_tot_all = pandas.merge(df_tot, df_y_tot, on=["harris_id"])
# df_t_all = pandas.merge(df_t, df_y_t, on=["harris_id"])
# df_ntat_all = pandas.merge(df_ntat, df_y_ntat, on=["harris_id"])

total_eff_tot = data_analysis(df_tot, bf_df, "ToT")
total_eff_t = data_analysis(df_t, bf_df, "T")
total_eff_ntat = data_analysis(df_ntat, bf_df, "nTaT")
#
# with pandas.ExcelWriter(outbook) as writer:
#     df_tot_all.to_excel(writer, sheet_name='Data_ToT')
#     print(f"Wrote: {outbook}, with sheet: Data_ToT")
#     df_t_all.to_excel(writer, sheet_name='Data_T')
#     print(f"Wrote: {outbook}, with sheet: Data_T")
#     df_ntat_all.to_excel(writer, sheet_name='Data_nTaT')
#     print(f"Wrote: {outbook}, with sheet: Data_nTaT")
