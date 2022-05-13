import pandas
import glob as glob
import ROOT
import numpy as np
from uncertainties import ufloat
from uncertainties import covariance_matrix, correlation_matrix
from uncertainties import unumpy
from uncertainties import ufloat_fromstr
import datetime
import os
ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")
ROOT.gStyle.SetPalette(ROOT.kBird)

RDF = ROOT.ROOT.RDataFrame

dm = ufloat(0.0938, 0.0016)
dsm = ufloat(0.0539, 0.0015)

d0bar = ufloat(0.0395, 0.00031)
d0bar4 = ufloat(0.0823, 0.0014)

dstmdmpi0 = ufloat(0.323, 0.005)
dstmd0pim = ufloat(0.677, 0.005)

B_B0_norm7 = ufloat(1.31e-03, 7.00e-05)
B_B0_norm8 = ufloat(1.07e-03, 1.10e-04)

# zs_data_fit_file = ROOT.TFile("../data_fit/data_files/corr_test_2_4.root")
# data_fit_file = ROOT.TFile("../data_fit/data_files/corr_test_2_12_fa.root")
# dws = data_fit_file.Get("test_2_12_fa")
# fitresult = dws.obj("fitresult_super_fit_Pdf_all_data_sets")

# dws_nds = data_fit_file.Get("test_2_4_nds")
# fitresult_nds = dws_nds.obj("fitresult_super_fit_Pdf_all_data_sets")

# def get_norm_param(norm_id, trigger_flag, year, df):
#
#     norm_id = norm_id.split("_")[0]
#     norm_base_file = ROOT.TFile(f"../normalization_fit/fit_files/{norm_id}_{year}_{trigger_flag}.root")
#     norm_ws = norm_base_file.Get(f"{norm_id}")
#     n_base = norm_ws.var(f"n_{norm_id}_signal")
#     fitresult = norm_ws.obj(f"fitresult_{norm_id}_fit_{norm_id}_data")
#     # errerr = n_base.getPropagatedError(fitresult)
#     # print(errerr, n_base.getError())
#     n_val = ufloat(n_base.getValV(), n_base.getError())
#     df.loc[(norm_id, int(year)),f"Data Yield {trigger_flag}"] = n_val
#     # return [n_val.n, n_val.s, trigger_flag, year]
#
# def ufloat_tuple(y, num, scheme_id):
#     bf_nom = num.n
#     bf_err = num.s
#     nom_obj = dws.obj(y)
#     y_nom = nom_obj.getValV()
#     y_err = nom_obj.getPropagatedError(fitresult)
#     return [y_nom, y_err, bf_nom, bf_err, scheme_id]

# def ufloat_ds_tuple(y, num, scheme_id):
#     bf_nom = num.n
#     bf_err = num.s
#     nom_obj = dws_ds.obj(y)
#     y_nom = nom_obj.getValV()
#     y_err = nom_obj.getPropagatedError(fitresult_ds)
#     return [y_nom, y_err, bf_nom, bf_err, scheme_id]

# a_01_z = ufloat_tuple("Z_m_p_01_yield",dm*dm, 1)
# a_02_z = ufloat_tuple("Z_m_p_02_yield",dstmdmpi0 * dm * dm, "2a,3a")
# a_03_z = ufloat_tuple("Z_m_p_03_yield",dstmdmpi0 * dm * dm, "2a,3a")
# a_04_z = ufloat_tuple("Z_m_p_04_yield",dstmdmpi0 * dstmdmpi0 * dm * dm, "4a")
#
# a_05_p = ufloat_tuple("P_z_p_05_yield",dm * d0bar, 5)
# a_06_p = ufloat_tuple("P_z_p_06_yield",dm * d0bar, 6)
# a_07_p = ufloat_tuple("P_z_p_07_yield",d0bar * dstmdmpi0 * dm, "7a")
# a_08_p = ufloat_tuple("P_z_p_08_yield",d0bar * dstmdmpi0 * dm, "8a")
#
# a_09_zz = ufloat_tuple("Z_z_z_09_yield",d0bar * d0bar, 9)
# a_10_zz = ufloat_tuple("Z_z_z_10_yield",d0bar * d0bar, 10)
# a_12_zz = ufloat_tuple("Z_z_z_12_yield",d0bar * d0bar, 12)
#
# a_13_s = ufloat_tuple("Zs_sm_p_13_yield",dsm * dm, 13)
# a_14_s = ufloat_tuple("Zs_sm_p_14_yield",dsm * dm, 14)
# a_15_s = ufloat_tuple("Zs_sm_p_15_yield",dstmdmpi0 * dsm * dm, 15)
# a_16_s = ufloat_tuple("Zs_sm_p_16_yield",dstmdmpi0 * dsm * dm, 16)
#
# bf_list = [a_01_z,
# a_02_z,
# a_03_z,
# a_04_z,
# a_05_p,
# a_06_p,
# a_07_p,
# a_08_p,
# a_09_zz,
# a_10_zz,
# a_12_zz,
# a_13_s,
# a_14_s,
# a_15_s,
# a_16_s,]
#
# # print(bf_list + bf_s)
#
# df_raw = pandas.DataFrame(bf_list, columns = ['Fit Yield', 'Fit Error', 'BF_factor', 'BF_error', 'Scheme ID'])
# df_raw = df_raw.set_index(['Scheme ID'])

# print(df_raw)

gen_file_list = glob.glob(f"build_root_effs/*gen_eff_txt")
gen_mcdf_cols = ['Scheme ID','Year','Number Accepted','$\epsilon_{geometrical}$']

def get_err(df_slice, column):
    a = df_slice[column].to_numpy()
    b = (a[0] + a[1] + a[2])/3
    # c = b.mean()
    return b

boot_file_list = glob.glob(f"boot_effs/*")
boot_mcdf_cols = ['Scheme ID','Year','Decay Description','Final Bootstraped TOS Events','Final Bootstraped TIS Events']

df_gen = pandas.DataFrame()
df_boot = pandas.DataFrame()

for gen_file in gen_file_list:
    frame_1 = pandas.read_csv(gen_file, usecols = gen_mcdf_cols)
    df_gen = df_gen.append(frame_1)
for boot_file in boot_file_list:
    frame_3 = pandas.read_csv(boot_file, usecols = boot_mcdf_cols)
    df_boot = df_boot.append(frame_3)

df_gen = df_gen.rename(columns={'$\epsilon_{geometrical}$': 'Generator (%)'})
df_gen = df_gen.set_index(['Scheme ID','Year'])
df_boot = df_boot.set_index(['Scheme ID','Year'])
df = pandas.concat([df_gen, df_boot], axis=1, join="inner")

# n7bf = d0bar * d0bar4
# n8bf = dm*d0bar4

# df["Data Yield T"] = ""
# df["Data Yield nTaT"] = ""
# df["Branching Fraction D1_X_D2"] = ""
# df.loc["norm7","Branching Fraction D1_X_D2"] = n7bf
# df.loc["norm8","Branching Fraction D1_X_D2"] = n8bf

# for n in ["norm7","norm8"]:
#     for trigger_flag in ["T", "nTaT"]:
#         for year in ["2016", "2017", "2018"]:
#             get_norm_param(n, trigger_flag, year, df)
#
# df_final = pandas.DataFrame()
#
# df_final["MC_TOS_Yield"] = df['Final Bootstraped TOS Events'].apply(ufloat_fromstr)
# df_final["MC_TIS_Yield"] = df['Final Bootstraped TIS Events'].apply(ufloat_fromstr)
# df_final["Gen Events"] = unumpy.uarray(df['Number Accepted'], df['Number Accepted'].apply(lambda x : np.sqrt(x)/x))
# #
# df_final["Gen Eff"] = (df['Generator (%)'].apply(ufloat_fromstr))/100
# #
# df_final["TOS Eff"] = (df_final["MC_TOS_Yield"]/df_final["Gen Events"])*df_final["Gen Eff"]
# df_final["TIS Eff"] = (df_final["MC_TIS_Yield"]/df_final["Gen Events"])*df_final["Gen Eff"]
#
# df_final["TOS Eff n"] = df_final["TOS Eff"].apply(lambda x : x.n)
# df_final["TOS Eff err"] = df_final["TOS Eff"].apply(lambda x : x.s)
# df_final["TIS Eff n"] = df_final["TIS Eff"].apply(lambda x : x.n)
# df_final["TIS Eff err"] = df_final["TIS Eff"].apply(lambda x : x.s)
#
# df_final["TOS Data n"] = df["Data Yield T"].apply(lambda x : x.n)
# df_final["TOS Data err"] = df["Data Yield T"].apply(lambda x : x.s)
# df_final["TIS Data n"] = df["Data Yield nTaT"].apply(lambda x : x.n)
# df_final["TIS Data err"] = df["Data Yield nTaT"].apply(lambda x : x.s)
# df_final["BF D1 x D2 n"] = df["Branching Fraction D1_X_D2"].apply(lambda x : x.n)
# df_final["BF D1 x D2 err"] = df["Branching Fraction D1_X_D2"].apply(lambda x : x.s)

# df_final["Eff"] = ((df_final["MC_TOS_Yield"] + df_final["MC_TIS_Yield"])/df_final["Gen Events"])*df_final["Gen Eff"]
#
# # df_final["TOS Eff"] = (df_final["MC_TOS_Yield"]/df_final["Gen Events"])*df_final["Gen Eff"]
# # df_final["TIS Eff"] = (df_final["MC_TIS_Yield"]/df_final["Gen Events"])*df_final["Gen Eff"]
# # df_final["TOS Eff n"] = df_final["TOS Eff"].apply(lambda x : x.n)
# # df_final["TOS Eff err"] = df_final["TOS Eff"].apply(lambda x : x.s)
# # df_final["TIS Eff n"] = df_final["TIS Eff"].apply(lambda x : x.n)
# # df_final["TIS Eff err"] = df_final["TIS Eff"].apply(lambda x : x.s)
#
# # df_group["Gen Eff"] = df_final.groupby(level = ['Scheme ID']).apply(get_err, "Gen Eff")
# # df_group["TOS"] = df_final.groupby(level = ['Scheme ID']).apply(get_err, "MC_TOS_Yield")
# # df_group["TIS"] = df_final.groupby(level = ['Scheme ID']).apply(get_err, "MC_TIS_Yield")
# # df_group["Generator Events"] = df_final.groupby(level = ['Scheme ID']).apply(get_err, "Gen Events")
#
# # df_group["Rec Eff"] = (df_group["TOS"] + df_group["TIS"]) / df_group["Generator Events"]
# df_group["Final Eff"] = df_final.groupby(level = ['Scheme ID']).apply(get_err, "Eff")
# df_group["Final Eff Nom"] = df_group["Final Eff"].apply(lambda x : x.n)
# df_group["Final Eff Err"] = df_group["Final Eff"].apply(lambda x : x.s)
#
# norm7_y = df_final.loc["norm7", ["TOS Eff n","TOS Eff err", "TIS Eff n", "TIS Eff err"]]
# norm8_y = df_final.loc["norm8", ["TOS Eff n","TOS Eff err", "TIS Eff n", "TIS Eff err"]]

# with pandas.ExcelWriter("Arvind_Norm_02152022_numbers.xlsx") as writer:
#     df_final.to_excel(writer, sheet_name='Norm Breakdown')
#     # df_group[["Final Eff Nom", "Final Eff Err"]].to_excel(writer, sheet_name='MC Effs')
#     # df_norm7.to_excel(writer, sheet_name='Norm7_fit')
#     # df_norm8.to_excel(writer, sheet_name='Norm8_fit')
#     # norm7_y.to_excel(writer, sheet_name='Norm7_eff')
#     # norm8_y.to_excel(writer, sheet_name='Norm8_eff')
#
with pandas.ExcelWriter("Harris_numbers.xlsx") as writer:
    df.to_excel(writer)
    # df_group[["Final Eff Nom", "Final Eff Err"]].to_excel(writer, sheet_name='MC Effs')
    # df_norm7.to_excel(writer, sheet_name='Norm7_fit')
    # df_norm8.to_excel(writer, sheet_name='Norm8_fit')
    # norm7_y.to_excel(writer, sheet_name='Norm7_eff')
    # norm8_y.to_excel(writer, sheet_name='Norm8_eff')
