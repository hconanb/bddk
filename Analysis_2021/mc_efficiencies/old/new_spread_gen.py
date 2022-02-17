import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

mcdf_cols = ['Scheme ID','Year','epsilon_geometrical','epsilon_final_ToT','epsilon_final_T','epsilon_final_nTaT']

def eff_spread():
    file_list = glob.glob(f"../build_rootfiles/eff_txt_files/*_test")
    df = pandas.DataFrame()
    index_col = []
    for file in file_list:
        frame = pandas.read_csv(file, usecols = mcdf_cols)
        df = df.append(frame)
    df_i= df.set_index(['Scheme ID','Year'])
    # print(df_i)
    df_i.sort_index(inplace=True)
    # df_i = df_i[[
    # "$\epsilon_{geometrical}$",
    # "$\epsilon_{final_ToT}$",
    # "$\epsilon_{final_T}$",
    # "$\epsilon_{final_nTaT}$"]]
    # # print(df_i)
    # df_i.columns = df_i.columns.str.replace(' ', '')
    # df_i.columns = df_i.columns.str.replace('$', '')
    # df_i.columns = df_i.columns.str.replace('\\', '')
    # df_i.columns = df_i.columns.str.replace('{', '')
    # df_i.columns = df_i.columns.str.replace('}', '')
    # # df_i = df_i.replace('$', '', regex = True)
    # df_i = df_i.replace("\$","", regex=True)
    # test = df_i["epsilon_geometrical"].str.split(" ", n = 2, expand=True)
    print(df_i)
    df_i.to_excel("text.xlsx")

eff_spread()
