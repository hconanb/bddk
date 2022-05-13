import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *

def eff_table(flag):

    for flag in ["trigger"]:

        file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_efficiencies/build_root_effs/*{flag}_eff_txt")

        if flag == "gen":
            mcdf_cols = ['Scheme ID','Year','Number Accepted','$\epsilon_{geometrical}$']
            caption_line = "MC Efficiencies: Generator"
            the_column_format = 'cccc'
            label_line_1 = "tab:genstrip_1"
            label_line_2 = "tab:genstrip_2"

        if flag == "can":
            mcdf_cols = ['Scheme ID','Year','$\epsilon_{stripping}$', '$\epsilon_{offline}$','$\epsilon_{D Window}$']
            caption_line = "MC Candidate Efficiencies"
            the_column_format = 'ccccc'
            label_line_1 = "tab:can_offdwin_1"
            label_line_2 = "tab:can_offdwin_2"

        if flag == "event":
            mcdf_cols = ['Scheme ID','Year','$\epsilon_{stripping}$', '$\epsilon_{offline}$','$\epsilon_{D Window}$']
            caption_line = "MC Event Efficiencies"
            the_column_format = 'ccccc'
            label_line_1 = "tab:event_offdwin_1"
            label_line_2 = "tab:event_offdwin_2"

        if flag == "trigger":
            file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_efficiencies/boot_effs/*")
            mcdf_cols = ['Scheme ID','Year','Final Bootstraped TOS Events','Final Bootstraped TIS Events']
            caption_line = "MC Efficiencies: Disjoint Trigger Regions"
            the_column_format = 'cccc'
            label_line_1 = "tab:tr_1"
            label_line_2 = "tab:tr_2"

        df = pandas.DataFrame()
        for file in file_list:
            frame = pandas.read_csv(file, usecols = mcdf_cols)
            df = df.append(frame)

        df_i = df.set_index(['Scheme ID','Year'])

        scheme = [1,'2a,3a','4a','4b','4c','4d',5,6,'7a','7b','7c','8a','8b','8c',9,10,12,"norm7","norm8"]

        df_i = df_i.reindex(scheme, level=0)

        df_i_1 = df_i.iloc[:39,:]
        df_i_2 = df_i.iloc[39:69,:]
        # df_i_3 = df_i.iloc[42:63,:]

        latex_1 = df_i_1.to_latex(index = True, multirow = True, escape=False, column_format=the_column_format, label = label_line_1, caption = caption_line, )
        latex_2 = df_i_2.to_latex(index = True, multirow = True, escape=False, column_format=the_column_format, label = label_line_2, caption = caption_line, )
        # latex_3 = df_i_3.to_latex(index = True, multirow = True, escape=False, column_format=the_column_format, label = label_line, caption = caption_line, )

        with open(f'tex_files/t_{flag}.tex','w') as tf:
            tf.write(latex_1)
            tf.write(latex_2)
            # tf.write(latex_3)
        # print(latex)
def build_d_window_dict():
    for d1_flag in ["mp", "z", "d0k3pi", "dst"]:
        if d1_flag == "dst":
            dst_flag = True
        else:
            dst_flag = False

        d1_mstart, d1_std  = get_dwindow_values("temp", d1_flag, "z", dst_flag, rflag = "table")

        d1minw = d1_mstart - 2*d1_std
        d1maxw = d1_mstart + 2*d1_std
        d1Lsbmin = d1_mstart - 5*d1_std
        d1Lsbmax = d1_mstart - 3*d1_std
        d1Rsbmin = d1_mstart + 3*d1_std
        d1Rsbmax = d1_mstart + 5*d1_std

        if d1_flag == "mp":
            d1_name = "\Dm"
        if d1_flag == "z":
            d1_name = "\Dzb"
        if d1_flag == "dst" :
            d1_name = "\Dstarm - \Dzb"
        if d1_flag == "d0k3pi":
            d1_name = r"\Dz\rightarrow\kaon\pion\pion\pion"

        d1_dict_temp = {
                        'D Candidate': f"${d1_name}$",
                        'Fit Mean [MeV]': f"${d1_mstart:.2f}$",
                        'Signal Window [MeV, MeV]': f'$[{d1minw:.2f}, {d1maxw:.2f}]$',
                        'Sideband Window [MeV, MeV]': f'$[{d1Lsbmin:.2f}, {d1Lsbmax:.2f}]$ and $[{d1Rsbmin:.2f}, {d1Rsbmax:.2f}]$'
                        }

        d1_df = pandas.DataFrame(d1_dict_temp, index=[0])
        D1_txt = f"d_window_txt_files/{d1_flag}.txt"
        d1_df.to_csv(D1_txt)
def build_d_window_latex():
    file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/build_latex_tables/d_window_txt_files/*")
    df = pandas.DataFrame()

    for file in file_list:
        df_cols = ['D Candidate','Fit Mean [MeV]','Signal Window [MeV, MeV]', 'Sideband Window [MeV, MeV]']
        frame = pandas.read_csv(file, usecols = df_cols)
        df = df.append(frame)
        # index_col.append(iname)
    # df.index = index_col
    # print(df)
    # df_i = pandas.MultiIndex.from_frame(df)
    caption_line = "Summary of Mass Window Cuts on D Mesons"
    # df_i = df_i.rename(columns={"Candidates in Z_z_z": "Candidates in ZZ", 'Candidates in P_z_pst': 'Candidates in ST', 'Number of Duplicate Candidates': 'Number of Shared Candidates'})
    df_i = df.set_index(['D Candidate'])
    # df_i['Spectrum ID'] = df_i['Spectrum ID'].astype("category")
    # df_i['Spectrum ID'].cat.set_categories(sorter)
    # df_i = df_i.sort_values('Spectrum ID')
    # df_i = df_i.reindex(sorter,level=0)
    label_line = "tab:dwindows"

    latex = df_i.to_latex(index = True, multirow = True, escape=False, column_format='ccccc', caption = caption_line, label = label_line)
    with open(f'tex_files/dwin.tex','w') as tf:
        tf.write(latex)
def build_lap_latex():
    file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/build_rootfiles/lap_txt_files/*")
    df = pandas.DataFrame()

    for file in file_list:
        df_cols = ['Spec','Year','Candidates in Z_z_z','Candidates in P_z_pst','Shared Candidates Removed from Z_z_z']
        frame = pandas.read_csv(file, usecols = df_cols)
        df = df.append(frame)
        # index_col.append(iname)
    # df.index = index_col
    # print(df)
    # df_i = pandas.MultiIndex.from_frame(df)
    caption_line = "Summary of Canidates removed from ZZ spectrum"
    label_line = "tab:lap"
    df = df.rename(columns={"Candidates in Z_z_z": "Candidates in ZZ", 'Candidates in P_z_pst': 'Candidates in ST', 'Shared Candidates Removed from Z_z_z': 'Number of Shared Candidates Removed from ZZ'})

    df_i = df.set_index(['Spec','Year'])
    df_i = df_i.rename(index={"Z_z_z": "Data", '04_Z_z_z_11198023': '4d', '07_Z_z_z_12197045': '7b', '08_Z_z_z_12197423': '8b'})

    latex = df_i.to_latex(index = True, multirow = True, escape=False, column_format='ccccc', caption = caption_line, label = label_line)
    with open(f'tex_files/lap.tex','w') as tf:
        tf.write(latex)
def build_can_latex():
    file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mult_can/mult_can_txt_files/*")
    df = pandas.DataFrame()

    for file in file_list:
        df_cols = ["Spectrum","Year","Canidates/Unique Event","Canidates Cut"]
        frame = pandas.read_csv(file, usecols = df_cols)
        df = df.append(frame)
        # index_col.append(iname)
    # df.index = index_col
    # print(df)
    # df_i = pandas.MultiIndex.from_frame(df)
    caption_line = "Multiple Canidates"
    label_line = "tab:can"

    df_i = df.set_index(["Spectrum","Year"])
    df_i = df_i.rename(index={'P_z_pst': 'ST', 'Z_z_z': 'ZZ', 'Z_m_p': 'Z', 'P_z_p': 'P', 'M_m_z': 'M', 'Zs_sm_p': 'Zs'})

    latex = df_i.to_latex(index = True, multirow = True, escape=False, column_format='ccccc', caption = caption_line, label = label_line)
    print(latex)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--mc_efficiency', action='store_true',
                    help="Gen MC_EFF Table")
    parser.add_argument('--eff_flag',default='gen', choices=['gen', 'can', 'event'])
    parser.add_argument('--d_window_dict', action='store_true',
                    help='Gen D_WINDOW Table')
    parser.add_argument('--d_window_latex', action='store_true',
                    help='Gen D_WINDOW latex')
    parser.add_argument('--lap_latex', action='store_true')
    parser.add_argument('--can_latex', action='store_true')
    # parser.add_argument('--charmless', action='store_true',
    #                 help='Print charmless background estimation')
    # parser.add_argument('--sigma', action='store_true',
    #                 help='Print sigma scale parameters')
    # parser.add_argument('--yields', action='store_true',
    #                 help='Print signal yields')
    # parser.add_argument('--params', action='store_true',
    #                 help='Print fit parameters')
    # parser.add_argument('--ratios', action='store_true',
    #                 help='Print fit parameters')
    # parser.add_argument('--expected', action='store_true',
    #                 help='Print expected yields')
    args = parser.parse_args()
    efficiency = args.mc_efficiency
    eff_flag = args.eff_flag
    d_window_dict = args.d_window_dict
    d_window_latex = args.d_window_latex
    lap_latex = args.lap_latex
    can_latex = args.can_latex
    if efficiency:
        eff_table(eff_flag)

    if d_window_dict:
        build_d_window_dict()

    if d_window_latex:
        build_d_window_latex()

    if lap_latex:
        build_lap_latex()

    if can_latex:
        build_can_latex()
