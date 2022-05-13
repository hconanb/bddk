import sys
import os

basedir = os.getcwd().split('latex_tables')[0]
sys.path.append(basedir)

from essential_functions import *

def get_ws_var(ws, name, fit_result):
    var = ws.var(f"{name}")
    if not var:
        return "NA"
    else:
        var_err = var.getPropagatedError(fit_result)
        return (ufloat(var.getValV(), var_err))
def build_d_window():
    for d1_flag in ["mp", "z", "d0k3pi", "dst"]:

        d1_mstart, d1_std  = get_dwindow_values(d1_flag, rflag = "table")

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

        d1_df = pd.DataFrame(d1_dict_temp, index=[0])
        D1_txt = f"d_window_txt/{d1_flag}.txt"
        d1_df.to_csv(D1_txt)

    file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/latex_tables/d_window_txt/*")
    df = pd.DataFrame()

    for file in file_list:
        df_cols = ['D Candidate','Fit Mean [MeV]','Signal Window [MeV, MeV]', 'Sideband Window [MeV, MeV]']
        frame = pd.read_csv(file, usecols = df_cols)
        df = df.append(frame)

    caption_line = "Summary of Mass Window Cuts on D Mesons"
    df_i = df.set_index(['D Candidate'])
    label_line = "tab:dwindows"

    latex = df_i.to_latex(index = True, multirow = True, escape=False, column_format='ccccc', caption = caption_line, label = label_line)
    with open(f'tex_files/dwin.tex','w') as tf:
        tf.write(latex)
def build_mc_shape_table():
    for spec, shape_flag in ids_to_bestfit_dict.items():
        fit_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/mc_files/{id_to_spec_dict[spec]}.root")
        ws = fit_file.Get("fit_ws")
        fit_result = ws.obj(f"fitresult_{spec}_fit_{shape_flag}_{spec}_events")

        # fit_pdf = ws.pdf(f"{spec}_fit_{shape_flag}")
        mean = get_ws_var(ws, f"mean_{spec}_{shape_flag}", fit_result)

        if shape_flag == "DG":
            width_A = get_ws_var(ws, f"width_a_{spec}_{shape_flag}", fit_result)
            width_B = get_ws_var(ws, f"width_b_{spec}_{shape_flag}", fit_result)
            alpha_1 = get_ws_var(ws, f"alpha_1_{spec}_{shape_flag}", fit_result)
            alpha_2 = get_ws_var(ws, f"alpha_2_{spec}_{shape_flag}", fit_result)
            frac = get_ws_var(ws, f"a_frac_{spec}_{shape_flag}", fit_result)
        if shape_flag == "BGEP":
            width_A = get_ws_var(ws, f"width_L_{spec}_{shape_flag}", fit_result)
            width_B = get_ws_var(ws, f"width_R_{spec}_{shape_flag}", fit_result)
            alpha_1 = get_ws_var(ws, f"alpha_1_{spec}_{shape_flag}", fit_result)
            alpha_2 = get_ws_var(ws, f"alpha_2_{spec}_{shape_flag}", fit_result)
            frac = get_ws_var(ws, f"a_frac_{spec}_{shape_flag}", fit_result)
        if "_fr" in shape_flag:
            g_width_A = get_ws_var(ws, f"width_{spec}_{shape_flag}_a", fit_result)
            g_alpha_A = get_ws_var(ws, f"alpha_{spec}_{shape_flag}_a", fit_result)
            g_width_B = get_ws_var(ws, f"width_{spec}_{shape_flag}_b", fit_result)
            g_alpha_B = get_ws_var(ws, f"alpha_{spec}_{shape_flag}_b", fit_result)
            frac = get_ws_var(ws, f"a_frac_{spec}_{shape_flag}", fit_result)
            if shape_flag == "GAddBGEP_fr" or shape_flag == "GEPAddBGEP_fr":
                bg_width_1 = get_ws_var(ws, f"width_L_{spec}_{shape_flag}_b", fit_result)
                bg_width_2 = get_ws_var(ws, f"width_R_{spec}_{shape_flag}_b", fit_result)
                bg_alpha_1 = get_ws_var(ws, f"alpha_1_{spec}_{shape_flag}_b", fit_result)
                bg_alpha_2 = get_ws_var(ws, f"alpha_2_{spec}_{shape_flag}_b", fit_result)
            else:
                bg_width_1 = get_ws_var(ws, f"width_1_{spec}_{shape_flag}_b", fit_result)
                bg_width_2 = get_ws_var(ws, f"width_2_{spec}_{shape_flag}_b", fit_result)
                bg_alpha_1 = get_ws_var(ws, f"alpha_1_{spec}_{shape_flag}_b", fit_result)
                bg_alpha_2 = get_ws_var(ws, f"alpha_2_{spec}_{shape_flag}_b", fit_result)

            stratintab = id_to_scheme_dict[id_to_spec_dict[spec]].split("_fr")[0]
            shapeintab1 = shape_flag.split("Add")[0]
            shapeintab2 = shape_flag.split("_fr")[0].split("Add")[1]
            dict_temp_1 = {
                            'MC Scheme': f"{stratintab}",
                            'Shape': f"{shapeintab1}",
                            'Mean': f'{mean}',
                            'Fraction of Fit': f'{frac}',
                            'Width 1': f'{g_width_A}',
                            'Width 2': f'{g_width_B}',
                            'Alpha 1': f'{g_alpha_A}',
                            'Alpha 2': f'{g_alpha_B}'
                            }
            dict_temp_2 = {
                            'MC Scheme': f"{stratintab}",
                            'Shape': f"{shapeintab2}",
                            'Mean': f'{mean}',
                            'Fraction of Fit': f'{1 - frac}',
                            'Width 1': f'{bg_width_1}',
                            'Width 2': f'{bg_width_2}',
                            'Alpha 1': f'{bg_alpha_1}',
                            'Alpha 2': f'{bg_alpha_2}'
                            }
            df = pd.DataFrame(dict_temp_1, index=[0])
            txt = f"mc_fit_txt/{spec}_1.txt"
            df.to_csv(txt)
            df = pd.DataFrame(dict_temp_2, index=[0])
            txt = f"mc_fit_txt/{spec}_2.txt"
            df.to_csv(txt)

        else:
            dict_temp = {
                            'MC Scheme': f"{id_to_scheme_dict[id_to_spec_dict[spec]]}",
                            'Shape': f"{shape_flag}",
                            'Mean': f'{mean}',
                            'Fraction of Fit': f'{frac}',
                            'Width 1': f'{width_A}',
                            'Width 2': f'{width_B}',
                            'Alpha 1': f'{alpha_1}',
                            'Alpha 2': f'{alpha_2}'
                            }

            df = pd.DataFrame(dict_temp, index=[0])
            txt = f"mc_fit_txt/{spec}.txt"
            df.to_csv(txt)

    file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/latex_tables/mc_fit_txt/*")

    df_1 = pd.DataFrame()
    frames = []

    for file in file_list:
        df_cols_1 = ['MC Scheme','Shape','Mean', 'Width 1', 'Width 2', 'Alpha 1', 'Alpha 2']
        frame = pd.read_csv(file, usecols = df_cols_1)
        print(frame)
        frames.append(frame)

    df_1 = pd.concat(frames)

    caption_line = "Summary of Fit PDFs for MC"
    df_1 = df_1.set_index(['MC Scheme','Shape'])

    label_line = "tab:mcfits"

    order_list = [1,'2a,3a','4a','4b','4c','4d',5,6,'7a','7b','7c','8a','8b','8c',9,10,12,"norm7","norm8"]
    df_1 = df_1.reindex(order_list, level = 0)

    latex_1 = df_1.style.to_latex(position="h", position_float="centering", hrules=True, multirow_align="t", multicol_align="r", column_format='ccccccc', caption = caption_line, label = f"{label_line}_1")

    print(df_1)

    with open(f'tex_files/mcfit_tab.tex','w') as tf:
        tf.write(latex_1)
        # tf.write(latex_2)
def build_mc_rw():

    file_list = glob.glob(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/sys_studies/txt_files/*rw.txt")
    frames = []

    pd.set_option('display.float_format', '{:.2E}'.format)

    for file in file_list:
        df_cols_1 = ['Data Spectrum','MC Scheme','Old Overall Eff','New Overall Eff','Relative Change (Systematic)']
        frame = pd.read_csv(file, usecols = df_cols_1)
        frames.append(frame)

    df = pd.concat(frames)

    caption_line = "Systematic Uncertiees for MC Efficiencis due to Mis-Modeling"
    df = df.set_index(['Data Spectrum','MC Scheme'])

    label_line = "tab:mcrw"

    order_list = ['Z','ZZ','P','ST',"N7",'N8']
    df = df.reindex(order_list, level = 0)

    print(df)


    s = df.style.format('{:.2E}')

    latex_1 = s.to_latex(position="h", position_float="centering", hrules=True, multirow_align="t", multicol_align="r", column_format='ccccccc', caption = caption_line, label = f"{label_line}_1")

    with open(f'tex_files/t_mcrw.tex','w') as tf:
        tf.write(latex_1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--d_window', action='store_true',
                    help='Gen D_WINDOW Table')
    parser.add_argument('--mc_fits', action='store_true',
                    help='Gen mc shape table')
    parser.add_argument('--mc_rw', action='store_true',
                    help='Gen mc_rw Table')

    args = parser.parse_args()
    d_window = args.d_window
    mc_fits = args.mc_fits
    mc_rw = args.mc_rw

    if d_window:
        build_d_window()
    if mc_fits:
        build_mc_shape_table()
    if mc_rw:
        build_mc_rw()
