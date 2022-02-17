import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *
from tabulate import tabulate

inbook = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_eff.xlsx"

Kst_factor = ufloat(2 / 3, 0)

dm = ufloat(0.0938, 0.0016)
dsm = ufloat(0.0539, 0.0015)

d0bar = ufloat(0.0395, 0.00031)
d0bar4 = ufloat(0.0823, 0.0014)

dstmdmpi0 = ufloat(0.323, 0.005)
dstmd0pim = ufloat(0.677, 0.005)

B_B0_norm7 = ufloat(1.31e-03, 7.00e-05)
B_B0_norm8 = ufloat(1.07e-03, 1.10e-04)

#"spec" : [Harris ID, "Shape", "B Mean Guess", "B_Mean_Guess_Window", "B Fit Window"] ...

ids_and_shapes_dp = {

                    "z" : [["01", "G", 5280, 10, 50],
                           ["02", "BG", 5130, 10, 50],
                           ["04", "GEP", 4975, 10, 50]],

                     "zz": [["09", "G", 5280, 10, 50],
                           ["0710", "BG", 5130, 15, 80],
                           ["040812", "BG", 4975, 20, 50]],

                     "p": [["05", "G", 5280, 10, 50],
                           ["020607","BGEP", 5130, 15, 80],
                           ["0408", "BG", 4975, 20, 50]],

                     "m": [["03", "BG", 5130, 15, 80],
                           ["04", "G", 4975, 20, 50]],

                     "st": [["07", "G", 5280, 10, 50],
                           ["0408","BGEP", 5130, 15, 80]],
                     }



def get_mc_shape(dws, spec, fit_strat, slist):

    hid = slist[0]
    shape = slist[1]
    hid_string = f"{spec}_{hid}_{shape}"

    mcws_base = ROOT.TFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_{hid_string}.root")
    mcws = mcws_base.Get(f"fit_{spec}_{hid}")
    mc_pdf = mcws.pdf(f"{spec}_{hid}_{shape}_fit")
    mc_pdf.SetName(f"{spec}_{hid}_fit")

    print (f"{spec}_{hid}_fit")

    dws.Import(mc_pdf)

    if "_c" in fit_strat:
        vars = mcws.allVars()
        for i in vars:
            if i.GetName != "B_DTF_M":
                temp = dws.var(i.GetName())
                temp.setConstant(True)
                print(f"{temp} is constant")

def makeeigenpars(label, parvals, covlists, debug=False, numin=-10, numax=10, expoExtrap=False):

    vpars = []
    vnupars = []

    n = len(parvals)

    if n == 1:
        ## handle the case we were given 1 parameter gracefully
        eigval = ROOT.TVectorD(1)
        eigvec = ROOT.TMatrixD(1, 1)
        eigval[0] = covlists[0][0]
        eigvec[0][0] = 1.0
    else:

        ## covlist should be a simple nested list with covariance
        covmat = ROOT.TMatrixD(n, n)
        for i in range(n):
            for j in range(n):
                covmat[i][j] = covlists[i][j]

        if debug:
            print(label, "covmat =")
            covmat.Print()

        eigval = ROOT.TVectorD()
        eigvec = covmat.EigenVectors(eigval)

    if debug:
        print(label, " eigenvals and vecs:")
        eigval.Print()
        eigvec.Print()

    for i in range(n):
        vnupars.append(
            ROOT.RooRealVar(label + "_nu{}".format(i), label, 0, numin, numax)
        )
        vnupars[-1].setError(1.0)

    vnupars_flist = ROOT.RooArgList()
    for j in vnupars:
        vnupars_flist.add(j)

    for i in range(n):

        if expoExtrap:
            expr = "{}".format(parvals[i])
            for j in range(n):
                expr += "*pow( ({v0} + {l}*{e})/{v0}, {label}_nu{j})".format(
                    v0=parvals[i], l=sqrt(eigval[j]), e=eigvec[i][j], label=label, j=j
                )
        else:
            expr = "{}".format(parvals[i])
            for j in range(n):
                expr += " + {}_nu{}*{}*{}".format(
                    label, j, sqrt(eigval[j]), eigvec[i][j]
                )
        if debug:
            print(expr)

        vpars.append(
            ROOT.RooFormulaVar(label + "_par{}".format(i), label, expr, vnupars_flist)
        )

    # vpars[-1].Print()
    return vpars, vnupars

def get_norm_param(norm_id, trigger_flag, year):

    norm_id = norm_id.split("_")[0]
    norm_base_file = ROOT.TFile(f"{analysis_path}/normalization_fit/base_norm_files/{norm_id}_2016_nTaT.root")
    norm_ws = norm_base_file.Get(f"{norm_id}")
    n_base = norm_ws.var(f"n_{norm_id}_signal")
    n_val = ufloat(n_base.getValV(), n_base.getError())/6

    #FIX /6 Once finishd new dataset

    if "norm8" in norm_id:
        cfactor =  (B_B0_norm8 * dm * d0bar4)
    if "norm7" in norm_id:
        cfactor = (B_B0_norm7 * d0bar * d0bar4)

    return (n_val, cfactor)

def get_eff(id, trigger_flag, year):
    df = pandas.read_excel(inbook, sheet_name=f"MC_eff_{trigger_flag}")
    n = df.loc[df["Front_ID"] == id, f"total_booteff_{year}"].values[0]
    s = df.loc[df["Front_ID"] == id, f"err_total_booteff_{year}"].values[0]
    eff = ufloat(n, s)
    return eff

df_n_rows_list = []
df_tre_rows_list = []

df_bfs_list = []

df_index_list = []
df_tre_index_list = []

def get_a_factor(id, norm_id, signal_d_bfs):

    eff = 0
    dict = {}
    dict_tre = {}
    final_bf_factor= 0

    for year in ["2016", "2017", "2018"]:
        for trigger_flag in ["T", "nTaT"]:

            #Norm_Yield is the yield of the noralization mode
            #b_norm_factors is the product of theBranching Fractions for B->D1normD2normK * D1_norm * D2_norm
            Norm_Yield, b_norm_factors = get_norm_param(norm_id, trigger_flag, year)

            #efs is the signal mode efficency
            efs = get_eff(id, trigger_flag, year)

            #efn is the normalizatopn mode efficency
            efn = get_eff(norm_id, trigger_flag, year)

            dict[f"{year}_{trigger_flag}"] = f"{efs:.3e}"

            #This is the final number for a given year, trigger condition
            #signal_bfs is the produut of the signal branching fractions for D1_sig*D2_sig*Kst0
            signal_bfs = signal_d_bfs*Kst_factor

            final_base_factor = Norm_Yield*(signal_bfs/b_norm_factors)*(efs/efn)

            dict_tre[f"{year}_{trigger_flag}"] = f"{final_base_factor:.3e}"

            #this is the sum of all final_base_factors
            final_bf_factor = final_bf_factor + final_base_factor

    df_bfs_list.append(signal_bfs/b_norm_factors)
    df_n_rows_list.append(dict)
    df_tre_rows_list.append(dict_tre)
    df_index_list.append(id)
    df_tre_index_list.append(id)

    #    add norm modes
    if id == "16_s":
        for norm_id in ["norm7_norm7", "norm8_norm8"]:
            dict = {}
            for year in ["2016", "2017", "2018"]:
                for trigger_flag in ["T", "nTaT"]:
                    efn = get_eff(norm_id, trigger_flag, year)
                    dict[f"{year}_{trigger_flag}"] = f"{efn:.3e}"
            df_n_rows_list.append(dict)
            df_index_list.append(norm_id)

    return final_bf_factor

def build_fit(dws, run_name, gc_onflag, fit_strat, specs, ids_and_shapes):
    #Build n params
    params = []

    a_01_z = get_a_factor("01_z", "norm8_norm8", dm*dm)
    a_02_z = get_a_factor("02_z", "norm8_norm8", dstmdmpi0 * dm * dm)
    a_03_z = get_a_factor("02_z", "norm8_norm8", dstmdmpi0 * dm * dm)
    a_04_z = get_a_factor("04_z", "norm8_norm8", dstmdmpi0 * dstmdmpi0 * dm * dm)

    params = params + [a_01_z, a_02_z, a_03_z, a_04_z]

    a_04_zz = get_a_factor("04_zz", "norm7_norm7", dstmd0pim * dstmd0pim * d0bar * d0bar)

    a_07_zz = get_a_factor("07_zz", "norm7_norm7", dstmd0pim * d0bar * d0bar)
    a_08_zz = get_a_factor("08_zz", "norm7_norm7", dstmd0pim * d0bar * d0bar)

    a_09_zz = get_a_factor("09_zz", "norm7_norm7", d0bar * d0bar)
    a_10_zz = get_a_factor("10_zz", "norm7_norm7", d0bar * d0bar)
    a_12_zz = get_a_factor("12_zz", "norm7_norm7", d0bar * d0bar)

    params = params + [a_04_zz, a_07_zz, a_08_zz, a_09_zz, a_10_zz, a_12_zz]

    a_02_p = get_a_factor("02_p", "norm7_norm7", dstmd0pim * dm * d0bar)
    a_04_p = get_a_factor("04_p", "norm7_norm7", dstmd0pim * dstmdmpi0 * d0bar * dm)

    a_05_p = get_a_factor("05_p", "norm7_norm7", dm * d0bar)
    a_06_p = get_a_factor("06_p", "norm7_norm7", dm * d0bar)

    a_07_p = get_a_factor("07_p", "norm7_norm7", d0bar * dstmdmpi0 * dm)
    a_08_p = get_a_factor("08_p", "norm7_norm7", d0bar * dstmdmpi0 * dm)

    params = params + [a_02_p, a_04_p, a_05_p, a_06_p, a_07_p, a_08_p]

    a_03_m = get_a_factor("03_m", "norm7_norm7", dstmd0pim * d0bar * dm)
    a_04_m = get_a_factor("04_m", "norm7_norm7", dstmdmpi0 * dstmd0pim * dm * d0bar)

    params = params + [a_03_m, a_04_m]

    a_04_st = get_a_factor("04_st", "norm7_norm7", dstmd0pim * dstmd0pim * d0bar * d0bar)
    a_07_st = get_a_factor("07_st", "norm7_norm7", dstmd0pim * d0bar * d0bar)
    a_08_st = get_a_factor("08_st", "norm7_norm7", dstmd0pim * d0bar * d0bar)

    params = params + [a_04_st, a_07_st, a_08_st]

    # a_13_s = get_a_factor("13_s", "norm8_norm8", dsm * dm)
    # a_14_s = get_a_factor("14_s", "norm8_norm8", dsm * dm)
    # a_15_s = get_a_factor("15_s", "norm8_norm8", dstmdmpi0 * dsm * dm)
    # a_16_s = get_a_factor("16_s", "norm8_norm8", dstmdmpi0 * dsm * dm)
    #
    # params = params + [a_13_s, a_14_s, a_15_s, a_16_s]

    df_bfs_dict = {"Known BF Product": df_bfs_list}
    df_bfs = pandas.DataFrame(df_bfs_dict, index = df_tre_index_list)
    # for p in params:
    df_params_dict = {"Final Parameter": params}
    df_params = pandas.DataFrame(df_params_dict, index = df_tre_index_list)

    # #params before multiplied by rel branching fractions
    df_effs = pandas.DataFrame(df_n_rows_list, index = df_index_list)
    df_prebfs = pandas.DataFrame(df_tre_rows_list, index = df_tre_index_list)

    original_stdout = sys.stdout # Save a reference to the original standard output

    with open('tables_eff.md', 'w') as f:
        sys.stdout = f # Change the standard output to the file we created.
        print(df_effs.to_markdown())
        print("\n")
        print(df_prebfs.to_markdown())
        print("\n")
        print(df_bfs.to_markdown())
        print("\n")
        print(df_params.to_markdown())
        sys.stdout = original_stdout
    #
    c = covariance_matrix(params)
    test_a, test_b = makeeigenpars(run_name, [p.n for p in params], c, debug=False)
    for a in test_a:
        if gc_onflag == 1:
            dws.Import(a)
        if gc_onflag == 0:
            new_a = a.evaluate()
            new_a_name = a.GetName()
            dws.factory(f"{new_a_name}[{new_a}]")
            tv = dws.var(f"{new_a_name}")
            tv.Print()

    if gc_onflag == 1:
        for i in range(0, len(params)):
            dws.factory(f"Gaussian::{run_name}_nu{i}_g({run_name}_nu{i}, 0, 1)")

    for i in range(1,17):
        dws.factory(f"bf_{i}[0.01,0,1]")

    for spec in specs:
        get_shapes_bkg(spec, "Exponential", dws)
        for s_list in ids_and_shapes[spec]:
            get_mc_shape(dws, spec, fit_strat, s_list)

    dws.factory(f"expr::n_01_z('bf_1*{run_name}_par0', bf_1, {run_name}_par0)")
    dws.factory(f"expr::n_02_z('bf_2*{run_name}_par1', bf_2, {run_name}_par1)")
    dws.factory(f"expr::n_03_z('bf_3*{run_name}_par2', bf_3, {run_name}_par2)")
    dws.factory(f"expr::n_04_z('bf_4*{run_name}_par3', bf_4, {run_name}_par3)")
    dws.factory(f"expr::n_0203_z('n_02_z + n_03_z', n_02_z, n_03_z)")

    dws.factory(f"expr::n_04_zz('bf_4*{run_name}_par4', bf_4, {run_name}_par4)")
    dws.factory(f"expr::n_07_zz('bf_7*{run_name}_par5', bf_7, {run_name}_par5)")
    dws.factory(f"expr::n_08_zz('bf_8*{run_name}_par6', bf_8, {run_name}_par6)")
    dws.factory(f"expr::n_09_zz('bf_9*{run_name}_par7', bf_9, {run_name}_par7)")
    dws.factory(f"expr::n_10_zz('bf_10*{run_name}_par8', bf_10, {run_name}_par8)")
    dws.factory(f"expr::n_12_zz('bf_12*{run_name}_par9', bf_12, {run_name}_par9)")

    dws.factory("expr::n_0710_zz('n_07_zz + n_10_zz', n_07_zz, n_10_zz)")
    dws.factory("expr::n_040812_zz('n_04_zz + n_08_zz + n_12_zz',n_04_zz,n_08_zz,n_12_zz)")

    dws.factory(f"expr::n_02_p('bf_2*{run_name}_par10', bf_2, {run_name}_par10)")
    dws.factory(f"expr::n_04_p('bf_4*{run_name}_par11', bf_4, {run_name}_par11)")
    dws.factory(f"expr::n_05_p('bf_5*{run_name}_par12', bf_5, {run_name}_par12)")
    dws.factory(f"expr::n_06_p('bf_6*{run_name}_par13', bf_6, {run_name}_par13)")
    dws.factory(f"expr::n_07_p('bf_7*{run_name}_par14', bf_7, {run_name}_par14)")
    dws.factory(f"expr::n_08_p('bf_8*{run_name}_par15', bf_8, {run_name}_par15)")

    dws.factory(f"expr::n_020607_p('n_02_p + n_06_p + n_07_p', n_02_p, n_06_p, n_07_p)")
    dws.factory(f"expr::n_0408_p('n_04_p + n_08_p', n_04_p, n_08_p)")

    dws.factory(f"expr::n_03_m('bf_3*{run_name}_par16', bf_3, {run_name}_par16)")
    dws.factory(f"expr::n_04_m('bf_4*{run_name}_par17', bf_4, {run_name}_par17)")

    dws.factory(f"expr::n_04_st('bf_4*{run_name}_par18', bf_4, {run_name}_par18)")
    dws.factory(f"expr::n_07_st('bf_7*{run_name}_par19', bf_7, {run_name}_par19)")
    dws.factory(f"expr::n_08_st('bf_8*{run_name}_par20', bf_8, {run_name}_par20)")
    dws.factory("expr::n_0408_st('n_04_st + n_08_st', n_04_st, n_08_st)")

    # dws.factory(f"expr::n_13_s('bf_13*{run_name}_par21', bf_13, {run_name}_par21)")
    # dws.factory(f"expr::n_14_s('bf_14*{run_name}_par22', bf_14, {run_name}_par22)")
    # dws.factory(f"expr::n_15_s('bf_15*{run_name}_par23', bf_15, {run_name}_par23)")
    # dws.factory(f"expr::n_16_s('bf_16*{run_name}_par24', bf_16, {run_name}_par24)")


    if "dp_f" in fit_strat:
        dws.factory("SUM::z_spectrum_all_fit(n_01_z*z_01_fit,n_0203_z*z_02_fit,n_04_z*z_04_fit,n_z_bkg[100,0,10000]*z_spectrum_bkg)")
        dws.factory("SUM::zz_spectrum_all_fit(n_09_zz*zz_09_fit,n_0710_zz*zz_0710_fit,n_040812_zz*zz_040812_fit,n_zz_bkg[100,0,100000]*zz_spectrum_bkg)")
        dws.factory("SUM::p_spectrum_all_fit(n_05_p*p_05_fit,n_020607_p*p_020607_fit,n_0408_p*p_0408_fit,n_p_bkg[100,0,10000]*p_spectrum_bkg)")
        dws.factory("SUM::m_spectrum_all_fit(n_03_m*m_03_fit,n_04_m*m_04_fit,n_m_bkg[100,0,10000]*m_spectrum_bkg)")
        dws.factory("SUM::st_spectrum_all_fit(n_07_st*st_07_fit,n_0408_st*st_0408_fit,n_st_bkg[100,0,1000]*st_spectrum_bkg)")
        # dws.factory("SUM::s_spectrum_all_fit(n_13_s*s_13_fit,n_14_s*s_14_fit,n_15_s*s_15_fit,n_16_s*s_16_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")

    if fit_strat == "dp_f_nn":
        for i in range(1, 18):
            dws.factory(f"nny_{i}[500,0,10000]")
        dws.factory("SUM::z_spectrum_all_fit(nny_1*z_01_fit,nny_2*z_02_fit,nny_3*z_04_fit,n_z_bkg[100,0,10000]*z_spectrum_bkg)")
        dws.factory("SUM::zz_spectrum_all_fit(nny_4*zz_09_fit,nny_5*zz_0710_fit,nny_6*zz_040812_fit,n_zz_bkg[100,0,100000]*zz_spectrum_bkg)")
        dws.factory("SUM::p_spectrum_all_fit(nny_7*p_05_fit,nny_8*p_020607_fit,nny_9*p_0408_fit,n_p_bkg[100,0,10000]*p_spectrum_bkg)")
        dws.factory("SUM::m_spectrum_all_fit(nny_10*m_03_fit,nny_11*m_04_fit,n_m_bkg[100,0,10000]*m_spectrum_bkg)")
        dws.factory("SUM::st_spectrum_all_fit(nny_12*st_07_fit,nny_13*st_0408_fit,n_st_bkg[100,0,1000]*st_spectrum_bkg)")
        # dws.factory("SUM::s_spectrum_all_fit(nny_14*s0_fit,nny_15*sL1_fit,nny_16*sR1_fit,nny_17*s2_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")

def build_data_ws(run_name, gc_onflag, fit_strat, specs):

    dws = ROOT.RooWorkspace(run_name)
    dws.factory(f"B_DTF_M[{bmin},{bmax}]")

    ids_and_shapes = ids_and_shapes_dp

    ###########################
    build_fit(dws, run_name, gc_onflag, fit_strat, specs, ids_and_shapes)
    ###########################

    b_dtf_m = dws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_dtf_m)
    ############################

    all_cats = ROOT.RooCategory("all_cats", "all_cats")

    all_cats.defineType("z_spectrum")
    z_data_file = ROOT.TFile(f"{root_basepath}DATA/z_none_ToT_spectrum_filtered.root")
    z_tree = z_data_file.Get("DecayTreeTuple")
    z_data_set = ROOT.RooDataSet("z_data", "z_data", z_tree, data_args)
    all_cats.defineType("zz_spectrum")
    zz_data_file = ROOT.TFile(f"{root_basepath}DATA/zz_none_ToT_spectrum_filtered.root")
    zz_tree = zz_data_file.Get("DecayTreeTuple")
    zz_data_set = ROOT.RooDataSet("zz_data", "zz_data", zz_tree, data_args)
    all_cats.defineType("p_spectrum")
    p_data_file = ROOT.TFile(f"{root_basepath}DATA/p_none_ToT_spectrum_filtered.root")
    p_tree = p_data_file.Get("DecayTreeTuple")
    p_data_set = ROOT.RooDataSet("p_data", "p_data", p_tree, data_args)
    all_cats.defineType("m_spectrum")
    m_data_file = ROOT.TFile(f"{root_basepath}DATA/m_none_ToT_spectrum_filtered.root")
    m_tree = m_data_file.Get("DecayTreeTuple")
    m_data_set = ROOT.RooDataSet("m_data", "m_data", m_tree, data_args)
    all_cats.defineType("st_spectrum")
    st_data_file = ROOT.TFile(f"{root_basepath}DATA/st_none_ToT_spectrum_filtered.root")
    st_tree = st_data_file.Get("DecayTreeTuple")
    st_data_set = ROOT.RooDataSet("st_data", "st_data", st_tree, data_args)
    all_cats.defineType("s_spectrum")
    s_data_file = ROOT.TFile(f"{root_basepath}DATA/s_none_ToT_spectrum_filtered.root")
    s_tree = s_data_file.Get("DecayTreeTuple")
    s_data_set = ROOT.RooDataSet("s_data", "s_data", s_tree, data_args)

    all_data_sets = ROOT.RooDataSet(
        "all_data_sets",
        "all_data_sets",
        data_args,
        ROOT.RooFit.Index(all_cats),
        ROOT.RooFit.Import("z_spectrum", z_data_set),
        ROOT.RooFit.Import("zz_spectrum", zz_data_set),
        ROOT.RooFit.Import("p_spectrum", p_data_set),
        ROOT.RooFit.Import("m_spectrum", m_data_set),
        ROOT.RooFit.Import("st_spectrum", st_data_set),
        ROOT.RooFit.Import("s_spectrum", s_data_set),
    )

    getattr(dws, "import")(all_cats)
    getattr(dws, "import")(all_data_sets)

    dws.writeToFile(f"{analysis_path}/data_fit/base_data_files/{run_name}.root")
    print(f"Wrote dws to: {analysis_path}/data_fit/base_data_files/{run_name}.root")
    dws.Print()

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--run_name')
parser.add_argument('--gc_onflag', type=int)
parser.add_argument('--fit_strat')
parser.add_argument('--specs', nargs='+')

args = parser.parse_args()
run_name = args.run_name
gc_onflag = args.gc_onflag
fit_strat = args.fit_strat
specs  = args.specs

build_data_ws(run_name, gc_onflag, fit_strat, specs)
