import ROOT as ROOT
from Analysis_2021.essentials import get_shapes_bkg, get_free_shapes, save_png
import glob as glob
analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"
from uncertainties import ufloat
bbins = 100
bmin = 4800
bmax = 5600

#########################################
#"spec" : [Harris ID, "Shape", "B Mean Guess", "B_Mean_Guess_Window", "B Fit Window"] ...
ids_and_shapes_dp = {

                    "Z_m_p" : [["01", "DG", 5280, 10, 50],
                           ["02", "BGEP", 5130, 10, 50],
                           ["04", "G", 4985, 15, 50]],

                     "Z_z_z": [["09", "DG", 5280, 10, 50],
                           ["0710", "BG", 5130, 15, 80],
                           ["040812", "BG", 4975, 20, 50]],

                     "P_z_p": [["05", "G", 5280, 10, 50],
                           ["020607","BGEP", 5130, 15, 80],
                           ["0408", "BG", 4975, 20, 50]],

                     "M_m_z": [["03", "BG", 5130, 15, 80],
                           ["04", "BG", 4975, 20, 50]],

                     "P_z_pst": [["07", "BG", 5280, 10, 50],
                           ["0408","DG", 5130, 15, 80]],

                     "Zs_sm_p" : [["13","G", 5370, 10, 50],
                                 ["14","G", 5225, 15, 80],
                                 ["15","G", 5225, 15, 80],
                                 ["16","BG", 5075, 15, 80]],
                     }

def get_mc_shape(dws, spec, fit_strat, slist):

    hid = slist[0]
    shape = slist[1]
    hid_string = f"{hid}_{spec}_{shape}"

    print(hid_string)
    if spec == "M_m_z" or spec == "P_z_pst" or spec == "Z_m_p" or spec == "Z_z_z" or spec == "Zs_sm_p":
        mcws_base = ROOT.TFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_MC_{hid_string}.root")
        mcws = mcws_base.Get(f"fit_MC_{hid}_{spec}")
        mc_pdf = mcws.pdf(f"MC_{hid}_{spec}_{shape}_fit")
    else:
        mcws_base = ROOT.TFile(f"{analysis_path}/mc_fit/fit_mc_files/fit_{hid_string}.root")
        mcws = mcws_base.Get(f"fit_{hid}_{spec}")
        mc_pdf = mcws.pdf(f"{hid}_{spec}_{shape}_fit")
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
    norm_base_file = ROOT.TFile(f"{analysis_path}/normalization_fit/base_norm_files/{norm_id}_{year}_{trigger_flag}.root")
    norm_ws = norm_base_file.Get(f"{norm_id}")
    n_base = norm_ws.var(f"n_{norm_id}_signal")
    n_val = ufloat(n_base.getValV(), n_base.getError())

    if "norm8" in norm_id:
        cfactor =  (B_B0_norm8 * dm * d0bar4)
        print(year, n_val)
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

def MakeSWeights(outfilename, outtreename, data, model, yields):

    """Determine s-weights from fit.

    arguments:
    outfilename -- name of .root file to create with `outtreename`
    outtreename -- name of TTree with s-weights to save in `outfilename`
    data -- RooDataSet to which `model` was fitted
    model -- fitted RooAbsPdf
    yields -- RooArgList of RooRealVars extracted from fitting `model` to `data`
    """
    from array import array

    print(f"using data '{data.GetName()}'")
    print(f"using model '{model.GetName()}'")
    print(f"using yields '{[x.GetName() for x in yields]}'")

    sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields)

    print("Check SWeights:")
    # for y in yields:
    #     oval = y.getVal()
    #     sval = sData.GetYieldFromSWeight(y.GetName())
    #     print(f"Yield of {y.GetName()} is {oval}")
    #     print(f"from sWeights it is {sval}")
    #     if not (0.9995 < oval / sval < 1.0005):
    #         raise Exception("sWeighted yield should match")
    # for i in range(10):
    #     for y in yields:
    #         print(f"    {y.GetName()} Weight {sData.GetSWeight(i, y.GetName())}")
    #     totw = sData.GetSumOfEventSWeight(i)
    #     print(f"Total Weight {totw}")
    #     if not (0.9995 < totw < 1.0005):
    #         raise Exception("sum of sWeight should be 1")
    swnames = sorted([f"{x.GetName()}_sw" for x in yields])
    print(f"weights: {swnames}")
    # create output file
    nf = ROOT.TFile.Open(outfilename, "recreate")
    # create directory hierarchy
    nd = nf
    if len(outtreename.split("/")) > 1:
        for d in outtreename.split("/")[:-1]:
            nd = nd.mkdir(d)

    nd.cd()
    # create output TTree
    nt = ROOT.TTree(outtreename.split("/")[-1], outtreename.split("/")[-1])
    # declare output branches
    swvals = [array("f", [0]) for x in swnames]
    for val, nm in zip(swvals, swnames):
        nt.Branch(nm, val, f"{nm}/F")
    # loop data
    for i in range(data.numEntries()):
        # get vars
        swvars = sorted(
            [x for x in data.get(i) if x.GetName() in swnames],
            key=lambda x: x.GetName(),
        )
        assert [x.GetName() for x in swvars] == swnames  # check sorting worked
        # set values
        for val, var in zip(swvals, swvars):
            val[0] = var.getVal()
        # fill values
        nt.Fill()
    nt.Write()
    nf.Close()

    # tltxt = ROOT.TLatex()

def get_a_factor(id, norm_id, signal_d_bfs, genfix = "None"):

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
            if genfix == "None":
                efs = get_eff(id, trigger_flag, year)

            if genfix != "None":
                pi_num = get_eff("04_Z_mst_p", trigger_flag, year)
                pi_denom = get_eff("04_P_z_p", trigger_flag, year)
                if "Z_z_z" in id:
                    efs = get_eff(genfix, trigger_flag, year)*(1 - (pi_num/pi_denom))
                if "P_z_pst" in id:
                    efs = get_eff(genfix, trigger_flag, year)*(pi_num/pi_denom)

            #efn is the normalizatopn mode efficency
            efn = get_eff(norm_id, trigger_flag, year)

            dict[f"{year}_{trigger_flag}"] = f"{efs:.3e}"

            #This is the final number for a given year, trigger condition
            #signal_bfs is the produut of the signal branching fractions for D1_sig*D2_sig*Kst0
            signal_bfs = signal_d_bfs*Kst_factor

            final_base_factor = Norm_Yield*(signal_bfs/b_norm_factors)*(efs/efn)

            # dict_tre[f"{year}_{trigger_flag}"] = f"{signal_bfs}"
            dict_tre[f"{year}_{trigger_flag}"] = f"{final_base_factor:.3e}"

            #this is the sum of all final_base_factors
            final_bf_factor = final_bf_factor + final_base_factor

    # df_bfs_list.append(signal_bfs/b_norm_factors)
    df_bfs_list.append(signal_bfs)
    df_n_rows_list.append(dict)
    df_tre_rows_list.append(dict_tre)
    df_index_list.append(id)
    df_tre_index_list.append(id)

    #  add norm modes
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

def build_nn_sw_fit(dws, run_name, gc_onflag, fit_strat, specs, ids_and_shapes):

    for spec in specs:
        get_shapes_bkg(spec, "Exponential", dws)
        for s_list in ids_and_shapes[spec]:
            get_mc_shape(dws, spec, fit_strat, s_list)

    for i in range(1, 18):
        dws.factory(f"nny_{i}[500,0,10000]")

    dws.factory("SUM::Z_m_p_spectrum_all_fit(nny_1*Z_m_p_01_fit,nny_2*Z_m_p_02_fit,nny_3*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
    # dws.factory("SUM::Z_z_z_spectrum_all_fit(nny_4*Z_z_z_09_fit,nny_5*Z_z_z_0710_fit,nny_6*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
    # dws.factory("SUM::P_z_p_spectrum_all_fit(nny_7*P_z_p_05_fit,nny_8*P_z_p_020607_fit,nny_9*P_z_p_0408_fit,n_P_z_P_z_p_bkg[100,0,10000]*P_z_p_spectrum_bkg)")
    # dws.factory("SUM::M_m_z_spectrum_all_fit(nny_10*M_m_z_03_fit,nny_11*M_m_z_04_fit,n_M_m_z_bkg[100,0,10000]*M_m_z_spectrum_bkg)")
    # dws.factory("SUM::P_z_pst_spectrum_all_fit(nny_12*P_z_pst_07_fit,nny_13*P_z_pst_0408_fit,n_P_z_st_bkg[100,0,1000]*P_z_pst_spectrum_bkg)")
    # dws.factory("SUM::s_spectrum_all_fit(nny_14*s0_fit,nny_15*sL1_fit,nny_16*sR1_fit,nny_17*s2_fit,n_s_bkg[100,0,1000]*s_spectrum_bkg)")
    dws.Print()

def build_fit(dws, run_name, gc_onflag, fit_strat, specs, ids_and_shapes):
    #Build n params
    params = []


    a_01_z = get_a_factor("01_Z_m_p", "norm8_norm8", dm*dm)
    a_02_z = get_a_factor("02_Z_m_p", "norm8_norm8", dstmdmpi0 * dm * dm)
    a_03_z = get_a_factor("02_Z_m_p", "norm8_norm8", dstmdmpi0 * dm * dm)
    a_04_z = get_a_factor("04_Z_m_p", "norm8_norm8", dstmdmpi0 * dstmdmpi0 * dm * dm)

    params = params + [a_01_z, a_02_z, a_03_z, a_04_z]

    a_04_zz = get_a_factor("04_Z_z_z", "norm7_norm7", dstmd0pim * dstmd0pim * d0bar * d0bar, "12_Z_z_z")

    a_07_zz = get_a_factor("07_Z_z_z", "norm7_norm7", dstmd0pim * d0bar * d0bar, "10_Z_z_z")
    a_08_zz = get_a_factor("08_Z_z_z", "norm7_norm7", dstmd0pim * d0bar * d0bar, "12_Z_z_z")

    a_09_zz = get_a_factor("09_Z_z_z", "norm7_norm7", d0bar * d0bar)
    a_10_zz = get_a_factor("10_Z_z_z", "norm7_norm7", d0bar * d0bar)
    a_12_zz = get_a_factor("12_Z_z_z", "norm7_norm7", d0bar * d0bar)

    params = params + [a_04_zz, a_07_zz, a_08_zz, a_09_zz, a_10_zz, a_12_zz]

    a_02_p = get_a_factor("02_P_z_p", "norm7_norm7", dstmd0pim * dm * d0bar)
    a_04_p = get_a_factor("04_P_z_p", "norm7_norm7", dstmd0pim * dstmdmpi0 * d0bar * dm)

    a_05_p = get_a_factor("05_P_z_p", "norm7_norm7", dm * d0bar)
    a_06_p = get_a_factor("06_P_z_p", "norm7_norm7", dm * d0bar)

    a_07_p = get_a_factor("07_P_z_p", "norm7_norm7", d0bar * dstmdmpi0 * dm)
    a_08_p = get_a_factor("08_P_z_p", "norm7_norm7", d0bar * dstmdmpi0 * dm)

    params = params + [a_02_p, a_04_p, a_05_p, a_06_p, a_07_p, a_08_p]

    a_03_m = get_a_factor("03_M_m_z", "norm7_norm7", dstmd0pim * d0bar * dm)
    a_04_m = get_a_factor("04_M_m_z", "norm7_norm7", dstmdmpi0 * dstmd0pim * dm * d0bar)

    params = params + [a_03_m, a_04_m]

    a_04_st = get_a_factor("04_P_z_pst", "norm7_norm7", dstmd0pim * dstmd0pim * d0bar * d0bar, "12_Z_z_z")
    a_07_st = get_a_factor("07_P_z_pst", "norm7_norm7", dstmd0pim * d0bar * d0bar, "10_Z_z_z")
    a_08_st = get_a_factor("08_P_z_pst", "norm7_norm7", dstmd0pim * d0bar * d0bar,"12_Z_z_z")

    params = params + [a_04_st, a_07_st, a_08_st]

    a_13_s = get_a_factor("13_Zs_sm_p", "norm8_norm8", dsm * dm)
    a_14_s = get_a_factor("14_Zs_sm_p", "norm8_norm8", dsm * dm)
    a_15_s = get_a_factor("15_Zs_sm_p", "norm8_norm8", dstmdmpi0 * dsm * dm)
    a_16_s = get_a_factor("16_Zs_sm_p", "norm8_norm8", dstmdmpi0 * dsm * dm)

    params = params + [a_13_s, a_14_s, a_15_s, a_16_s]

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
        dws.factory(f"bf_{i}[0.001,0,0.01]")

    for spec in specs:
        get_shapes_bkg(spec, "Exponential", dws)
        for s_list in ids_and_shapes[spec]:
            get_mc_shape(dws, spec, fit_strat, s_list)

    dws.factory(f"expr::n_01_Z_m_p('bf_1*{run_name}_par0', bf_1, {run_name}_par0)")
    dws.factory(f"expr::n_02_Z_m_p('bf_2*{run_name}_par1', bf_2, {run_name}_par1)")
    dws.factory(f"expr::n_03_Z_m_p('bf_3*{run_name}_par2', bf_3, {run_name}_par2)")
    dws.factory(f"expr::n_04_Z_m_p('bf_4*{run_name}_par3', bf_4, {run_name}_par3)")
    dws.factory(f"expr::n_0203_Z_m_p('n_02_Z_m_p + n_03_Z_m_p', n_02_Z_m_p, n_03_Z_m_p)")

    dws.factory(f"expr::n_04_Z_z_z('bf_4*{run_name}_par4', bf_4, {run_name}_par4)")
    dws.factory(f"expr::n_07_Z_z_z('bf_7*{run_name}_par5', bf_7, {run_name}_par5)")
    dws.factory(f"expr::n_08_Z_z_z('bf_8*{run_name}_par6', bf_8, {run_name}_par6)")
    dws.factory(f"expr::n_09_Z_z_z('bf_9*{run_name}_par7', bf_9, {run_name}_par7)")
    dws.factory(f"expr::n_10_Z_z_z('bf_10*{run_name}_par8', bf_10, {run_name}_par8)")
    dws.factory(f"expr::n_12_Z_z_z('bf_12*{run_name}_par9', bf_12, {run_name}_par9)")

    dws.factory("expr::n_0710_Z_z_z('n_07_Z_z_z + n_10_Z_z_z', n_07_Z_z_z, n_10_Z_z_z)")
    dws.factory("expr::n_040812_Z_z_z('n_04_Z_z_z + n_08_Z_z_z + n_12_Z_z_z',n_04_Z_z_z,n_08_Z_z_z,n_12_Z_z_z)")

    dws.factory(f"expr::n_02_P_z_p('bf_2*{run_name}_par10', bf_2, {run_name}_par10)")
    dws.factory(f"expr::n_04_P_z_p('bf_4*{run_name}_par11', bf_4, {run_name}_par11)")
    dws.factory(f"expr::n_05_P_z_p('bf_5*{run_name}_par12', bf_5, {run_name}_par12)")
    dws.factory(f"expr::n_06_P_z_p('bf_6*{run_name}_par13', bf_6, {run_name}_par13)")
    dws.factory(f"expr::n_07_P_z_p('bf_7*{run_name}_par14', bf_7, {run_name}_par14)")
    dws.factory(f"expr::n_08_P_z_p('bf_8*{run_name}_par15', bf_8, {run_name}_par15)")

    dws.factory(f"expr::n_020607_P_z_p('n_02_P_z_p + n_06_P_z_p + n_07_P_z_p', n_02_P_z_p, n_06_P_z_p, n_07_P_z_p)")
    dws.factory(f"expr::n_0408_P_z_p('n_04_P_z_p + n_08_P_z_p', n_04_P_z_p, n_08_P_z_p)")

    dws.factory(f"expr::n_03_M_m_z('bf_3*{run_name}_par16', bf_3, {run_name}_par16)")
    dws.factory(f"expr::n_04_M_m_z('bf_4*{run_name}_par17', bf_4, {run_name}_par17)")

    dws.factory(f"expr::n_04_P_z_pst('bf_4*{run_name}_par18', bf_4, {run_name}_par18)")
    dws.factory(f"expr::n_07_P_z_pst('bf_7*{run_name}_par19', bf_7, {run_name}_par19)")
    dws.factory(f"expr::n_08_P_z_pst('bf_8*{run_name}_par20', bf_8, {run_name}_par20)")
    dws.factory("expr::n_0408_P_z_pst('n_04_P_z_pst + n_08_P_z_pst', n_04_P_z_pst, n_08_P_z_pst)")

    dws.factory(f"expr::n_13_Zs_sm_p('bf_13*{run_name}_par21', bf_13, {run_name}_par21)")
    dws.factory(f"expr::n_14_Zs_sm_p('bf_14*{run_name}_par22', bf_14, {run_name}_par22)")
    dws.factory(f"expr::n_15_Zs_sm_p('bf_15*{run_name}_par23', bf_15, {run_name}_par23)")
    dws.factory(f"expr::n_16_Zs_sm_p('bf_16*{run_name}_par24', bf_16, {run_name}_par24)")

    # if fit_strat == "dp_f":
        # for i in range(1, 18):
        #     dws.factory(f"nny_{i}[500,0,10000]")
        # if len(specs) == 1:
        #     dws.factory("SUM::Z_m_p_spectrum_all_fit(n_01_Z_m_p*Z_m_p_01_fit,n_02_Z_m_p*Z_m_p_02_fit,n_04_Z_m_p*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
        # if len(specs) == 2:
        #     dws.factory("SUM::Z_m_p_spectrum_all_fit(n_01_Z_m_p*Z_m_p_01_fit,n_02_Z_m_p*Z_m_p_02_fit,n_04_Z_m_p*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
        #     dws.factory("SUM::Z_z_z_spectrum_all_fit(n_09_Z_z_z*Z_z_z_09_fit,nny_3*Z_z_z_0710_fit,nny_3*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
        # if len(specs) == 3:
        #     dws.factory("SUM::Z_m_p_spectrum_all_fit(n_01_Z_m_p*Z_m_p_01_fit,n_02_Z_m_p*Z_m_p_02_fit,n_04_Z_m_p*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
        #     dws.factory("SUM::Z_z_z_spectrum_all_fit(n_09_Z_z_z*Z_z_z_09_fit,nny_3*Z_z_z_0710_fit,nny_3*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
        #     dws.factory("SUM::P_z_p_spectrum_all_fit(n_05_P_z_p*P_z_p_05_fit,nny_5*P_z_p_020607_fit,nny_6*P_z_p_0408_fit,n_P_z_P_z_p_bkg[100,0,10000]*P_z_p_spectrum_bkg)")
        # else:
    dws.factory("SUM::Z_m_p_spectrum_all_fit(n_01_Z_m_p*Z_m_p_01_fit,n_0203_Z_m_p*Z_m_p_02_fit,n_04_Z_m_p*Z_m_p_04_fit,n_Z_m_p_bkg[100,0,10000]*Z_m_p_spectrum_bkg)")
    dws.factory("SUM::Z_z_z_spectrum_all_fit(n_09_Z_z_z*Z_z_z_09_fit,n_0710_Z_z_z*Z_z_z_0710_fit,n_040812_Z_z_z*Z_z_z_040812_fit,n_Z_z_z_bkg[100,0,100000]*Z_z_z_spectrum_bkg)")
    dws.factory("SUM::P_z_p_spectrum_all_fit(n_05_P_z_p*P_z_p_05_fit,n_020607_P_z_p*P_z_p_020607_fit,n_0408_P_z_p*P_z_p_0408_fit,n_P_z_p_bkg[100,0,10000]*P_z_p_spectrum_bkg)")
    dws.factory("SUM::M_m_z_spectrum_all_fit(n_03_M_m_z*M_m_z_03_fit,n_04_M_m_z*M_m_z_04_fit,n_M_m_z_bkg[100,0,10000]*M_m_z_spectrum_bkg)")
    dws.factory("SUM::P_z_pst_spectrum_all_fit(n_07_P_z_pst*P_z_pst_07_fit,n_0408_P_z_pst*P_z_pst_0408_fit,n_P_z_pst_bkg[100,0,1000]*P_z_pst_spectrum_bkg)")
    dws.factory("SUM::Zs_sm_p_spectrum_all_fit(n_13_Zs_sm_p*Zs_sm_p_13_fit,n_14_Zs_sm_p*Zs_sm_p_14_fit,n_15_Zs_sm_p*Zs_sm_p_15_fit,n_16_Zs_sm_p*Zs_sm_p_16_fit,n_Zs_sm_p_bkg[100,0,1000]*Zs_sm_p_spectrum_bkg)")

def build_data_ws(run_name, gc_onflag, fit_strat, specs):

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


    dws = ROOT.RooWorkspace(run_name)
    dws.factory(f"B_DTF_M[{bmin},{bmax}]")

    ids_and_shapes = ids_and_shapes_dp
    if fit_strat == "dp_f":
        ###########################
        build_fit(dws, run_name, gc_onflag, fit_strat, specs, ids_and_shapes)
        ###########################
    if fit_strat == "dp_f_nn":
        ###########################
        build_nn_sw_fit(dws, run_name, gc_onflag, fit_strat, specs, ids_and_shapes)
        ###########################

    b_dtf_m = dws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_dtf_m)
    ############################

    all_cats = ROOT.RooCategory("all_cats", "all_cats")

    dspectrum_list = ["Z_m_p_spectrum", "Z_z_z_spectrum", "P_z_p_spectrum","M_m_z_spectrum", "P_z_pst_spectrum", "Zs_sm_p_spectrum"]

    data_set_list = []

    for dspec in dspectrum_list:
        all_cats.defineType(dspec)
        pathname = dspec.split("_spectrum")[0]
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/{pathname}_postdcuts.root")
        tchain = ROOT.TChain("DecayTreeTuple")
        print(tchain.GetEntries())
        for file_name in file_list:
            print(file_name)
            tchain.Add(file_name)
        print(tchain.GetEntries())
        data_set_list.append(ROOT.RooDataSet(f"{dspec}_data", f"{dspec}_data", tchain, data_args))

    all_data_sets = ROOT.RooDataSet(
        "all_data_sets",
        "all_data_sets",
        data_args,
        ROOT.RooFit.Index(all_cats),
        ROOT.RooFit.Import("Z_m_p_spectrum", data_set_list[0]),
        ROOT.RooFit.Import("Z_z_z_spectrum", data_set_list[1]),
        ROOT.RooFit.Import("P_z_p_spectrum", data_set_list[2]),
        ROOT.RooFit.Import("M_m_z_spectrum", data_set_list[3]),
        ROOT.RooFit.Import("P_z_pst_spectrum", data_set_list[4]),
        ROOT.RooFit.Import("Zs_sm_p_spectrum", data_set_list[5])
    )

    getattr(dws, "import")(all_cats)
    getattr(dws, "import")(all_data_sets)

    dws.writeToFile(f"{analysis_path}/data_fit/base_data_files/{run_name}.root")
    print(f"Wrote dws to: {analysis_path}/data_fit/base_data_files/{run_name}.root")
    dws.Print()

def run_data_fit(run_name, gc_onflag, specs):

    lp = 0
    if "Z_m_p" in specs:
        lp = lp + 4
    if "Z_z_z" in specs:
        lp = lp + 6
    if "P_z_p" in specs:
        lp = lp + 6
    if "M_m_z" in specs:
        lp = lp + 2
    if "P_z_pst" in specs:
        lp = lp + 3
    if "Zs_sm_p" in specs:
        lp = lp + 4

    dws_base_file = ROOT.TFile(f"{analysis_path}/data_fit/base_data_files/{run_name}.root")
    dws = dws_base_file.Get(f"{run_name}")

    if gc_onflag == 1:
        glist = ROOT.RooArgSet()
        for i in range(0, lp):
            gvar = dws.pdf(f"{run_name}_nu{i}_g")
            glist.add(gvar)

    all_data_sets = dws.data("all_data_sets")
    all_cats = dws.cat("all_cats")
    for i in all_cats:
        print(i)
    b_dtf_m = dws.var("B_DTF_M")
    all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)

    dws.Print()

    if "Z_m_p" in specs:
        z_model = dws.pdf("Z_m_p_spectrum_all_fit")
        all_fit.addPdf(z_model, "Z_m_p_spectrum")
    if "Z_z_z" in specs:
        zz_model = dws.pdf("Z_z_z_spectrum_all_fit")
        all_fit.addPdf(zz_model, "Z_z_z_spectrum")
    if "P_z_p" in specs:
        p_model = dws.pdf("P_z_p_spectrum_all_fit")
        all_fit.addPdf(p_model, "P_z_p_spectrum")
    if "M_m_z" in specs:
        m_model = dws.pdf("M_m_z_spectrum_all_fit")
        all_fit.addPdf(m_model, "M_m_z_spectrum")
    if "P_z_pst" in specs:
        st_model = dws.pdf("P_z_pst_spectrum_all_fit")
        all_fit.addPdf(st_model, "P_z_pst_spectrum")
    if "Zs_sm_p" in specs:
        s_model = dws.pdf("Zs_sm_p_spectrum_all_fit")
        all_fit.addPdf(s_model, "Zs_sm_p_spectrum")
    if gc_onflag == 1:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(glist), ROOT.RooFit.Save())
    if gc_onflag == 0:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    fitws = ROOT.RooWorkspace(f"fit_{run_name}")
    fitws.Import(all_data_sets)
    fitws.Import(all_cats)
    fitws.Import(all_fit)
    fitws.Import(nall_fit)
    fitws.Print()
    fitws.writeToFile(f"{analysis_path}/data_fit/fit_data_files/fit_{run_name}.root")

    fit_test = fitws.pdf("Z_m_p_spectrum_all_fit")
    data_test = fitws.data("all_data_sets")
    # data_test_2 = ROOT.RooDataSet("Z_m_p_data", "Z_m_p_data", data_test, ROOT.RooArgSet(b_dtf_m), "all_cats==all_cats::Z_m_p_spectrum")
    data_test2 = data_test.reduce(ROOT.RooFit.Cut(f"all_cats==all_cats::Z_m_p_spectrum"))
    nyield_1 = fitws.obj(f"nny_1")
    nyield_2 = fitws.obj(f"nny_2")
    nyield_3 = fitws.obj(f"nny_3")
    n_Z_m_p_bkg = fitws.obj(f"n_Z_m_p_bkg")

    yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, n_Z_m_p_bkg)

    MakeSWeights(f"{analysis_path}/data_fit/sw_files/sw_Z_m_p.root", "SW_tree", data_test2, fit_test, yields)

def plot_data(run_name, fit_strat, specs, data_only_flag = 0):

    print(run_name)
    fws_base_plot_file = ROOT.TFile(f"{analysis_path}/data_fit/fit_data_files/fit_{run_name}.root")
    fws = fws_base_plot_file.Get(f"fit_{run_name}")
    fws.Print()

    b_dtf_m = fws.var("B_DTF_M")
    all_cats = fws.cat("all_cats")

    all_data_sets = fws.data("all_data_sets")
    all_fit = fws.pdf("super_fit_Pdf")

    all_cats_Set = ROOT.RooArgSet(all_cats)
    fitresult = fws.obj("fitresult_super_fit_Pdf_all_data_sets")

    bf_string_dict = {
        "1": "B^{0} -> D- D+ K*0",
        "2": "B^{0} -> D*- D+ K*0",
        "3": "B^{0} -> D- D*+ K*0",
        "4": "B^{0} -> D*- D*+ K*0",
        "5": "B^{+} -> #bar{D0} D+ K*0",
        "6": "B^{+} -> #bar{D*0} D+ K*0",
        "7": "B^{+} -> #bar{D0} D*+ K*0",
        "8": "B^{+} -> #bar{D*0} D*+ K*0",
        "9": "B^{0} -> D0 D0 K*0",
        "10": "B^{0} -> D*0 D0 K*0 + B^{0} -> D0 D*0 K*0",
        "12": "B^{0} -> D*0 D*0 K*0",
        "13": "B_{s}^{0} -> D_{s}- D+ K*0",
        "14": "B_{s}^{0} -> D_{s}*- D+ K*0",
        "15": "B_{s}^{0} -> D_{s}- D*+ K*0",
        "16": "B_{s}^{0} -> D_{s}*- D*+ K*0",
    }

    # file = open("bf.txt","w")
    # for i in range(1,13):
    #     if i != 11:
    #         bf = fws.var(f"bf_{i}")
    #         bf_val = bf.getValV()
    #         bf_err = bf.getPropagatedError(fitresult)
    #         bf_print = ufloat(bf_val, bf_err)
    #         string = bf_string_dict[str(i)]
    #         file.write(f"{string} : {bf_print:.2E} \n")
    # file.close()


    for spec in specs:

        spectrum = f"{spec}_spectrum"
        frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}_spectrum"), ROOT.RooFit.Bins(bbins))
        all_data_sets.plotOn(frame, ROOT.RooFit.Cut(f"all_cats==all_cats::{spec}_spectrum"), ROOT.RooFit.Name("data"))

        if spec == "Z_m_p":
            ylist = ["01","0203","04"]
            title = "D^{-} D^{+} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{+}"

        if spec == "Z_z_z":
            title = "#bar{D^{0}} D^{0} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{0}"
            ylist = ["09", "0710", "040812"]

        if spec == "P_z_p":
            title = "#bar{D^{0}} D^{+} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{+}"
            ylist = ["05", "020607", "0408"]

        if spec == "M_m_z" :
            title = "D^{-} D^{0} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{0}"
            ylist = ["03", "04"]

        if spec == "P_z_pst" :
            title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
            d1 = "D^{-}"
            d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
            ylist = ["07", "0408"]

        if spec == "Zs_sm_p" :
            ylist = ["13","14","15","16"]
            title = "D_{s}^{-} D^{+} K^{*0}"
            d1 = "D_{s}^{-}"
            d2 = "D^{+}"

        if data_only_flag == 0:

            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name("data"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_all_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))

            if spec == "Z_m_p":
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_01_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_02_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                ylist = ["1","2","3"]
                title = "D^{-} D^{+} K^{*0}"
                d1 = "D^{-}"
                d2 = "D^{+}"

            if spec == "Z_z_z":
                title = "#bar{D^{0}} D^{0} K^{*0}"
                d1 = "#bar{D^{0}}"
                d2 = "D^{0}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_09_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0710_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_040812_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                ylist = ["09", "0710", "040812"]

            if spec == "P_z_p":
                title = "#bar{D^{0}} D^{+} K^{*0}"
                d1 = "#bar{D^{0}}"
                d2 = "D^{+}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_05_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_020607_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0408_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                ylist = ["05", "020607", "0408"]

            if spec == "M_m_z" :
                title = "D^{-} D^{0} K^{*0}"
                d1 = "D^{-}"
                d2 = "D^{0}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_03_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_04_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                ylist = ["03", "04"]

            if spec == "P_z_pst" :
                title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
                d1 = "D^{-}"
                d2 = "(D^{*+} #rightarrow D^{0} #pi+)"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_07_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_0408_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                ylist = ["07", "0408"]

            if spec == "Zs_sm_p" :
                ylist = ["13","14","15","16"]
                title = "D_{s}^{-} D^{+} K^{*0}"
                d1 = "D_{s}^{-}"
                d2 = "D^{+}"
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_13_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fit0"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_14_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kTeal), ROOT.RooFit.Name("fit1"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_15_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("fit2"))
                all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_16_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet), ROOT.RooFit.Name("fit3"))


            d_nyield_bkg = fws.obj(f"n_{spec}_bkg")
            d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

        p = ROOT.TCanvas("p1","p1")
        p.cd()

        frame.GetXaxis().SetTitle(f" m({title}) [MeV]")
        frame.Draw()

        if data_only_flag == 0:

            legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
            legend.AddEntry("pdf","Total Fit","lp")
            legend.AddEntry("bkg","Background","l")
            legend.AddEntry("fit0", "Fit 0 Missing Particles","l")
            legend.AddEntry("fit1", "Fit 1 Missing Particles","l")
            if spec != "M_m_z":
                if spec != "P_z_pst":
                    legend.AddEntry("fit2","Fit 2 Missing Particles","l")

            legend.SetTextSize(0.030)
            legend.Draw()

            d_chi2 = frame.chiSquare(f"pdf", f"data")
            if spec != "Z_z_z":
                dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
            if spec == "Z_z_z":
                dtpave = ROOT.TPaveText(0.50, 0.65, 0.70, 0.85, "NB NDC")

            dtpave.SetFillStyle(0)
            dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
            dtpave.AddText(f"Yield BKG: {round(d_nyield_bkg.getValV(),3)} #pm {round(d_nyield_bkg_err,3)}")

            if len(ylist) == 3:
                for y,name in zip(ylist,["Yield Right", "Yield Middle", "Yield Left"]):
                    d_nyield = fws.obj(f"nny_{y}")
                    d_nyield_err = d_nyield.getPropagatedError(fitresult)
                    dtpave.AddText(f"{name}: {round(d_nyield.getValV(),3)} #pm {round(d_nyield_err,3)}")
            if len(ylist) == 2:
                for y,name in zip(ylist,["Yield Right", "Yield Left"]):
                    d_nyield = fws.obj(f"n_{y}_{spec}")
                    d_nyield_err = d_nyield.getPropagatedError(fitresult)
                    dtpave.AddText(f"{name}: {round(d_nyield.getValV(),3)} #pm {round(d_nyield_err,3)}")
            if len(ylist) == 4:
                for y,name in zip(ylist,["Yield Right", "Yield Middle 1","Yield Middle 2", "Yield Left"]):
                    d_nyield = fws.obj(f"n_{y}_{spec}")
                    d_nyield_err = d_nyield.getPropagatedError(fitresult)
                    dtpave.AddText(f"{name}: {round(d_nyield.getValV(),3)} #pm {round(d_nyield_err,3)}")

            dtpave.Draw()
            frame.addObject(dtpave)

        save_png(p, f"fit_tests", f"{spectrum}_{run_name}", rpflag = 0)
