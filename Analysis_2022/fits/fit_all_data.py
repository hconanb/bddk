import sys
import os
import ROOT
import uproot
import pandas as pd
import argparse

basedir = os.getcwd().split('fits')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

def otherget_mc_shape(dws, data_spec, mc_specs, fix_flag, smear_flag, varoi):

    for mc_tuple in mc_specs:

        mc_spec = mc_tuple[0]
        mc_small_spec = mc_spec.split("_1")[0]
        shape_flag = mc_tuple[1]

        mcws_base = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/mc_files/{mc_spec}.root")
        mcws = mcws_base.Get(f"fit_ws")

        vars = mcws.allVars()

        if fix_flag:
            for i in vars:
                if "width_L_02_Z_m_p" not in i.GetName() and "width_R_02_Z_m_p" not in i.GetName():
                    if i.GetName() != varoi and shape_flag in i.GetName():
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")

        mc_pdf = mcws.pdf(f"{mc_small_spec}_fit_{shape_flag}")

        if smear_flag:
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable(f"{mc_small_spec}_fit_{shape_flag}", f"{mc_small_spec}_fit_{shape_flag}_pc"), ROOT.RooFit.RenameVariable("B_DTF_M", f"{varoi}"))
            dws.factory(
                f"FCONV::{mc_small_spec}_fit({varoi}, {mc_small_spec}_fit_{shape_flag}_pc, {data_spec}_gsmear)"
            )
        if not smear_flag:
            #and mc_spec not in ["Z_z_z_0710_fit","Z_z_z_040812_fit","P_z_p_020607_fit","P_z_p_0408_fit"]
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable(f"{mc_small_spec}_fit_{shape_flag}", f"{mc_small_spec}_fit"), ROOT.RooFit.RenameVariable("B_DTF_M", f"{varoi}"))
def get_mc_shapes(dws, base_spec, split_flag, fix_flag, smear_flag, gsmear_strat):
    slist = ids_and_shapes[f"{base_spec}"]
    for tuple in slist:
        mc_spec = tuple[0]
        imp_spec = tuple[0].split("_fit")[0]
        strat = tuple[1]
        mcws_base = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_fit/fit_mc_files/{mc_spec}_{strat}.root")
        mcws = mcws_base.Get(f"fit_ws")
        if fix_flag:
            vars = mcws.allVars()
            for i in vars:
                if not fix_mean_flag:
                    if i.GetName() != "B_DTF_M" and "mean" not in i.GetName():
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")
                else:
                    if i.GetName() != "B_DTF_M":
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")
        if smear_flag:
            mc_pdf = mcws.pdf(f"{mc_spec}")
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable(mc_spec, f"{mc_spec}_pc"))
            if gsmear_strat == "SPLIT" or "P_z_pst" in mc_spec:
                dws.factory(
                    f"FCONV::{mc_spec}(B_DTF_M, {mc_spec}_pc, {base_spec}_gsmear)"
                )
            if gsmear_strat == "ALL" and "P_z_pst" not in mc_spec:
                dws.factory(
                    f"FCONV::{mc_spec}(B_DTF_M, {mc_spec}_pc, a_gsmear)"
                )
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
                    v0=parvals[i], l=snp.sqrt(eigval[j]), e=eigvec[i][j], label=label, j=j
                )
        else:
            expr = "{}".format(parvals[i])
            for j in range(n):
                expr += " + {}_nu{}*{}*{}".format(
                    label, j, np.sqrt(eigval[j]), eigvec[i][j]
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
    norm_base_file = ROOT.TFile(f"data_files/{norm_id}_{year}_{trigger_flag}.root")
    norm_ws = norm_base_file.Get(f"{norm_id}")
    n_base = norm_ws.var(f"n_{norm_id}_signal")
    n_val = ufloat(n_base.getValV(), n_base.getError())
    if "norm8" in norm_id:
        cfactor =  (B_B0_norm8 * dm * d0bar4)
    if "norm7" in norm_id:
        cfactor = (B_B0_norm7 * d0bar * d0bar4)
    return (n_val, cfactor)
def get_eff(id, trigger_flag, year):
    df = pandas.read_excel("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_efficiencies/Harris_numbers.xlsx", index_col=[0,1])
    df.index = df.index.set_levels([df.index.levels[0].astype(str), df.index.levels[1].astype(str)])
    if trigger_flag == "T":
        num = ufloat_fromstr(df.loc[id, year]["Final Bootstraped TOS Events"])
    if trigger_flag == "nTaT":
        num = ufloat_fromstr(df.loc[id, year]["Final Bootstraped TIS Events"])
    den_n = int(df.loc[id, year]["Number Accepted"])
    den = ufloat(den_n, np.sqrt(den_n))
    geneff = ufloat_fromstr(df.loc[id, year]['Generator (%)'])
    eff = (geneff/100)*num/den
    return eff
def get_a_factor(id, norm_id, signal_d_bfs):

    final_bf_factor= 0
    for year in ["2016", "2017", "2018"]:
        for trigger_flag in ["T", "nTaT"]:
            #Norm_Yield is the yield of the noralization mode
            #b_norm_factors is the product of theBranching Fractions for B->D1normD2normK * D1_norm * D2_norm
            Norm_Yield, b_norm_factors = get_norm_param(norm_id, trigger_flag, year)
            print(Norm_Yield, b_norm_factors)
            #efs is the signal mode efficency
    #         efs = get_eff(id, trigger_flag, year)
    #         efn = get_eff(norm_id, trigger_flag, year)
    #
    #         #This is the final number for a given year, trigger condition
    #         #signal_bfs is the produut of the signal branching fractions for D1_sig*D2_sig*Kst0
    #         signal_bfs = signal_d_bfs*Kst_factor
    #
    #         final_base_factor = Norm_Yield*(signal_bfs/b_norm_factors)*(efs/efn)
    #
    #
    # #        #this is the sum of all final_base_factors
    #         final_bf_factor = final_bf_factor + final_base_factor
    # #


    return final_bf_factor
def get_params(dws):

    params = []

    a_01_z = get_a_factor("1", "norm8", dm*dm)
    a_02_z = get_a_factor("2a,3a", "norm8", dstmdmpi0 * dm * dm)
    a_03_z = get_a_factor("2a,3a", "norm8", dstmdmpi0 * dm * dm)
    a_04_z = get_a_factor("4a", "norm8", dstmdmpi0 * dstmdmpi0 * dm * dm)

    params = params + [a_01_z, a_02_z, a_03_z, a_04_z]

    # a_04_zz = get_a_factor("4c", "norm7", dstmd0pim * dstmd0pim * d0bar * d0bar)
    #
    # a_07_zz = get_a_factor("7b", "norm7", dstmd0pim * d0bar * d0bar)
    # a_08_zz = get_a_factor("8b", "norm7", dstmd0pim * d0bar * d0bar)
    #
    # a_09_zz = get_a_factor("9", "norm7", d0bar * d0bar)
    # a_10_zz = get_a_factor("10", "norm7", d0bar * d0bar)
    # a_12_zz = get_a_factor("12", "norm7", d0bar * d0bar)
    #
    # params = params + [a_04_zz, a_07_zz, a_08_zz, a_09_zz, a_10_zz, a_12_zz]
    #
    # a_02_p = get_a_factor("2b,3b", "norm7", dstmd0pim * dm * d0bar)
    # a_04_p = get_a_factor("4b", "norm7", dstmd0pim * dstmdmpi0 * d0bar * dm)
    #
    # a_05_p = get_a_factor("5", "norm7", dm * d0bar)
    # a_06_p = get_a_factor("6", "norm7", dm * d0bar)
    #
    # a_07_p = get_a_factor("7b", "norm7", d0bar * dstmdmpi0 * dm)
    # a_08_p = get_a_factor("8b", "norm7", d0bar * dstmdmpi0 * dm)
    #
    # params = params + [a_02_p, a_04_p, a_05_p, a_06_p, a_07_p, a_08_p]
    #
    # a_03_m = get_a_factor("2b,3b", "norm7", dstmd0pim * d0bar * dm)
    # a_04_m = get_a_factor("4b", "norm7", dstmdmpi0 * dstmd0pim * dm * d0bar)
    #
    # params = params + [a_03_m, a_04_m]
    #
    # a_04_st = get_a_factor("4d", "norm7", dstmd0pim * dstmd0pim * d0bar * d0bar)
    # a_07_st = get_a_factor("7c", "norm7", dstmd0pim * d0bar * d0bar)
    # a_08_st = get_a_factor("8c", "norm7", dstmd0pim * d0bar * d0bar)
    #
    # params = params + [a_04_st, a_07_st, a_08_st]

    return params
def fit_data(fix_flag, smear_flag, gc_onflag):

    dws = ROOT.RooWorkspace("fit_ws")
    bmin = 4800
    bmax = 5600

    dws.factory(f"B_DTF_M[{bmin},{bmax}]")

    b_dtf_m = dws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_dtf_m)

    params = get_params(dws)
    c = covariance_matrix(params)

    test_a, test_b = makeeigenpars(run_name, [p.n for p in params], c, debug=False)

    for a in test_a:
        print(a)
        if gc_onflag:
            dws.Import(a)
        if not gc_onflag:
            new_a = a.evaluate()
            new_a_name = a.GetName()
            dws.factory(f"{new_a_name}[{new_a}]")
            tv = dws.var(f"{new_a_name}")
            tv.Print()
    if gc_onflag:
        for i in range(0, len(params)):
            dws.factory(f"Gaussian::{run_name}_nu{i}_g({run_name}_nu{i}, 0, 1)")
    for i in range(1,13):
        dws.factory(f"bf_{i}[0.001,0,0.01]")

    dws.factory(f"expr::Z_m_p_01_yield('bf_1*{run_name}_par0', bf_1, {run_name}_par0)")
    dws.factory(f"expr::Z_m_p_02_yield('bf_2*{run_name}_par1', bf_2, {run_name}_par1)")
    dws.factory(f"expr::Z_m_p_03_yield('bf_3*{run_name}_par2', bf_3, {run_name}_par2)")
    dws.factory(f"expr::Z_m_p_04_yield('bf_4*{run_name}_par3', bf_4, {run_name}_par3)")

    dws.factory(f"expr::Z_z_z_04_yield('bf_4*{run_name}_par4', bf_4, {run_name}_par4)")
    dws.factory(f"expr::Z_z_z_07_yield('bf_7*{run_name}_par5', bf_7, {run_name}_par5)")
    dws.factory(f"expr::Z_z_z_08_yield('bf_8*{run_name}_par6', bf_8, {run_name}_par6)")
    dws.factory(f"expr::Z_z_z_09_yield('bf_9*{run_name}_par7', bf_9, {run_name}_par7)")
    dws.factory(f"expr::Z_z_z_10_yield('bf_10*{run_name}_par8', bf_10, {run_name}_par8)")
    dws.factory(f"expr::Z_z_z_12_yield('bf_12*{run_name}_par9', bf_12, {run_name}_par9)")

    dws.factory(f"expr::P_z_p_02_yield('bf_2*{run_name}_par10', bf_2, {run_name}_par10)")
    dws.factory(f"expr::P_z_p_04_yield('bf_4*{run_name}_par11', bf_4, {run_name}_par11)")
    dws.factory(f"expr::P_z_p_05_yield('bf_5*{run_name}_par12', bf_5, {run_name}_par12)")
    dws.factory(f"expr::P_z_p_06_yield('bf_6*{run_name}_par13', bf_6, {run_name}_par13)")
    dws.factory(f"expr::P_z_p_07_yield('bf_7*{run_name}_par14', bf_7, {run_name}_par14)")
    dws.factory(f"expr::P_z_p_08_yield('bf_8*{run_name}_par15', bf_8, {run_name}_par15)")

    dws.factory(f"expr::M_m_z_03_yield('bf_3*{run_name}_par16', bf_3, {run_name}_par16)")
    dws.factory(f"expr::M_m_z_04_yield('bf_4*{run_name}_par17', bf_4, {run_name}_par17)")

    dws.factory(f"expr::P_z_pst_04_yield('bf_4*{run_name}_par18', bf_4, {run_name}_par18)")
    dws.factory(f"expr::P_z_pst_07_yield('bf_7*{run_name}_par19', bf_7, {run_name}_par19)")
    dws.factory(f"expr::P_z_pst_08_yield('bf_8*{run_name}_par20', bf_8, {run_name}_par20)")


    # dws.factory(f"width_P_z_pst_gsmear[1, 0.0001, 2]")
    # gsmear_strat = "ALL"
    # if gsmear_strat == "SPLIT":
    #     for spec in specs:
    #         if spec != "P_z_pst":
    #             dws.factory(f"width_{spec}_gsmear[1, 0.001, 5]")
    #         dws.factory(f"mean_{spec}_gsmear[0]")
    #         dws.factory(f"Gaussian::{spec}_gsmear(B_DTF_M, mean_{spec}_gsmear, width_{spec}_gsmear)")

    # if gsmear_strat == "ALL":
    #     dws.factory(f"width_gsmear[1, 0.001, 5]")
    #     dws.factory(f"mean_gsmear[0]")
    #     dws.factory(f"Gaussian::a_gsmear(B_DTF_M, mean_gsmear, width_gsmear)")
    #     dws.factory(f"Gaussian::P_z_pst_gsmear(B_DTF_M, mean_gsmear, width_P_z_pst_gsmear)")
    #
    #
    # for spec in specs:
    #     get_shapes_bkg(spec, "Exponential", dws)
    #     # get_mc_shape_corr(dws, spec, split_flag = True, fix_flag = True, smear_flag = True, gsmear_strat = "SPLIT")
    #     get_mc_shape_corr(dws, spec, split_flag = True, fix_flag = True, smear_flag = True, gsmear_strat = "ALL")
    #
    # dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_yield*Z_m_p_01_fit, Z_m_p_02_yield*Z_m_p_02_fit, Z_m_p_03_yield*Z_m_p_03_fit, Z_m_p_04_yield*Z_m_p_04_fit, Z_m_p_bkg_yield[100,0,100000]* Z_m_p_spectrum_bkg)")
    # dws.factory("SUM::Z_z_z_spectrum_all_fit(Z_z_z_09_yield*Z_z_z_09_fit, Z_z_z_07_yield*Z_z_z_07_fit, Z_z_z_10_yield*Z_z_z_10_fit, Z_z_z_04_yield*Z_z_z_04_fit, Z_z_z_08_yield*Z_z_z_08_fit, Z_z_z_12_yield*Z_z_z_12_fit, Z_z_z_bkg_yield[100,0,100000]*Z_z_z_spectrum_bkg)")
    # dws.factory("SUM::P_z_p_spectrum_all_fit(P_z_p_05_yield*P_z_p_05_fit,P_z_p_02_yield*P_z_p_02_fit,P_z_p_06_yield*P_z_p_06_fit,P_z_p_07_yield*P_z_p_07_fit,P_z_p_04_yield*P_z_p_04_fit,P_z_p_08_yield*P_z_p_08_fit,P_z_p_bkg_yield[100,0,10000]*P_z_p_spectrum_bkg)")
    # dws.factory("SUM::M_m_z_spectrum_all_fit(M_m_z_03_yield*M_m_z_03_fit,M_m_z_04_yield*M_m_z_04_fit,M_m_z_bkg_yield[100,0,10000]*M_m_z_spectrum_bkg)")
    # dws.factory("SUM::P_z_pst_spectrum_all_fit(P_z_pst_07_yield*P_z_pst_07_fit,P_z_pst_04_yield*P_z_pst_04_fit,P_z_pst_08_yield*P_z_pst_08_fit,P_z_pst_bkg_yield[100,0,1000]*P_z_pst_spectrum_bkg)")
    # # dws.factory("SUM::Zs_sm_p_spectrum_all_fit(Zs_sm_p_13_yield*Zs_sm_p_13_fit, Zs_sm_p_14_yield*Zs_sm_p_14_fit, Zs_sm_p_15_yield*Zs_sm_p_15_fit, Zs_sm_p_16_yield*Zs_sm_p_16_fit, Zs_sm_p_bkg_yield[100,0,10000]*Zs_sm_p_spectrum_bkg)")
    #
    # all_cats = ROOT.RooCategory("all_cats", "all_cats")
    # dspectrum_list = ["Z_m_p_spectrum", "Z_z_z_spectrum", "P_z_p_spectrum","M_m_z_spectrum", "P_z_pst_spectrum"]
    #
    # data_set_list = []
    #
        # file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{id_to_spec_dict[data_spec]}.root"
        # DecayTree_List = ["DecayTreeTuple_SIG"]
        # tchain = grab_files_and_chain(file_path, DecayTree_List)
    # for dspec in dspectrum_list:
    #     all_cats.defineType(dspec)
    #     pathname = dspec.split("_spectrum")[0]
    #     file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{pathname}.root")
    #     tchain = ROOT.TChain("DecayTreeTuple_SIG")
    #     for file_name in file_list:
    #         tchain.Add(file_name)
    #     data_set_list.append(ROOT.RooDataSet(f"{dspec}_data", f"{dspec}_data", tchain, data_args))
    #
    # all_data_sets = ROOT.RooDataSet(
    #     "all_data_sets",
    #     "all_data_sets",
    #     data_args,
    #     ROOT.RooFit.Index(all_cats),
    #     ROOT.RooFit.Import("Z_m_p_spectrum", data_set_list[0]),
    #     ROOT.RooFit.Import("Z_z_z_spectrum", data_set_list[1]),
    #     ROOT.RooFit.Import("P_z_p_spectrum", data_set_list[2]),
    #     ROOT.RooFit.Import("M_m_z_spectrum", data_set_list[3]),
    #     ROOT.RooFit.Import("P_z_pst_spectrum", data_set_list[4]),
    # )
    #
    # getattr(dws, "import")(all_cats)
    # getattr(dws, "import")(all_data_sets)
    #
    # if gc_onflag == 1:
    #     lp = 0
    #     if "Z_m_p" in specs:
    #         lp = lp + 4
    #     if "Z_z_z" in specs:
    #         lp = lp + 6
    #     if "P_z_p" in specs:
    #         lp = lp + 6
    #     if "M_m_z" in specs:
    #         lp = lp + 2
    #     if "P_z_pst" in specs:
    #         lp = lp + 3
    #     glist = ROOT.RooArgSet()
    #     for i in range(0, lp):
    #         gvar = dws.pdf(f"{run_name}_nu{i}_g")
    #         glist.add(gvar)
    #
    # all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)
    # z_model = dws.pdf("Z_m_p_spectrum_all_fit")
    # all_fit.addPdf(z_model, "Z_m_p_spectrum")
    # zz_model = dws.pdf("Z_z_z_spectrum_all_fit")
    # all_fit.addPdf(zz_model, "Z_z_z_spectrum")
    # p_model = dws.pdf("P_z_p_spectrum_all_fit")
    # all_fit.addPdf(p_model, "P_z_p_spectrum")
    # m_model = dws.pdf("M_m_z_spectrum_all_fit")
    # all_fit.addPdf(m_model, "M_m_z_spectrum")
    # st_model = dws.pdf("P_z_pst_spectrum_all_fit")
    # all_fit.addPdf(st_model, "P_z_pst_spectrum")
    # # s_model = dws.pdf("Zs_sm_p_spectrum_all_fit")
    # # all_fit.addPdf(s_model, "Zs_sm_p_spectrum")
    #
    # if gc_onflag == 1:
    #     nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(glist), ROOT.RooFit.Save())
    # if gc_onflag == 0:
    #     nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())
    #
    # # dws.Print()
    #
    # file = open("bf.txt","w")
    # for i in range(1,13):
    #     if i != 11:
    #         bf = dws.var(f"bf_{i}")
    #         bf_val = bf.getValV()
    #         bf_err = bf.getPropagatedError(nall_fit)
    #         bf_print = ufloat(bf_val, bf_err)
    #         string = bf_string_dict[str(i)]
    #         file.write(f"{string} : {bf_val:.2E} {bf_err:.2E} \n")
    # file.close()
    #
    # dws.Import(nall_fit)
    # dws.Import(all_fit)
    #
    # output_base_data = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/corr_{run_name}.root"
    # dws.writeToFile(output_base_data)
    # print(f"Wrote dws to: {output_base_data}")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Fit a PDF to all data")
    parser.add_argument('--FIT', action='store_true')
    parser.add_argument('--PLOT', action='store_true')

    parser.add_argument('--FIX', action='store_true')
    parser.add_argument('--SMEAR', action='store_true')
    parser.add_argument('--GCON', action='store_true')

    args = parser.parse_args()

    fit_flag = args.FIT
    plot_flag = args.PLOT
    fix_flag = args.FIX
    smear_flag = args.SMEAR
    gc_onflag = args.GCON



    if fit_flag:
        fit_data(fix_flag, smear_flag, gc_onflag)
    if plot_flag:
        plot_data()
