import ROOT as ROOT
from Analysis_2021.essentials import get_shapes_bkg, get_free_shapes, save_png, save_pdf, MakeSWeights
import glob as glob
from Analysis_2021 import rootutils
from uncertainties import ufloat
from uncertainties import covariance_matrix, correlation_matrix
from uncertainties import unumpy, ufloat_fromstr
import math as math
import pandas
import numpy as np
RDF = ROOT.ROOT.RDataFrame

#00876c
#84b76e
#f4e07f
#ef9250
#d43d51

col_hex_list = ["#00876c", "#f4e07f", "#ef9250", "#d43d51", "#84b76e"]
col_base_list = []
for c in col_hex_list:
    col_base_list.append(ROOT.TColor.GetColor(c))

sig_col = col_base_list[0]
bkg_col = col_base_list[1]
p0_col = col_base_list[2]
p1_col = col_base_list[3]
p2_col = col_base_list[4]

bf_string_dict = {
    "1": "B^{0} -> D^{-} D^{+} K^{*0}",
    "2": "B^{0} -> D^{*-} D^{+} K^{*0}",
    "3": "B^{0} -> D^{-} D^{*+} K^{*0}",
    "4": "B^{0} -> D^{*-} D^{*+} K^{*0}",
    "5": "B^{+} -> #bar{D}^{0} D^{+} K^{*0}",
    "6": "B^{+} -> #bar{D}^{*0} D^{+} K^{*0}",
    "7": "B^{+} -> #bar{D}^{0} D^{*+} K^{*0}",
    "8": "B^{+} -> #bar{D}^{*0} D^{*+} K^{*0}",
    "9": "B^{0} -> #bar{D}^{0} D^{0} K^{*0}",
    "10": "B^{0} -> #bar{D}^{*0} D^{0} K^{*0} or #bar{D}^{0} D^{*0} K^{*0}",
    "12": "B^{0} -> #bar{D}^{*0} D^{*0} K^{*0}",
    "13": "B_{s}^{0} -> D_{s}- D+ K^{*0}",
    "14": "B_{s}^{0} -> D_{s}*- D+ K^{*0}",
    "15": "B_{s}^{0} -> D_{s}- D*+ K^{*0}",
    "16": "B_{s}^{0} -> D_{s}*- D*+ K^{*0}",
}

fix_mean_flag = True

Kst_factor = ufloat(2 / 3, 0)

dm = ufloat(0.0938, 0.0016)
dsm = ufloat(0.0539, 0.0015)

d0bar = ufloat(0.0395, 0.00031)
d0bar4 = ufloat(0.0823, 0.0014)

dstmdmpi0 = ufloat(0.323, 0.005)
dstmd0pim = ufloat(0.677, 0.005)

B_B0_norm7 = ufloat(1.31e-03, 7.00e-05)
B_B0_norm8 = ufloat(1.07e-03, 1.10e-04)

analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"

ids_and_shapes = {  "Z_m_p_no_split_B_M" : [("Z_m_p_01_fit", "DG"), ("Z_m_p_0203_fit" , "BG"), ("Z_m_p_04_fit" , "BGEP")],
                    "Z_m_p_no_split_B_DTF_M" : [("Z_m_p_01_fit", "DG"), ("Z_m_p_0203_fit" , "BG"), ("Z_m_p_04_fit" , "BGEP")],
                    "Z_m_p" : [("Z_m_p_01_fit", "DG"), ("Z_m_p_02_fit" , "BGEP"), ("Z_m_p_03_fit" , "BGEP"), ("Z_m_p_04_fit" , "BGEP")],
                    ####
                    "Z_z_z_no_split_B_M" : [("Z_z_z_09_fit", "DG"), ("Z_z_z_0710_fit" , "GAddBGEP_fr"), ("Z_z_z_040812_fit" , "GAddBGEP_fr")],
                    "Z_z_z_no_split_B_DTF_M" : [("Z_z_z_09_fit", "DG"), ("Z_z_z_0710_fit" , "GAddBGEP_fr"), ("Z_z_z_040812_fit" , "GAddBGEP_fr")],
                    "Z_z_z" : [("Z_z_z_09_fit", "DG"), ("Z_z_z_07_fit" ,  "BGEP"), ("Z_z_z_10_fit" , "GEPAddBGEP_fr"), ("Z_z_z_04_fit" , "BGEP"), ("Z_z_z_08_fit" , "GAddBGEP_fr"), ("Z_z_z_12_fit" , "GAddBGEP_fr")],
                    ####
                    "P_z_p_no_split_B_M" : [("P_z_p_020607_fit", "GEPAddBGEP_fr"), ("P_z_p_0408_fit", "GAddBGEP_fr"), ("P_z_p_05_fit", "DG")],
                    "P_z_p_no_split_B_DTF_M" : [("P_z_p_020607_fit", "GEPAddBGEP_fr"), ("P_z_p_0408_fit", "GAddBGEP_fr"), ("P_z_p_05_fit", "DG")],
                    "P_z_p" : [("P_z_p_02_fit", "BGEP"), ("P_z_p_04_fit", "BGEP"), ("P_z_p_05_fit", "DG"), ("P_z_p_06_fit", "GEPAddBGEP_fr"), ("P_z_p_07_fit", "BGEP"), ("P_z_p_08_fit", "GAddBGEP_fr")],
                    ###
                    "M_m_z_no_split_B_M" : [("M_m_z_03_fit", "BGEP"), ("M_m_z_04_fit", "BGEP")],
                    "M_m_z_no_split_B_DTF_M" : [("M_m_z_03_fit", "BGEP"), ("M_m_z_04_fit", "BGEP")],
                    "M_m_z" : [("M_m_z_03_fit", "BGEP"), ("M_m_z_04_fit", "BGEP")],
                    ###
                    "P_z_pst_no_split_B_M" : [("P_z_pst_07_fit", "DG"), ("P_z_pst_0408_fit",  "GAddBGEP_fr")],
                    "P_z_pst_no_split_B_DTF_M" : [("P_z_pst_07_fit", "DG"), ("P_z_pst_0408_fit", "GAddBGEP_fr")],
                    "P_z_pst" : [("P_z_pst_07_fit", "DG"), ("P_z_pst_04_fit", "BGEP"), ("P_z_pst_08_fit", "BGEP")],
                  }

def get_mc_shape_nn(dws, base_spec, split_flag, fix_flag, smear_flag, varoi):
    slist = ids_and_shapes[f"{base_spec}_no_split_{varoi}"]
    for tuple in slist:

        mc_spec = tuple[0]
        imp_spec = tuple[0].split("_fit")[0]
        strat = tuple[1]

        mcws_base = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_fit/fit_mc_files/{mc_spec}_{strat}.root")
        mcws = mcws_base.Get(f"fit_ws")
        mc_pdf = mcws.pdf(f"{mc_spec}")

        if mc_spec in ["Z_z_z_0710_fit","Z_z_z_040812_fit","P_z_p_020607_fit","P_z_p_0408_fit","P_z_pst_0408_fit"]:
            vars = mcws.allVars()
            for i in vars:
                if i.GetName() != "B_DTF_M" and i.GetName() != "B_M" and "frac" not in i.GetName():
                    temp = mcws.var(i.GetName())
                    temp.setConstant(True)
                    print(f"{temp} is constant")
            fix_flag = False
            smear_flag = True
        if mc_spec in ["Z_m_p_01_fit", "Z_m_p_04_fit","Z_z_z_09_fit","P_z_p_05_fit","M_m_z_03_fit","M_m_z_04_fit","P_z_pst_07_fit"]:

            fix_flag = True
            smear_flag = True

        if mc_spec in ["Z_m_p_0203_fit"]:

            fix_flag = False
            smear_flag = False

        if fix_flag:
            vars = mcws.allVars()
            for i in vars:
                if not fix_mean_flag:
                    if i.GetName() != "B_DTF_M" and i.GetName() != "B_M" and "mean" not in i.GetName():
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")
                else:
                    if i.GetName() != "B_DTF_M" and i.GetName() != "B_M":
                        temp = mcws.var(i.GetName())
                        temp.setConstant(True)
                        print(f"{temp} is constant")
        # print(smear_flag)
        if smear_flag:
            mc_pdf = mcws.pdf(f"{mc_spec}")
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable(mc_spec, f"{mc_spec}_pc"), ROOT.RooFit.RenameVariable("B_DTF_M", f"{varoi}"))
            dws.factory(
                f"FCONV::{mc_spec}({varoi}, {mc_spec}_pc, {base_spec}_gsmear)"
            )

        if not smear_flag:
            #and mc_spec not in ["Z_z_z_0710_fit","Z_z_z_040812_fit","P_z_p_020607_fit","P_z_p_0408_fit"]
            dws.Import(mc_pdf, ROOT.RooFit.RenameVariable("B_DTF_M", f"{varoi}"))

def get_mc_shape_corr(dws, base_spec, split_flag, fix_flag, smear_flag, gsmear_strat):
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
    norm_base_file = ROOT.TFile(f"{analysis_path}/normalization_fit/fit_files/{norm_id}_{year}_{trigger_flag}.root")
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
            efn = get_eff(norm_id, trigger_flag, year)
            # print(efs, efn)
            #This is the final number for a given year, trigger condition
            #signal_bfs is the produut of the signal branching fractions for D1_sig*D2_sig*Kst0
            signal_bfs = signal_d_bfs*Kst_factor

            final_base_factor = Norm_Yield*(signal_bfs/b_norm_factors)*(efs/efn)
            # print(final_base_factor)
    #        # dict_tre[f"{year}_{trigger_flag}"] = f"{signal_bfs}"
    #         dict_tre[f"{year}_{trigger_flag}"] = f"{final_base_factor:.3e}"
    #
    #        #this is the sum of all final_base_factors
            final_bf_factor = final_bf_factor + final_base_factor
    #
    # df_bfs_list.append(signal_bfs/b_norm_factors)
    # df_bfs_list.append(signal_bfs)
    # df_n_rows_list.append(dict)
    # df_tre_rows_list.append(dict_tre)
    # df_index_list.append(id)
    # df_tre_index_list.append(id)
    #
    # #  add norm modes
    # if id == "16_s":
    #     for norm_id in ["norm7_norm7", "norm8_norm8"]:
    #         dict = {}
    #         for year in ["2016", "2017", "2018"]:
    #             for trigger_flag in ["T", "nTaT"]:
    #                 efn = get_eff(norm_id, trigger_flag, year)
    #                 dict[f"{year}_{trigger_flag}"] = f"{efn:.3e}"
    #         df_n_rows_list.append(dict)
    #         df_index_list.append(norm_id)

    return final_bf_factor
def get_params(dws):
    #Build n params
    params = []

    a_01_z = get_a_factor("1", "norm8", dm*dm)
    a_02_z = get_a_factor("2a,3a", "norm8", dstmdmpi0 * dm * dm)
    a_03_z = get_a_factor("2a,3a", "norm8", dstmdmpi0 * dm * dm)
    a_04_z = get_a_factor("4a", "norm8", dstmdmpi0 * dstmdmpi0 * dm * dm)

    params = params + [a_01_z, a_02_z, a_03_z, a_04_z]

    a_04_zz = get_a_factor("4c", "norm7", dstmd0pim * dstmd0pim * d0bar * d0bar)

    a_07_zz = get_a_factor("7b", "norm7", dstmd0pim * d0bar * d0bar)
    a_08_zz = get_a_factor("8b", "norm7", dstmd0pim * d0bar * d0bar)

    a_09_zz = get_a_factor("9", "norm7", d0bar * d0bar)
    a_10_zz = get_a_factor("10", "norm7", d0bar * d0bar)
    a_12_zz = get_a_factor("12", "norm7", d0bar * d0bar)

    params = params + [a_04_zz, a_07_zz, a_08_zz, a_09_zz, a_10_zz, a_12_zz]

    a_02_p = get_a_factor("2b,3b", "norm7", dstmd0pim * dm * d0bar)
    a_04_p = get_a_factor("4b", "norm7", dstmd0pim * dstmdmpi0 * d0bar * dm)

    a_05_p = get_a_factor("5", "norm7", dm * d0bar)
    a_06_p = get_a_factor("6", "norm7", dm * d0bar)

    a_07_p = get_a_factor("7b", "norm7", d0bar * dstmdmpi0 * dm)
    a_08_p = get_a_factor("8b", "norm7", d0bar * dstmdmpi0 * dm)

    params = params + [a_02_p, a_04_p, a_05_p, a_06_p, a_07_p, a_08_p]

    a_03_m = get_a_factor("2b,3b", "norm7", dstmd0pim * d0bar * dm)
    a_04_m = get_a_factor("4b", "norm7", dstmdmpi0 * dstmd0pim * dm * d0bar)

    params = params + [a_03_m, a_04_m]

    a_04_st = get_a_factor("4d", "norm7", dstmd0pim * dstmd0pim * d0bar * d0bar)
    a_07_st = get_a_factor("7c", "norm7", dstmd0pim * d0bar * d0bar)
    a_08_st = get_a_factor("8c", "norm7", dstmd0pim * d0bar * d0bar)

    params = params + [a_04_st, a_07_st, a_08_st]

    # a_13_s = get_a_factor("13", "norm8", dsm * dm)
    # a_14_s = get_a_factor("14", "norm8", dsm * dm)
    # a_15_s = get_a_factor("15", "norm8", dstmdmpi0 * dsm * dm)
    # a_16_s = get_a_factor("16", "norm8", dstmdmpi0 * dsm * dm)
    #
    # params = params + [a_13_s, a_14_s, a_15_s, a_16_s]

    return params

def build_nn_sw_fit(run_name, spec, split_flag, fix_flag, smear_flag, varoi = "B_DTF_M"):

    print(varoi)
    dws = ROOT.RooWorkspace(run_name)
    bmin = 4800
    bmax = 5600
    dws.factory(f"{varoi}[{bmin},{bmax}]")

    dws.factory(f"width_{spec}_gsmear[10, 0.01, 30]")
    dws.factory(f"mean_{spec}_gsmear[0]")
    dws.factory(f"Gaussian::{spec}_gsmear({varoi}, mean_{spec}_gsmear, width_{spec}_gsmear)")

    fix_flag = False
    get_mc_shape_nn(dws, spec, split_flag, fix_flag, smear_flag, varoi)
    get_shapes_bkg(spec, "Exponential", dws, varoi)

    if "Z_m_p" in spec:
        dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_yield[500,0,10000]*Z_m_p_01_fit, Z_m_p_0203_yield[500,0,10000]*Z_m_p_0203_fit, Z_m_p_04_yield[500,0,10000]*Z_m_p_04_fit, Z_m_p_bkg_yield[100,0,100000]* Z_m_p_spectrum_bkg)")
    if "Z_z_z" in spec:
        dws.factory("SUM::Z_z_z_spectrum_all_fit(Z_z_z_09_yield[500,0,10000]*Z_z_z_09_fit,Z_z_z_0710_yield[500,0,10000]*Z_z_z_0710_fit,Z_z_z_040812_yield[500,0,10000]*Z_z_z_040812_fit,Z_z_z_bkg_yield[100,0,100000]*Z_z_z_spectrum_bkg)")
    if "P_z_p" in spec and "P_z_pst" not in spec:
        dws.factory("SUM::P_z_p_spectrum_all_fit(P_z_p_05_yield[500,0,10000]*P_z_p_05_fit,P_z_p_020607_yield[500,0,10000]*P_z_p_020607_fit,P_z_p_0408_yield[500,0,10000]*P_z_p_0408_fit,P_z_p_bkg_yield[100,0,10000]*P_z_p_spectrum_bkg)")
    if "M_m_z" in spec:
        dws.factory("SUM::M_m_z_spectrum_all_fit(M_m_z_03_yield[500,0,10000]*M_m_z_03_fit,M_m_z_04_yield[500,0,10000]*M_m_z_04_fit,M_m_z_bkg_yield[100,0,10000]*M_m_z_spectrum_bkg)")
    if "P_z_pst" in spec:
        dws.factory("SUM::P_z_pst_spectrum_all_fit(P_z_pst_07_yield[500,0,10000]*P_z_pst_07_fit,P_z_pst_0408_yield[500,0,10000]*P_z_pst_0408_fit,P_z_pst_bkg_yield[100,0,1000]*P_z_pst_spectrum_bkg)")

    b_dtf_m = dws.var(f"{varoi}")
    data_args = ROOT.RooArgSet(b_dtf_m)
    #############################

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{spec}.root")
    tchain = ROOT.TChain("DecayTreeTuple_SIG")
    for file_name in file_list:
        tchain.Add(file_name)

    data = ROOT.RooDataSet(f"{spec}_final_data", f"{spec}_final_data", tchain, data_args)
    model = dws.pdf(f"{spec}_spectrum_all_fit")
    fit = model.fitTo(data, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    dws.Import(data)
    dws.Import(fit)

    output_base_data = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/{spec}_{run_name}_{varoi}.root"
    dws.writeToFile(output_base_data)
    print(f"Wrote dws to: {output_base_data}")
    dws.Print()

    if varoi == "B_DTF_M":
        if spec == "Z_m_p":
            nyield_1 = dws.var(f"Z_m_p_01_yield")
            nyield_2 = dws.var(f"Z_m_p_0203_yield")
            nyield_3 = dws.var(f"Z_m_p_04_yield")
            nyield_bkg = dws.var(f"Z_m_p_bkg_yield")
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, nyield_bkg)
        if spec == "Z_z_z":
            nyield_1 = dws.var(f"Z_z_z_09_yield")
            nyield_2 = dws.var(f"Z_z_z_0710_yield")
            nyield_3 = dws.var(f"Z_z_z_040812_yield")
            nyield_bkg = dws.var(f"Z_z_z_bkg_yield")
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, nyield_bkg)
        if spec == "P_z_p":
            nyield_1 = dws.var(f"P_z_p_05_yield")
            nyield_2 = dws.var(f"P_z_p_020607_yield")
            nyield_3 = dws.var(f"P_z_p_0408_yield")
            nyield_bkg = dws.var(f"P_z_p_bkg_yield")
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, nyield_bkg)
        if spec == "M_m_z":
            nyield_1 = dws.var(f"M_m_z_03_yield")
            nyield_2 = dws.var(f"M_m_z_04_yield")
            nyield_bkg = dws.var(f"M_m_z_bkg_yield")
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_bkg)
        if spec == "P_z_pst":
            nyield_1 = dws.var(f"P_z_pst_0408_yield")
            nyield_2 = dws.var(f"P_z_pst_07_yield")
            nyield_bkg = dws.var(f"P_z_pst_bkg_yield")
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_bkg)
        if spec == "Zs_sm_p":
            nyield_1 = dws.var(f"Zs_sm_p_13_yield")
            nyield_2 = dws.var(f"Zs_sm_p_14_yield")
            nyield_3 = dws.var(f"Zs_sm_p_15_yield")
            nyield_4 = dws.var(f"Zs_sm_p_16_yield")
            nyield_bkg = dws.var(f"Zs_sm_p_bkg_yield")
            yields = ROOT.RooArgSet(nyield_1, nyield_2, nyield_3, nyield_4, nyield_bkg)

        MakeSWeights(f"{analysis_path}/data_fit/sw_files/sw_{spec}_{run_name}.root", "SW_tree", data, model, yields)
def build_comp_ws_fit(run_name, specs, split_flag, fix_flag, smear_flag, gc_onflag):

    dws = ROOT.RooWorkspace(run_name)
    bmin = 4800
    bmax = 5600
    dws.factory(f"B_DTF_M[{bmin},{bmax}]")
    b_dtf_m = dws.var("B_DTF_M")
    data_args = ROOT.RooArgSet(b_dtf_m)

    params = get_params(dws)
    c = covariance_matrix(params)

    test_a, test_b = makeeigenpars(run_name, [p.n for p in params], c, debug=False)
    print(test_a)
    print(test_b)
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

    # dws.factory(f"LinearVar::Z_m_p_01_yield(bf_1, {run_name}_par0, 0)")
    # dws.factory(f"LinearVar::Z_m_p_02_yield(bf_2, {run_name}_par1, 0)")
    # dws.factory(f"LinearVar::Z_m_p_03_yield(bf_3, {run_name}_par2, 0)")
    # dws.factory(f"LinearVar::Z_m_p_04_yield(bf_4, {run_name}_par3, 0)")

    # dws.factory(f"expr::n_0203_Z_m_p( 'n_02_Z_m_p + n_03_Z_m_p', n_02_Z_m_p, n_03_Z_m_p)")

    dws.factory(f"expr::Z_z_z_04_yield('bf_4*{run_name}_par4', bf_4, {run_name}_par4)")
    dws.factory(f"expr::Z_z_z_07_yield('bf_7*{run_name}_par5', bf_7, {run_name}_par5)")
    dws.factory(f"expr::Z_z_z_08_yield('bf_8*{run_name}_par6', bf_8, {run_name}_par6)")
    dws.factory(f"expr::Z_z_z_09_yield('bf_9*{run_name}_par7', bf_9, {run_name}_par7)")
    dws.factory(f"expr::Z_z_z_10_yield('bf_10*{run_name}_par8', bf_10, {run_name}_par8)")
    dws.factory(f"expr::Z_z_z_12_yield('bf_12*{run_name}_par9', bf_12, {run_name}_par9)")

    # dws.factory(f"LinearVar::Z_z_z_04_yield(bf_4, {run_name}_par4, 0)")
    # dws.factory(f"LinearVar::Z_z_z_07_yield(bf_7, {run_name}_par5, 0)")
    # dws.factory(f"LinearVar::Z_z_z_08_yield(bf_8, {run_name}_par6, 0)")
    # dws.factory(f"LinearVar::Z_z_z_09_yield(bf_9, {run_name}_par7, 0)")
    # dws.factory(f"LinearVar::Z_z_z_10_yield(bf_10, {run_name}_par8, 0)")
    # dws.factory(f"LinearVar::Z_z_z_12_yield(bf_12, {run_name}_par9, 0)")

    # # dws.factory("expr::n_0710_Z_z_z('n_07_Z_z_z + n_10_Z_z_z', n_07_Z_z_z, n_10_Z_z_z)")
    # # dws.factory("expr::n_040812_Z_z_z('n_04_Z_z_z + n_08_Z_z_z + n_12_Z_z_z',n_04_Z_z_z,n_08_Z_z_z,n_12_Z_z_z)")
    #
    dws.factory(f"expr::P_z_p_02_yield('bf_2*{run_name}_par10', bf_2, {run_name}_par10)")
    dws.factory(f"expr::P_z_p_04_yield('bf_4*{run_name}_par11', bf_4, {run_name}_par11)")
    dws.factory(f"expr::P_z_p_05_yield('bf_5*{run_name}_par12', bf_5, {run_name}_par12)")
    dws.factory(f"expr::P_z_p_06_yield('bf_6*{run_name}_par13', bf_6, {run_name}_par13)")
    dws.factory(f"expr::P_z_p_07_yield('bf_7*{run_name}_par14', bf_7, {run_name}_par14)")
    dws.factory(f"expr::P_z_p_08_yield('bf_8*{run_name}_par15', bf_8, {run_name}_par15)")

    # dws.factory(f"LinearVar::P_z_p_02_yield(bf_2, {run_name}_par10, 0)")
    # dws.factory(f"LinearVar::P_z_p_04_yield(bf_4, {run_name}_par11, 0)")
    # dws.factory(f"LinearVar::P_z_p_05_yield(bf_5, {run_name}_par12, 0)")
    # dws.factory(f"LinearVar::P_z_p_06_yield(bf_6, {run_name}_par13, 0)")
    # dws.factory(f"LinearVar::P_z_p_07_yield(bf_7, {run_name}_par14, 0)")
    # dws.factory(f"LinearVar::P_z_p_08_yield(bf_8, {run_name}_par15, 0)")

    # # dws.factory(f"expr::n_020607_P_z_p('n_02_P_z_p + n_06_P_z_p + n_07_P_z_p', n_02_P_z_p, n_06_P_z_p, n_07_P_z_p)")
    # # dws.factory(f"expr::n_0408_P_z_p('n_04_P_z_p + n_08_P_z_p', n_04_P_z_p, n_08_P_z_p)")
    #
    dws.factory(f"expr::M_m_z_03_yield('bf_3*{run_name}_par16', bf_3, {run_name}_par16)")
    dws.factory(f"expr::M_m_z_04_yield('bf_4*{run_name}_par17', bf_4, {run_name}_par17)")

    # dws.factory(f"LinearVar::M_m_z_03_yield(bf_3, {run_name}_par16, 0)")
    # dws.factory(f"LinearVar::M_m_z_04_yield(bf_4, {run_name}_par17, 0)")

    dws.factory(f"expr::P_z_pst_04_yield('bf_4*{run_name}_par18', bf_4, {run_name}_par18)")
    dws.factory(f"expr::P_z_pst_07_yield('bf_7*{run_name}_par19', bf_7, {run_name}_par19)")
    dws.factory(f"expr::P_z_pst_08_yield('bf_8*{run_name}_par20', bf_8, {run_name}_par20)")

    # dws.factory("expr::n_0408_P_z_pst('n_04_P_z_pst + n_08_P_z_pst', n_04_P_z_pst, n_08_P_z_pst)")

    # dws.factory(f"LinearVar::P_z_pst_04_yield(bf_4, {run_name}_par18, 0)")
    # dws.factory(f"LinearVar::P_z_pst_07_yield(bf_7, {run_name}_par19, 0)")
    # dws.factory(f"LinearVar::P_z_pst_08_yield(bf_8, {run_name}_par20, 0)")

    # dws.factory(f"expr::Zs_sm_p_13_yield('bf_13*{run_name}_par21', bf_13, {run_name}_par21)")
    # dws.factory(f"expr::Zs_sm_p_14_yield('bf_14*{run_name}_par22', bf_14, {run_name}_par22)")
    # dws.factory(f"expr::Zs_sm_p_15_yield('bf_15*{run_name}_par23', bf_15, {run_name}_par23)")
    # dws.factory(f"expr::Zs_sm_p_16_yield('bf_16*{run_name}_par24', bf_16, {run_name}_par24)")

    # dws.factory(f"LinearVar::Zs_sm_p_13_yield(bf_13, {run_name}_par21, 0)")
    # dws.factory(f"LinearVar::Zs_sm_p_14_yield(bf_14, {run_name}_par22, 0)")
    # dws.factory(f"LinearVar::Zs_sm_p_15_yield(bf_15, {run_name}_par23, 0)")
    # dws.factory(f"LinearVar::Zs_sm_p_16_yield(bf_16, {run_name}_par24, 0)")
    dws.factory(f"width_P_z_pst_gsmear[1, 0.0001, 2]")
    gsmear_strat = "ALL"
    if gsmear_strat == "SPLIT":
        for spec in specs:
            if spec != "P_z_pst":
                dws.factory(f"width_{spec}_gsmear[1, 0.001, 5]")
            dws.factory(f"mean_{spec}_gsmear[0]")
            dws.factory(f"Gaussian::{spec}_gsmear(B_DTF_M, mean_{spec}_gsmear, width_{spec}_gsmear)")
    if gsmear_strat == "ALL":
        dws.factory(f"width_gsmear[1, 0.001, 5]")
        dws.factory(f"mean_gsmear[0]")
        dws.factory(f"Gaussian::a_gsmear(B_DTF_M, mean_gsmear, width_gsmear)")
        dws.factory(f"Gaussian::P_z_pst_gsmear(B_DTF_M, mean_gsmear, width_P_z_pst_gsmear)")


    for spec in specs:
        get_shapes_bkg(spec, "Exponential", dws)
        # get_mc_shape_corr(dws, spec, split_flag = True, fix_flag = True, smear_flag = True, gsmear_strat = "SPLIT")
        get_mc_shape_corr(dws, spec, split_flag = True, fix_flag = True, smear_flag = True, gsmear_strat = "ALL")

    dws.factory("SUM::Z_m_p_spectrum_all_fit(Z_m_p_01_yield*Z_m_p_01_fit, Z_m_p_02_yield*Z_m_p_02_fit, Z_m_p_03_yield*Z_m_p_03_fit, Z_m_p_04_yield*Z_m_p_04_fit, Z_m_p_bkg_yield[100,0,100000]* Z_m_p_spectrum_bkg)")
    dws.factory("SUM::Z_z_z_spectrum_all_fit(Z_z_z_09_yield*Z_z_z_09_fit, Z_z_z_07_yield*Z_z_z_07_fit, Z_z_z_10_yield*Z_z_z_10_fit, Z_z_z_04_yield*Z_z_z_04_fit, Z_z_z_08_yield*Z_z_z_08_fit, Z_z_z_12_yield*Z_z_z_12_fit, Z_z_z_bkg_yield[100,0,100000]*Z_z_z_spectrum_bkg)")
    dws.factory("SUM::P_z_p_spectrum_all_fit(P_z_p_05_yield*P_z_p_05_fit,P_z_p_02_yield*P_z_p_02_fit,P_z_p_06_yield*P_z_p_06_fit,P_z_p_07_yield*P_z_p_07_fit,P_z_p_04_yield*P_z_p_04_fit,P_z_p_08_yield*P_z_p_08_fit,P_z_p_bkg_yield[100,0,10000]*P_z_p_spectrum_bkg)")
    dws.factory("SUM::M_m_z_spectrum_all_fit(M_m_z_03_yield*M_m_z_03_fit,M_m_z_04_yield*M_m_z_04_fit,M_m_z_bkg_yield[100,0,10000]*M_m_z_spectrum_bkg)")
    dws.factory("SUM::P_z_pst_spectrum_all_fit(P_z_pst_07_yield*P_z_pst_07_fit,P_z_pst_04_yield*P_z_pst_04_fit,P_z_pst_08_yield*P_z_pst_08_fit,P_z_pst_bkg_yield[100,0,1000]*P_z_pst_spectrum_bkg)")
    # dws.factory("SUM::Zs_sm_p_spectrum_all_fit(Zs_sm_p_13_yield*Zs_sm_p_13_fit, Zs_sm_p_14_yield*Zs_sm_p_14_fit, Zs_sm_p_15_yield*Zs_sm_p_15_fit, Zs_sm_p_16_yield*Zs_sm_p_16_fit, Zs_sm_p_bkg_yield[100,0,10000]*Zs_sm_p_spectrum_bkg)")

    all_cats = ROOT.RooCategory("all_cats", "all_cats")
    dspectrum_list = ["Z_m_p_spectrum", "Z_z_z_spectrum", "P_z_p_spectrum","M_m_z_spectrum", "P_z_pst_spectrum"]

    data_set_list = []

    for dspec in dspectrum_list:
        all_cats.defineType(dspec)
        pathname = dspec.split("_spectrum")[0]
        file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{pathname}.root")
        tchain = ROOT.TChain("DecayTreeTuple_SIG")
        for file_name in file_list:
            tchain.Add(file_name)
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
    )

    getattr(dws, "import")(all_cats)
    getattr(dws, "import")(all_data_sets)

    if gc_onflag == 1:
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
        glist = ROOT.RooArgSet()
        for i in range(0, lp):
            gvar = dws.pdf(f"{run_name}_nu{i}_g")
            glist.add(gvar)

    all_fit = ROOT.RooSimultaneous("super_fit_Pdf", "super_fit_Pdf", all_cats)
    z_model = dws.pdf("Z_m_p_spectrum_all_fit")
    all_fit.addPdf(z_model, "Z_m_p_spectrum")
    zz_model = dws.pdf("Z_z_z_spectrum_all_fit")
    all_fit.addPdf(zz_model, "Z_z_z_spectrum")
    p_model = dws.pdf("P_z_p_spectrum_all_fit")
    all_fit.addPdf(p_model, "P_z_p_spectrum")
    m_model = dws.pdf("M_m_z_spectrum_all_fit")
    all_fit.addPdf(m_model, "M_m_z_spectrum")
    st_model = dws.pdf("P_z_pst_spectrum_all_fit")
    all_fit.addPdf(st_model, "P_z_pst_spectrum")
    # s_model = dws.pdf("Zs_sm_p_spectrum_all_fit")
    # all_fit.addPdf(s_model, "Zs_sm_p_spectrum")

    if gc_onflag == 1:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.ExternalConstraints(glist), ROOT.RooFit.Save())
    if gc_onflag == 0:
        nall_fit = all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

    # dws.Print()

    file = open("bf.txt","w")
    for i in range(1,13):
        if i != 11:
            bf = dws.var(f"bf_{i}")
            bf_val = bf.getValV()
            bf_err = bf.getPropagatedError(nall_fit)
            bf_print = ufloat(bf_val, bf_err)
            string = bf_string_dict[str(i)]
            file.write(f"{string} : {bf_val:.2E} {bf_err:.2E} \n")
    file.close()

    dws.Import(nall_fit)
    dws.Import(all_fit)

    output_base_data = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/corr_{run_name}.root"
    dws.writeToFile(output_base_data)
    print(f"Wrote dws to: {output_base_data}")

def plot_nn_data(run_name, spec, split_flag, varoi = "B_DTF_M"):

    fws_base_plot_file = ROOT.TFile(f"{analysis_path}/data_fit/data_files/{spec}_{run_name}_{varoi}.root")
    fws = fws_base_plot_file.Get(f"{run_name}")

    b_dtf_m = fws.var(f"{varoi}")
    fit_pdf = fws.pdf(f"{spec}_spectrum_all_fit")
    data_set = fws.data(f"{spec}_final_data")
    fitresult = fws.obj(f"fitresult_{spec}_spectrum_all_fit_{spec}_final_data")
    frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}_spectrum"))

    if spec == "Z_m_p":
        ylist = ["01","0203","04"]
        title = "D^{-} D^{+} K^{*0}"
        d1 = "D^{-}"
        d2 = "D^{+}"
    if spec == "Z_z_z":
        if not split_flag:
            ylist = ["09","0710","040812"]
        if split_flag:
            ylist = ["09","07", "10","04","08","12"]
        title = "#bar{D^{0}} D^{0} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{0}"
    if spec == "P_z_p":
        if not split_flag:
            ylist = ["05","020607","0408"]
        if split_flag:
            ylist = ["05","02","06","07","04","08"]
        title = "#bar{D^{0}} D^{+} K^{*0}"
        d1 = "#bar{D^{0}}"
        d2 = "D^{0}"
    if spec == "M_m_z":
        ylist = ["03","04"]
        title = "D^{-} D^{0} K^{*0}"
    if spec == "P_z_pst":
        title = "#bar{D^{0}} (D^{*+} #rightarrow     D^{0} #pi^{+}) K^{*0}"
        if not split_flag:
            ylist = ["07","0408"]
        if split_flag:
            ylist = ["07","04","08"]
    if spec == "Zs_sm_p":
        ylist = ["13","14","15","16"]
        title = "D^{-} D^{0} K^{*0}"
    if spec == "norm7":
        ylist = ["norm7"]
        title = "D^{-} D^{0} K^{*0}"
    if spec == "norm8":
        ylist = ["norm8"]
        title = "D^{-} D^{0} K^{*0}"

    data_set.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.MarkerSize(1))
    fit_pdf.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kSolid),  ROOT.RooFit.LineColor(sig_col),ROOT.RooFit.LineWidth(3),  ROOT.RooFit.Name("pdf"))
    fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(bkg_col), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name("bkg"))
    color_list = [p0_col, p1_col, p2_col]
    for y, color in zip(ylist,color_list):
        if "norm" not in spec:
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_{y}_fit"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(color), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name(f"fit_{y}"))
        if "norm" in spec:
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit"), ROOT.RooFit.LineStyle(2), ROOT.RooFit.LineColor(color), ROOT.RooFit.LineWidth(3), ROOT.RooFit.Name(f"fit_{y}"))
    p = ROOT.TCanvas("p1","p1")
    p.cd()

    d_chi2 = frame.chiSquare(f"pdf", f"data")

    if spec == "Z_m_p":
        legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.905)
        name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield","Missing 2 Particle Yield"]

    if spec == "Z_z_z":
        legend = ROOT.TLegend(0.60, 0.55, 0.85, 0.85)
        name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield","Missing 2 Particle Yield"]

    if spec == "P_z_p":
        legend = ROOT.TLegend(0.65, 0.625, 0.90, 0.925)
        name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield","Missing 2 Particle Yield"]

    if spec == "M_m_z":
        legend = ROOT.TLegend(0.60, 0.55, 0.85, 0.85)
        name_list =  ["Missing 1 Particle Yield","Missing 2 Particle Yield"]

    if spec == "P_z_pst":
        legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.875)
        name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield"]

    if spec == "Zs_sm_p":
        legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.875)
        name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield", "Missing 1 Particle Yield", "Missing 2 Particle Yield"]

    legend.SetFillStyle(1001)

    d_nyield_bkg = fws.obj(f"{spec}_bkg_yield")
    d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

    legend.AddEntry(frame.findObject("data"),"Run 2 Data","ep")
    legend.AddEntry(frame.findObject("pdf"),"Total Fit PDF","l")
    legend.AddEntry(frame.findObject("bkg"),f"Background PDF : {round(d_nyield_bkg.getValV(),1)} #pm {round(d_nyield_bkg_err,1)}","l")
    for y, name in zip(ylist, name_list):
        d_nyield = fws.obj(f"{spec}_{y}_yield")
        d_nyield_err = d_nyield.getPropagatedError(fitresult)
        legend.AddEntry(frame.findObject(f"fit_{y}"), f"{name}: {round(d_nyield.getValV(),1)} #pm {round(d_nyield_err,1)}", "l")
    legend.SetTextSize(0.025)

    p = rootutils.residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)

    frame.Draw()
    legend.Draw()
    p.pb.cd()

    hpull = frame.pullHist(f"data", f"pdf")
    pull = b_dtf_m.frame()
    pull.addPlotable(hpull, "P")
    pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
    pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize() * p.blabelratio)
    pull.GetYaxis().SetLabelSize(
        ROOT.gStyle.GetLabelSize("Y") * (1 + p.padratio) / (2 * p.padratio)
    )
    pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
    pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
    pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset() / p.blabelratio)
    pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength() * p.blabelratio)
    pull.GetYaxis().SetTitle("Pull")
    pull.GetYaxis().SetNdivisions(202, False)
    pull.GetYaxis().SetRangeUser(-3, 3)
    pull.GetYaxis().CenterTitle()
    pull.Draw("AP")

    pull.GetXaxis().SetTitle(f" m({title}) [MeV]")
    # frame.GetXaxis().SetTitle(f" m({title}) [MeV]")

    save_pdf(p, f"fit_tests_{varoi}", f"{spec}_{run_name}_{varoi}", rpflag = 1)

def plot_coor_data(run_name, specs):
    fws_base_plot_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/data_files/corr_{run_name}.root")
    fws = fws_base_plot_file.Get(f"{run_name}")

    b_dtf_m = fws.var("B_DTF_M")
    all_cats = fws.cat("all_cats")
    all_data_sets = fws.data("all_data_sets")
    all_fit = fws.pdf("super_fit_Pdf")
    all_cats_Set = ROOT.RooArgSet(all_cats)
    fitresult = fws.obj("fitresult_super_fit_Pdf_all_data_sets")

    for spec in specs:

        spectrum = f"{spec}_spectrum"
        #ROOT.RooFit.Bins(bbins)
        frame = b_dtf_m.frame(ROOT.RooFit.Title(f"{spec}_spectrum"))

        if spec == "Z_m_p":
            ylist = ["01","02","03","04"]
            title = "D^{-} D^{+} K^{*0}"
            d1 = "D^{-}"
            d2 = "D^{+}"
        if spec == "Z_z_z":
            ylist = ["09","07", "10","04","08","12"]
            title = "#bar{D^{0}} D^{0} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{0}"
        if spec == "P_z_p":
            ylist = ["05","02","06","07","04","08"]
            title = "#bar{D^{0}} D^{+} K^{*0}"
            d1 = "#bar{D^{0}}"
            d2 = "D^{0}"
        if spec == "M_m_z":
            ylist = ["03","04"]
            title = "D^{-} D^{0} K^{*0}"
        if spec == "P_z_pst":
            title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi^{+}) K^{*0}"
            ylist = ["07","04","08"]
        if spec == "Zs_sm_p":
            title = "D_{s}^{-} D^{+} K^{*0}"
            ylist = ["13","14","15","16"]


        all_data_sets.plotOn(frame, ROOT.RooFit.Cut(f"all_cats==all_cats::{spec}_spectrum"), ROOT.RooFit.Name("data"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_all_fit"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("pdf"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_spectrum_bkg"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name("bkg"))
        for y, color in zip(ylist, [ROOT.kBlue, ROOT.kCyan, ROOT.kViolet, ROOT.kAzure, ROOT.kMagenta, ROOT.kGreen]):
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, spectrum), ROOT.RooFit.Components(f"{spec}_{y}_fit"), ROOT.RooFit.ProjWData(all_cats_Set, all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(color), ROOT.RooFit.Name(f"fit_{y}"))


        p = ROOT.TCanvas("p1","p1")
        p.cd()
        frame.GetXaxis().SetTitle(f" m({title}) [MeV]")

        if spec == "Z_m_p":
            legend = ROOT.TLegend(0.20, 0.5, 0.40, 0.900)
            name_list =  [bf_string_dict["1"], bf_string_dict["2"], bf_string_dict["3"], bf_string_dict["4"]]

        if spec == "Z_z_z":
            legend = ROOT.TLegend(0.515, 0.40, 0.95, 0.85)
            name_list = [bf_string_dict["9"], bf_string_dict["7"], bf_string_dict["10"], bf_string_dict["4"], bf_string_dict["8"], bf_string_dict["12"]]

        if spec == "P_z_p":
            legend = ROOT.TLegend(0.65, 0.425, 0.90, 0.825)
            name_list = [bf_string_dict["5"], bf_string_dict["2"], bf_string_dict["6"], bf_string_dict["7"], bf_string_dict["4"], bf_string_dict["8"]]

        if spec == "M_m_z":
            legend = ROOT.TLegend(0.60, 0.55, 0.85, 0.85)
            name_list =  [bf_string_dict["3"], bf_string_dict["4"]]

        if spec == "P_z_pst":
            legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.875)
            name_list =  [bf_string_dict["7"], bf_string_dict["4"], bf_string_dict["8"]]

        if spec == "Zs_sm_p":
            legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.875)
            name_list =  [bf_string_dict["13"], bf_string_dict["14"], bf_string_dict["15"], bf_string_dict["16"]]


        legend.SetTextSize(0.0375)
        # legend.SetEntrySeparation(1)


        legend.SetFillStyle(1001)

        d_nyield_bkg = fws.obj(f"{spec}_bkg_yield")
        d_nyield_bkg_err = d_nyield_bkg.getPropagatedError(fitresult)

        legend.AddEntry(frame.findObject("data"),"Run 2 Data","ep")
        legend.AddEntry(frame.findObject("pdf"),"Total Fit PDF","l")
        legend.AddEntry(frame.findObject("bkg"),f"Background PDF","l")
        #{round(d_nyield_bkg.getValV(),1)} #pm {round(d_nyield_bkg_err,1)}
        for y, name in zip(ylist, name_list):
            d_nyield = fws.obj(f"{spec}_{y}_yield")
            d_nyield_err = d_nyield.getPropagatedError(fitresult)
            legend.AddEntry(frame.findObject(f"fit_{y}"), f"{name}", "l")
            #: {round(d_nyield.getValV(),1)} #pm {round(d_nyield_err,1)}"

        frame.GetXaxis().SetTitle(f" m({title}) [MeV]")

        frame.addObject(legend)
        frame.Draw()
        legend.Draw()
        save_pdf(p, f"fit_tests_Corr", f"corr_{run_name}_{spec}", rpflag = 0)
