import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

def get_pdg_values(dws):
    xls = ExcelFile("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_mc_base.xlsx")
    for sheet in xls.sheet_names:
        if sheet == "BF_Values_PDG":
            block = xls.parse(sheet)
            dict = block.to_dict()
            for myid, bf_value, err_bf_value in zip(dict["myid"].values(),dict["bf_value"].values(),dict["err_bf_value"].values()):
                if gc_onflag == 1:
                    if err_bf_value != 0:
                        dws.factory(f"{myid}[{bf_value}, 0, 1]")
                        dws.factory(f"Gaussian::{myid}_g({myid}, {bf_value}, {err_bf_value})")
                else:
                    dws.factory(f"{myid}[{bf_value}]")
                print(myid, bf_value, err_bf_value)
def get_eff_values(dws, sheetname):
    xls = ExcelFile("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/spreadsheet/b20c_2021_eff.xlsx")
    for sheet in xls.sheet_names:
        if sheet == sheetname:
            block = xls.parse(sheet)
            dict = block.to_dict()
            for myid, eff, eff_err in zip(dict["my_eff_id"].values(),dict["total_eff"].values(),dict["err_total_eff"].values()):
                if pandas.isna(id) != True and pandas.isna(eff) != True and pandas.isna(eff_err) != True:
                    if gc_onflag == 1:
                        dws.factory(f"{myid}[{eff}, 0, 1]")
                        dws.factory(f"Gaussian::{myid}_g({myid}, {eff}, {eff_err})")
                        temppdf = dws.pdf(f"{myid}_g")
                    else:
                        dws.factory(f"{myid}[{eff}]")
                print(myid, eff, eff_err)
                
def build_data_ws(name, eff_sheet_name, tt_flag):

    dws = ROOT.RooWorkspace(name)
    get_pdg_values(dws)
    get_eff_values(dws, eff_sheet_name)

    norm7_ws_file = ROOT.TFile(f"normalization/norm7_{tt_flag}.root","READ")
    norm8_ws_file = ROOT.TFile(f"normalization/norm8_{tt_flag}.root","READ")

    norm7_ws = norm7_ws_file.Get(f"norm7_{tt_flag}")
    norm8_ws = norm8_ws_file.Get(f"norm8_{tt_flag}")

    getattr(dws, 'import')(norm7_ws.var("n_norm_signal"), ROOT.RooFit.RenameVariable("n_norm_signal", "n_norm7"))
    getattr(dws, 'import')(norm8_ws.var("n_norm_signal"), ROOT.RooFit.RenameVariable("n_norm_signal", "n_norm8"))

    dws.factory(f"B_DTF_M[{bmin},{bmax}]")
    dws.factory("D1_M[0,6000]")
    dws.factory("D2_M[0,6000]")
    dws.factory("D1_DIRA_ORIVX[-5,5]")
    dws.factory("D2_DIRA_ORIVX[-5,5]")
    dws.factory("D1_FDCHI2_ORIVX[0,1000]")
    dws.factory("D2_FDCHI2_ORIVX[0,1000]")
    dws.factory("D2st_M[0,10000]")

    dws.factory("bf_1[0.00001,0,1]")
    dws.factory("bf_2[0.00001,0,1]")
    dws.factory("bf_3[0.00001,0,1]")
    dws.factory("bf_4[0.00001,0,1]")
    dws.factory("bf_5[0.00001,0,1]")
    dws.factory("bf_6[0.00001,0,1]")
    dws.factory("bf_7[0.00001,0,1]")
    dws.factory("bf_8[0.00001,0,1]")
    dws.factory("bf_9[0.00001,0,1]")
    dws.factory("bf_10[0.00001,0,1]")
    dws.factory("bf_11[0.00001,0,1]")
    dws.factory("bf_12[0.00001,0,1]")
    dws.factory("bf_13[0.00001,0,1]")
    dws.factory("bf_14[0.00001,0,1]")
    dws.factory("bf_15[0.00001,0,1]")
    dws.factory("bf_16[0.00001,0,1]")
    dws.factory("eff_spi[0.3, 0, 1]")
    #
    dws.factory("expr::norm8_factor('(0.66666*n_norm8)/(b0norm8*eff_norm8*dm*d0bar4)', n_norm8, b0norm8, eff_norm8, dm, d0bar4)")
    dws.factory("expr::norm7_factor('(0.66666*n_norm7)/(bpnorm7*eff_norm7*d0bar*d0bar4)', n_norm7, bpnorm7, eff_norm7, d0bar, d0bar4)")

    dws.factory("expr::n_01_z_base('bf_1*dm*dm*eff_01', bf_1, dm, eff_01)")

    dws.factory("expr::n_02_z_base('bf_2*dstmdmpi0*dm*dm*eff_2a', bf_2, dstmdmpi0, dm, eff_2a)")
    dws.factory("expr::n_02_p_base('bf_2*dstmd0pim*d0bar*dm*eff_2b', bf_2, dstmd0pim, d0bar, dm, eff_2b)")

    dws.factory("expr::n_03_z_base('bf_3*dstmdmpi0*dm*dm*eff_3a', bf_3, dstmdmpi0, dm, eff_3a)")
    dws.factory("expr::n_03_m_base('bf_3*dstmdmpi0*d0bar*dm*eff_3b', bf_3, dstmdmpi0, d0bar, dm, eff_3b)")

    dws.factory("expr::n_04_z_base('bf_4*dstmdmpi0*dstmdmpi0*dm*dm*eff_4a', bf_4, dstmdmpi0, dm, eff_4a)")
    dws.factory("expr::n_04_p_base('bf_4*dstmdmpi0*dstmd0pim*dm*d0bar*eff_4b', bf_4, dstmdmpi0, dstmd0pim, dm, d0bar, eff_4b)")
    dws.factory("expr::n_04_m_base('bf_4*dstmdmpi0*dstmd0pim*dm*d0bar*eff_4c', bf_4, dstmdmpi0, dstmd0pim, dm, d0bar, eff_4c)")
    dws.factory("expr::n_04_zz_base('bf_4*dstmd0pim*dstmd0pim*d0bar*d0bar*eff_4d', bf_4, dstmd0pim, d0bar, eff_4d)")
    dws.factory("expr::n_04_st_base('bf_4*dstmd0pim*dstmd0pim*d0bar*d0bar*eff_4d*eff_spi', bf_4, dstmd0pim, d0bar, eff_4d, eff_spi)")

    dws.factory("expr::n_05_p_base('bf_5*d0bar*dm*eff_05', bf_5, d0bar, dm, eff_05)")
    dws.factory("expr::n_06_p_base('bf_6*d0bar*dm*eff_06', bf_6, d0bar, dm, eff_06)")

    dws.factory("expr::n_07_p_base('bf_7*dstmdmpi0*d0bar*dm*eff_7a', bf_7, dstmdmpi0, d0bar, dm, eff_7a)")
    dws.factory("expr::n_07_zz_base('bf_7*dstmd0pim*d0bar*d0bar*eff_7b', bf_7, dstmd0pim, d0bar, eff_7b)")
    dws.factory("expr::n_07_st_base('bf_7*dstmd0pim*d0bar*d0bar*eff_7b*eff_spi', bf_7, dstmd0pim, d0bar, eff_7b, eff_spi)")

    dws.factory("expr::n_08_p_base('bf_8*dstmdmpi0*d0bar*dm*eff_8a', bf_8, dstmdmpi0, d0bar, dm, eff_8a)")
    dws.factory("expr::n_08_zz_base('bf_8*dstmd0pim*d0bar*d0bar*eff_8b', bf_8, dstmd0pim, d0bar, eff_8b)")
    dws.factory("expr::n_08_st_base('bf_8*dstmd0pim*d0bar*d0bar*eff_8b*eff_spi', bf_8, dstmd0pim, d0bar, eff_8b, eff_spi)")

    dws.factory("expr::n_09_zz_base('bf_9*d0bar*d0bar*eff_09', bf_9, d0bar, eff_09)")
    dws.factory("expr::n_10_zz_base('bf_10*d0bar*d0bar*eff_10', bf_10, d0bar, eff_10)")
    dws.factory("expr::n_11_zz_base('bf_11*d0bar*d0bar*eff_11', bf_11, d0bar, eff_11)")
    dws.factory("expr::n_12_zz_base('bf_12*d0bar*d0bar*eff_12', bf_12, d0bar, eff_12)")

    dws.factory("expr::n_13_base('bf_13*dsm*dm*eff_13', bf_13, dsm, dm, eff_13)")
    dws.factory("expr::n_14_base('bf_14*dsm*dm*eff_14', bf_14, dsm, dm, eff_14)")
    dws.factory("expr::n_15_base('bf_15*dstmdmpi0*dsm*dm*eff_15', bf_15, dstmdmpi0, dsm, dm, eff_15)")
    dws.factory("expr::n_16_base('bf_16*dstmdmpi0*dsm*dm*eff_16', bf_16, dstmdmpi0, dsm, dm, eff_16)")

    dws.factory("expr::n_01_z('n_01_z_base*norm8_factor', n_01_z_base, norm8_factor)")

    dws.factory("expr::n_02_z('n_02_z_base*norm8_factor', n_02_z_base, norm8_factor)")
    dws.factory("expr::n_03_z('n_03_z_base*norm8_factor', n_03_z_base, norm8_factor)")
    dws.factory("expr::n_23_z('n_02_z + n_03_z', n_02_z, n_03_z)")

    dws.factory("expr::n_04_z('n_04_z_base*norm8_factor', n_04_z_base, norm8_factor)")

    dws.factory("expr::n_09_zz('n_09_zz_base*norm7_factor', n_09_zz_base, norm7_factor)")

    dws.factory("expr::n_07_zz('n_07_zz_base*norm7_factor', n_07_zz_base, norm7_factor)")
    dws.factory("expr::n_10_zz('n_10_zz_base*norm7_factor', n_10_zz_base, norm7_factor)")
    dws.factory("expr::n_11_zz('n_11_zz_base*norm7_factor', n_11_zz_base, norm7_factor)")
    dws.factory("expr::n_71011_zz('n_07_zz + n_10_zz + n_11_zz', n_07_zz, n_10_zz, n_11_zz)")

    dws.factory("expr::n_04_zz('n_04_zz_base*norm7_factor', n_04_zz_base, norm7_factor)")
    dws.factory("expr::n_08_zz('n_08_zz_base*norm7_factor', n_08_zz_base, norm7_factor)")
    dws.factory("expr::n_12_zz('n_12_zz_base*norm7_factor', n_12_zz_base, norm7_factor)")
    dws.factory("expr::n_4812_zz('n_04_zz + n_08_zz + n_12_zz',n_04_zz,n_08_zz,n_12_zz)")

    dws.factory("expr::n_05_p('n_05_p_base*norm7_factor', n_05_p_base, norm7_factor)")

    dws.factory("expr::n_02_p('n_02_p_base*norm7_factor', n_02_p_base, norm7_factor)")
    dws.factory("expr::n_06_p('n_06_p_base*norm7_factor', n_06_p_base, norm7_factor)")
    dws.factory("expr::n_07_p('n_07_p_base*norm7_factor', n_07_p_base, norm7_factor)")
    dws.factory("expr::n_267_p('n_02_p + n_06_p + n_07_p', n_02_p, n_06_p, n_07_p)")

    dws.factory("expr::n_04_p('n_04_p_base*norm7_factor', n_04_p_base, norm7_factor)")
    dws.factory("expr::n_08_p('n_08_p_base*norm7_factor', n_08_p_base, norm7_factor)")
    dws.factory("expr::n_48_p('n_04_p + n_08_p', n_04_p, n_08_p)")

    dws.factory("expr::n_03_m('n_03_m_base*norm7_factor',n_03_m_base, norm7_factor)")
    dws.factory("expr::n_04_m('n_04_m_base*norm7_factor',n_04_m_base, norm7_factor)")

    dws.factory("expr::n_07_st('n_07_st_base*norm8_factor', n_07_st_base, norm8_factor)")

    dws.factory("expr::n_04_st('n_04_st_base*norm8_factor', n_04_st_base, norm8_factor)")
    dws.factory("expr::n_08_st('n_08_st_base*norm8_factor', n_08_st_base, norm8_factor)")
    dws.factory("expr::n_48_st('n_04_st + n_08_st', n_04_st, n_08_st)")

    dws.factory("expr::n_13_s('n_13_base*norm8_factor',n_13_base, norm8_factor)")
    dws.factory("expr::n_14_s('n_14_base*norm8_factor',n_14_base, norm8_factor)")
    dws.factory("expr::n_15_s('n_15_base*norm8_factor',n_15_base, norm8_factor)")
    dws.factory("expr::n_1415_s('n_14_s + n_15_s', n_14_s, n_15_s)")
    dws.factory("expr::n_16_s('n_16_base*norm8_factor',n_16_base, norm8_factor)")

    b_dtf_m = dws.var("B_DTF_M")
    b_dtf_m.setRange("myrange", bmin, bmax)
    d1_mass = dws.var("D1_M")
    d2_mass = dws.var("D2_M")
    d1_dira = dws.var("D1_DIRA_ORIVX")
    d2_dira = dws.var("D2_DIRA_ORIVX")
    d1_fdx2 = dws.var("D1_FDCHI2_ORIVX")
    d2_fdx2 = dws.var("D2_FDCHI2_ORIVX")
    d2st_mass = dws.var("D2st_M")

    data_args = ROOT.RooArgSet(d1_mass, d2_mass, d1_dira, d2_dira, d1_fdx2, d2_fdx2, b_dtf_m)
    st_data_args = ROOT.RooArgSet(d1_mass, d2_mass, d1_dira, d2_dira, d1_fdx2, d2_fdx2, b_dtf_m, d2st_mass)

    all_cats = ROOT.RooCategory("all_cats","all_cats")
    all_cats.defineType("z_spectrum")
    all_cats.defineType("zz_spectrum")
    all_cats.defineType("p_spectrum")
    all_cats.defineType("m_spectrum")
    all_cats.defineType("st_spectrum")
    all_cats.defineType("s_spectrum")

    dwindow_file = ROOT.TFile(f"dmass/d_mass_fits.root","READ")
    dwindow_ws = dwindow_file.Get("d_mass_fits")

    z_m_cut = get_dwindow_values(dwindow_ws, "z")
    zz_m_cut = get_dwindow_values(dwindow_ws, "zz")
    p_m_cut = get_dwindow_values(dwindow_ws, "p")
    m_m_cut = get_dwindow_values(dwindow_ws, "m")
    st_m_cut = get_dwindow_values(dwindow_ws, "st")

    # dws.Print()
    print("done")
    if tt_flag == "ToT":
        z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
        zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
        p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")
        m_data_file = ROOT.TFile(data_basepath+"m_spectrum.root")
        st_data_file = ROOT.TFile(data_basepath+"st_spectrum_newfd.root")
        # s_data_file = ROOT.TFile(data_basepath+"s_spectrum.root")
        z_tree = z_data_file.Get("DecayTreeTuple")
        zz_tree = zz_data_file.Get("DecayTreeTuple")
        p_tree = p_data_file.Get("DecayTreeTuple")
        m_tree = m_data_file.Get("DecayTreeTuple")
        st_tree = st_data_file.Get("DecayTreeTuple")
        # s_tree = s_data_file.Get("DecayTreeTuple")
    if tt_flag == "T" or tt_flag == "nTaT":
        z_data_file = ROOT.TFile(TT_data_basepath+"TT_z_spectrum.root")
        zz_data_file = ROOT.TFile(TT_data_basepath+"TT_zz_spectrum.root")
        p_data_file = ROOT.TFile(TT_data_basepath+"TT_p_spectrum.root")
        m_data_file = ROOT.TFile(TT_data_basepath+"TT_m_spectrum.root")
        s_data_file = ROOT.TFile(TT_data_basepath+"TT_s_spectrum.root")
        st_data_file_t = ROOT.TFile(TT_data_basepath+"st_spectrum_newfd_T.root")
        st_data_file_ntat= ROOT.TFile(TT_data_basepath+"st_spectrum_newfd_nTaT.root")
    if tt_flag == "T":
        z_tree = z_data_file.Get("DecayTreeTuple_T")
        zz_tree = zz_data_file.Get("DecayTreeTuple_T")
        p_tree = p_data_file.Get("DecayTreeTuple_T")
        m_tree = m_data_file.Get("DecayTreeTuple_T")
        st_tree = st_data_file_t.Get("DecayTreeTuple")
        # s_tree = s_data_file.Get("DecayTreeTuple_T")
    if tt_flag == "nTaT":
        z_tree = z_data_file.Get("DecayTreeTuple_nTaT")
        zz_tree = zz_data_file.Get("DecayTreeTuple_nTaT")
        p_tree = p_data_file.Get("DecayTreeTuple_nTaT")
        m_tree = m_data_file.Get("DecayTreeTuple_nTaT")
        st_tree = st_data_file_ntat.Get("DecayTreeTuple")
        # s_tree = s_data_file.Get("DecayTreeTuple_nTaT")
    #
    z_data_set = ROOT.RooDataSet("z_data","z_data", z_tree, data_args, z_m_cut)
    zz_data_set = ROOT.RooDataSet("zz_data","zz_data", zz_tree, data_args, zz_m_cut)
    p_data_set = ROOT.RooDataSet("p_data","p_data", p_tree, data_args, p_m_cut)
    m_data_set = ROOT.RooDataSet("m_data","m_data", m_tree, data_args, m_m_cut)
    st_data_set = ROOT.RooDataSet("st_data","st_data", st_tree, st_data_args, st_m_cut)
    # s_data_set = ROOT.RooDataSet("s_data","s_data", s_tree, data_args, s_m_cut)

    all_data_sets = ROOT.RooDataSet("all_data_sets", "all_data_sets", data_args, ROOT.RooFit.Index(all_cats),
                ROOT.RooFit.Import(
                            "z_spectrum",
                            z_data_set),
                ROOT.RooFit.Import(
                            "zz_spectrum",
                            zz_data_set),
                ROOT.RooFit.Import(
                            "p_spectrum",
                            p_data_set),
                ROOT.RooFit.Import(
                            "m_spectrum",
                            m_data_set),
                ROOT.RooFit.Import(
                            "st_spectrum",
                            st_data_set),
                # ROOT.RooFit.Import(
                #             "s_spectrum",
                #             s_data_set),
                   )

    getattr(dws,'import')(all_cats)
    getattr(dws,'import')(all_data_sets)

    dws.writeToFile(f"signal/{name}.root")
    dws.Print()


build_data_ws("dws_ToT_test", "MC_eff_ToT", "ToT")
# build_data_ws("dws_T", "MC_eff_T", "T")
# build_data_ws("dws_nTaT", "MC_eff_nTaT", "nTaT")
