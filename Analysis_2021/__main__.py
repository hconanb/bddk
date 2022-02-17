import os
import sys
import datetime
import ROOT
from uncertainties import ufloat
from uncertainties import covariance_matrix, correlation_matrix
from uncertainties import unumpy

import xlrd
import shutil
import pprint

from pandas import *
from pandas.api.types import CategoricalDtype

import numpy as np
import math as math
from array import array

import root_pandas as rp
import itertools
RDF = ROOT.ROOT.RDataFrame

import ast
import argparse
from argparse import RawTextHelpFormatter

ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")
ROOT.gStyle.SetPalette(ROOT.kBird)

from Analysis_2021 import rootutils
from Analysis_2021.data_fit import data_fit_v1
from Analysis_2021.mc_fit import mc_ws_build_script, mc_res_test
from Analysis_2021.dalitz_reweighting import bdt_test

if __name__ == '__main__':

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument('--run_name', default = "test")
    parser.add_argument('--specs', default = ["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"], nargs='+')

    parser.add_argument('--data_build', action='store_true')
    parser.add_argument('--data_plot',action='store_true')
    parser.add_argument('--split_flag',action='store_true')

    parser.add_argument('--data_corr_build', action = 'store_true')
    parser.add_argument('--data_corr_plot', action = 'store_true')
    parser.add_argument('--data_post_d_plot', action = 'store_true')

    parser.add_argument('--fit_strat', default = "dp_f")
    parser.add_argument('--gc_onflag', default = 0, type=int)
    parser.add_argument('--fix_flag', action='store_true')
    parser.add_argument('--smear_flag', action='store_true')

    parser.add_argument('--norm_build',action = 'store_true')
    parser.add_argument('--norm_plot',action = 'store_true')

    parser.add_argument('--mc_build', action='store_true')
    parser.add_argument('--mc_data_comp', action = 'store_true')

    parser.add_argument('--mc_build_pizy', action='store_true')


    parser.add_argument('--mc_fit', action='store_true')
    parser.add_argument('--mc_plot',action='store_true')
    parser.add_argument('--mc_res',action='store_true')


    parser.add_argument('--bdt_rw', action='store_true')

    args = parser.parse_args()
    run_name = args.run_name
    specs = args.specs

    fit_strat = args.fit_strat
    gc_onflag = args.gc_onflag
    data_build = args.data_build
    # data_fit = args.data_fit
    data_plot = args.data_plot

    data_corr_build = args.data_corr_build
    data_corr_plot = args.data_corr_plot

    data_post_d_plot = args.data_post_d_plot

    split_flag = args.split_flag
    fix_flag = args.fix_flag
    smear_flag = args.smear_flag

    norm_build = args.norm_build
    norm_plot = args.norm_plot

    mc_build = args.mc_build
    mc_fit = args.mc_fit
    mc_plot = args.mc_plot
    mc_data_comp = args.mc_data_comp
    mc_res = args.mc_res

    mc_build_pizy = args.mc_build_pizy

    bdt_rw = args.bdt_rw

    mc_event_shape_list = [
        # # spectrum_id, [List of MC files], [List of shapes to fit], "B Mean Guess", "B_Mean_Guess_Window", "B Fit Window"
        # ["Z_m_p_01", ["01_Z_m_p_11198006"], ["G","DG","BG","GEP","BGEP"], 5280, 10, 25],
        # # ["Z_m_p_02", ["02_Z_m_p_11198400"], ["G","DG","BG","GEP","BGEP"], 5130, 10, 50],
        # # ["Z_m_p_03", ["02_Z_m_p_11198400"], ["G","DG","BG","GEP","BGEP"], 5130, 10, 50],
        # ["Z_m_p_04", ["04_Z_m_p_11198401"], ["G","DG","BG","GEP","BGEP"], 4985, 20, 50],
        # ["Z_m_p_0203", ["02_Z_m_p_11198400"], ["G","DG","BG","GEP","BGEP"], 5130, 10, 50],

            # ["Z_z_z_04",["04_Z_z_z_11198023"], ["G","DG","BG","GEP","BGEP"], 4970, 10, 50],
            # # ["Z_z_z_07",["07_Z_z_z_12197045"], ["G","DG","BG","GEP","BGEP"], 5125, 15, 50],
            # ["Z_z_z_08",["08_Z_z_z_12197423"], ["G","DG","BG","GEP","BGEP"], 4980, 10, 65],
            # # ["Z_z_z_09",["09_Z_z_z_11196019"], ["G","DG","BG","GEP","BGEP"], 5280, 10, 30],
            # # ["Z_z_z_10",["10_Z_z_z_11196413"], ["GEPAddBG"], 5130, 30, 60],
            # ["Z_z_z_12",["12_Z_z_z_11196414"], ["G","DG","BG","GEP","BGEP"], 4980, 10, 80],
        # #
        # ["Z_z_z_0710", ["07_Z_z_z_12197045", "10_Z_z_z_11196413"], ["GAddBGEP"],  5125, 30, 50],
        # # ["Z_z_z_040812", ["04_Z_z_z_11198023", "08_Z_z_z_12197423", "12_Z_z_z_11196414"], ["G","DG","BG","GEP","BGEP"],  4980, 10, 80],
        # #
        # # ["P_z_p_02", ["02_P_z_p_11198005"], ["G","DG","BG","GEP","BGEP"], 5130, 10, 60],
        # # ["P_z_p_04", ["04_P_z_p_11198410"], ["G","DG","BG","GEP","BGEP"], 4985, 10, 60],
        # # ["P_z_p_05", ["05_P_z_p_12197023"], ["G","DG","BG","GEP","BGEP"], 5280, 10, 40],
        # # ["P_z_p_06", ["06_P_z_p_12197410"], ["GEPAddBGEP"], 5130, 15, 80],
        # # ["P_z_p_07", ["07_P_z_p_12197400"], ["G","DG","BG","GEP","BGEP"], 5130, 15, 60],
        # # ["P_z_p_08", ["08_P_z_p_12197401"], ["G","DG","BG","GEP","BGEP"], 4985, 20, 70],
        # # # # #
        # ["P_z_p_020607", ["02_P_z_p_11198005", "06_P_z_p_12197410", "07_P_z_p_12197400"], ["GAddBGEP"], 5130, 20, 50],
        # ["P_z_p_0408", ["04_P_z_p_11198410", "08_P_z_p_12197401"], ["GAddBGEP"], 4980, 30, 70],
        # #
        # ["M_m_z_03", ["02_P_z_p_11198005"], ["G","DG","BG","GEP","BGEP"], 5130, 10, 60],
        # ["M_m_z_04", ["04_P_z_p_11198410"], ["G","DG","BG","GEP","BGEP"], 4985, 10, 60],
        # #
        # ["P_z_pst_07", ["07_P_z_pst_12197045"], ["G","DG","BG","GEP","BGEP"], 5280, 10, 30],
        # ["P_z_pst_04", ["04_P_z_pst_11198023"], ["G","DG","BG","GEP","BGEP"], 5130, 30, 50],
        # ["P_z_pst_08", ["08_P_z_pst_12197423"], ["G","DG","BG","GEP","BGEP"], 5130, 30, 50],
        # ["P_z_pst_0408", ["04_P_z_pst_11198023", "08_P_z_pst_12197423"], ["GAddBGEP"], 5130, 30, 50],

        ["Zs_sm_p_13", ["13_Zs_sm_p_13198040"], ["G","DG","BG","GEP","BGEP"], 5367, 30, 25],
        ["Zs_sm_p_14", ["14_Zs_sm_p_13198200"], ["G","DG","BG","GEP","BGEP"], 5220, 30, 100],
        ["Zs_sm_p_15", ["15_Zs_sm_p_13198400"], ["G","DG","BG","GEP","BGEP"], 5220, 30, 50],
        ["Zs_sm_p_16", ["16_Zs_sm_p_13198600"], ["G","DG","BG","GEP","BGEP"], 5075, 20, 100],

        # ["norm7", ["norm7_norm7_12197008"], ["G","DG","BG","GEP","BGEP"], 5280, 30, 50],
        # ["norm8", ["norm8_norm8_11198007"], ["G","DG","BG","GEP","BGEP"], 5280, 30, 50],
    ]

    # Build no correlation fit, no splits
    if data_build:
        for spec in specs:
            data_fit_v1.build_nn_sw_fit(run_name, spec, split_flag = False, fix_flag = False, smear_flag = False)
    if data_plot:
        for spec in specs:
            data_fit_v1.plot_nn_data(run_name, spec, split_flag)


    if data_corr_build:
        data_fit_v1.build_comp_ws_fit(run_name, ["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"], split_flag = True, fix_flag = True, smear_flag = True, gc_onflag = 0)
    if data_corr_plot:
        data_fit_v1.plot_coor_data(run_name, specs)
        #"Z_z_z","P_z_p","P_z_pst
    if mc_build:
        for nms in mc_event_shape_list:
            # el = " ".join(str(x) for x in nms[1])
            # sl = " ".join(str(x) for x in nms[2])
            name = nms[0]
            file_list = nms[1]
            shape_list = nms[2]
            mean_start = nms[3]
            mean_window = nms[4]
            fit_window = nms[5]
            mc_ws_build_script.build_mc_ws(name, file_list, shape_list, mean_start, mean_window, fit_window)
    if mc_fit:
        for nms in mc_event_shape_list:
            shape_list = nms[2]
            mc_ws_build_script.run_mc_fit(nms[0], shape_list)
    if mc_plot:
        for nms in mc_event_shape_list:
            shape_list = nms[2]
            mc_ws_build_script.plot_mc(nms[0], shape_list)
    if mc_build_pizy:
        spec = "Z_z_z_0710"
        piz_spec_list = ["07_Z_z_z_12197045"]
        y_spec_list = ["10_Z_z_z_11196413"]
        mc_ws_build_script.build_mc_pizy_ws(spec, piz_spec_list, y_spec_list, 5125, 20, 65)
        spec = "Z_z_z_040812"
        piz_spec_list = ["04_Z_z_z_11198023"]
        y_spec_list = ["08_Z_z_z_12197423", "12_Z_z_z_11196414"]
        mc_ws_build_script.build_mc_pizy_ws(spec, piz_spec_list, y_spec_list, 4980, 30, 70)
        spec = "P_z_p_0408"
        piz_spec_list = ["04_P_z_p_11198410"]
        y_spec_list = ["08_P_z_p_12197401"]
        mc_ws_build_script.build_mc_pizy_ws(spec, piz_spec_list, y_spec_list, 4980, 30, 70)
        spec = "P_z_p_020607"
        piz_spec_list = ["02_P_z_p_11198005", "07_P_z_p_12197400"]
        y_spec_list = ["06_P_z_p_12197410"]
        mc_ws_build_script.build_mc_pizy_ws(spec, piz_spec_list, y_spec_list, 5130, 30, 70)
    # if data_post_d_plot:
    #     for spec in ["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]:
    #         data_fit_v1.plot_post_d_data(spec)
    # if norm_build:
    #     for spec in ["norm7","norm8"]:
    #         data_fit_v1.build_norm_fit(run_name, spec, split_flag, fix_flag, smear_flag)
    # if norm_plot:
    #     for spec in ["norm7","norm8"]:
    #         data_fit_v1.plot_norm(run_name, spec)
    # if mc_res:
    #     mc_res_test.build_tmc_ws([("02_Z_m_p_11198400", 1), ("04_Z_m_p_11198401", 2)])
    #     # for mc_id, mc_file, mp in zip(["02","04"],["02_Z_m_p_11198400","04_Z_m_p_11198401"],[1,2]):
    #     #     mc_res_test.create_new_tree("Z_m_p", "01", mc_id, mc_file, mp)
    #     #     mc_res_test.compare_res("Z_m_p", mc_id, mc_file, mp)
    #     # mc_res_test.build_tmc_ws([("07_Z_z_z_12197045", 1),
    #     #                           ("10_Z_z_z_11196413", 1),
    #     #                           ("04_Z_z_z_11198023", 2),
    #     #                           ("08_Z_z_z_12197423", 2),
    #     #                           ("12_Z_z_z_11196414", 2)])
    #     # for mc_id, mc_file, mp in zip(["07","10","04","08","12"],["07_Z_z_z_12197045","10_Z_z_z_11196413","04_Z_z_z_11198023", "08_Z_z_z_12197423","12_Z_z_z_11196414"],[1,1,2,2,2]):
    #     #     mc_res_test.create_new_tree("Z_z_z", "09", mc_id, mc_file, mp)
    #     #     mc_res_test.compare_res("Z_z_z", mc_id, mc_file, mp)
    #     # mc_res_test.build_tmc_ws([("02_P_z_p_11198005", 1),
    #     #                           ("06_P_z_p_12197410", 1),
    #     #                           ("07_P_z_p_12197400", 1),
    #     #                           ("04_P_z_p_11198410", 2),
    #     #                           ("08_P_z_p_12197401", 2)])
    #     # for mc_id, mc_file, mp in zip(["02","06","07","04","08"],["02_P_z_p_11198005","06_P_z_p_12197410","07_P_z_p_12197400","04_P_z_p_11198410","08_P_z_p_12197401"],[1,1,1,2,2]):
    #     #     mc_res_test.create_new_tree("P_z_p", "05", mc_id, mc_file, mp)
    #     #     mc_res_test.compare_res("P_z_p", mc_id, mc_file, mp)
    #     # mc_res_test.build_tmc_ws([("04_P_z_pst_11198023", 1),
    #     # #                           ("08_P_z_pst_12197423", 1)])
    #     # for mc_id, mc_file, mp in zip(["04","08"],["04_P_z_pst_11198023", "08_P_z_pst_12197423"],[1,1]):
    #     #     mc_res_test.create_new_tree("P_z_pst", "07", mc_id, mc_file, mp)
    #     #     mc_res_test.compare_res("P_z_pst", mc_id, mc_file, mp)
    if bdt_rw:
        # ("Z_m_p", "01")]
        #[("Z_m_p", "02"), ("Z_m_p", "04")]
        #("Z_z_z", "09"),
        #("Z_z_z", "07"), ("Z_z_z", "10")
        #("Z_z_z", "04"), ("Z_z_z", "08"), ("Z_z_z", "12")
        # mc_list = [("Z_m_p", "01", [20e6, 12e6, 12e6]), ("Z_m_p", "02", [18.5e6, 11e6, 11e6]), ("Z_m_p", "04", [17e6, 10e6, 10e6])]
        mc_list = [("Z_z_z", "09", [50e6, 50e6, 50e6]), ("Z_z_z", "07", [50e6, 50e6, 50e6]), ("Z_z_z", "10", [50e6, 50e6, 50e6])]
        # mc_list = [("Z_z_z", "04", [50e6, 50e6, 50e6]), ("Z_z_z", "08", [50e6, 50e6, 50e6]), ("Z_z_z", "12", [50e6, 50e6, 50e6])]
        # mc_list = [("Z_m_p", "01", [20e6, 12e6, 12e6]), ("Z_m_p", "02", [18.5e6, 11e6, 11e6]), ("Z_m_p", "04", [17e6, 10e6, 10e6])]
        # mc_list = [("Z_m_p", "01", [20e6, 12e6, 12e6]), ("Z_m_p", "02", [18.5e6, 11e6, 11e6]), ("Z_m_p", "04", [17e6, 10e6, 10e6])]
        # mc_list = [("Z_m_p", "01", [20e6, 12e6, 12e6]), ("Z_m_p", "02", [18.5e6, 11e6, 11e6]), ("Z_m_p", "04", [17e6, 10e6, 10e6])]
        # mc_list = [("Z_m_p", "01", [20e6, 12e6, 12e6]), ("Z_m_p", "02", [18.5e6, 11e6, 11e6]), ("Z_m_p", "04", [17e6, 10e6, 10e6])]


        bdt_test.bdt_test(run_name, mc_list, "all")
