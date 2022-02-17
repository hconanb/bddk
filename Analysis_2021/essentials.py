import ROOT, os, sys
import datetime

from uncertainties import ufloat, covariance_matrix,correlation_matrix
from uncertainties import unumpy

import glob
import xlrd
import shutil
import pprint

from pandas import *
from pandas.api.types import CategoricalDtype

import numpy as np
import math as math
from array import array
from rootutils import *

import root_pandas as rp
import itertools
RDF = ROOT.ROOT.RDataFrame

import ast
import argparse
from argparse import RawTextHelpFormatter

ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")
ROOT.gStyle.SetPalette(ROOT.kBird)


def arg_as_list(s):
    v = " ".join(str(x) for x in s)
    return v

a_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"

#Constants
dwindow = 25.0
dmaxbkg = 85.0
dminbkg = 35.0
B0mass = 5279.61
Bpmass = 5279.29
Bsmass = 5366.79

dstpmass = 2010.27
pi0mass = 134.9766

bbins = 100
bmin = 4800
bmax = 5600

p_bbins = 100

p0_n8_bmin = 5230
p0_n8_bmax = 5330

p0_bmin = 5250
p0_bmax = 5310

p1_bmin = 5080
p1_bmax = 5180

p2_bmin = 4935
p2_bmax = 5035

p0s_bmin = 5340
p0s_bmax = 5390

dpmass = 1869.62
d0mass = 1864.84
dsmass = 1968.47

delta_st_m_cut = 200
#
root_basepath = "/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/"
# data_basepath = "/mnt/c/Users/Harris/Desktop/rootfiles/data_run2/"
analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"

##### MC ID's #####
B0_ID = 511
Bp_ID = 521
Bs_ID = 531

Dp_ID = 411
D0_ID = 421
Ds_ID = 431

Dpst_ID = 413
D0st_ID = 423
Dsst_ID = 433

k_ID = 321
pi_ID = 211
kst0_ID = 313

def DrawStack(
    pad,
    hlist,
    legend=(0.3, 0.21, 0.3, 0.21),
    XaxisTitle=True,
    YaxisTitle=True,
    drawopts="nostack plc pmc",
    legfillstyle=0,
    legfillcoloralpha=None,
):
    """Draw hlist on pad, return THStack.

    arguments:
    pad -- TPad on which to draw histograms
    hlist -- list of histograms to draw
    legend -- if bool(legend), draws with given coordinates using pad.BuildLegend
    XaxisTitle -- if True, sets title of Xaxis to that of hlist[0]
        if str, sets to given value
        if not bool(XaxisTitle), no title is added
    YaxisTitle -- like XaxisTitle, but for Yaxis
    drawopts -- passed to THStack.Draw
    legfillstyle -- passed to SetFillStyle; 0 = transparent
    legfillcoloralpha -- passed to SetFillColorAlpha
    """

    pad.cd()
    hs = ROOT.THStack()
    for h in hlist:
        hs.Add(h)
    hs.Draw(drawopts)
    if XaxisTitle:
        if isinstance(XaxisTitle, str):
            hs.GetXaxis().SetTitle(XaxisTitle)
        else:
            hs.GetXaxis().SetTitle(hlist[0].GetXaxis().GetTitle())
    if YaxisTitle:
        if isinstance(YaxisTitle, str):
            hs.GetYaxis().SetTitle(YaxisTitle)
        else:
            hs.GetYaxis().SetTitle(hlist[0].GetYaxis().GetTitle())
    if legend:
        leg = pad.BuildLegend(*legend)
        leg.SetFillStyle(legfillstyle)
        if legfillcoloralpha:
            leg.SetFillColorAlpha(*legfillcoloralpha)
    pad.Update()
    return hs


def save_pdf(thing, folder, name, rpflag):
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}/{folder}/'):
        os.makedirs(f'plots/{now.month}_{now.day}/{folder}/')
    thing.SaveAs(f"plots/{now.month}_{now.day}/{folder}/{name}.pdf")
    print(f"saved: {name}.pdf in plots/{now.month}_{now.day}/{folder}/")

def save_png(thing, folder, name, rpflag = 0):
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}/{folder}/'):
        os.makedirs(f'plots/{now.month}_{now.day}/{folder}/')
    if rpflag == 1:
        thing.save(f"plots/{now.month}_{now.day}/{folder}/{name}.png")
    else:
        thing.SaveAs(f"plots/{now.month}_{now.day}/{folder}/{name}.png")
    print(f"saved: {name}.png in plots/{now.month}_{now.day}/{folder}/")


def get_dwindow_values(spec, d1_flag, d2_flag, dst_flag = False, rflag = "apply"):
    if rflag == "print":
        d1window_file = ROOT.TFile(f"../build_rootfiles/d_window_root_files/d_{d1_flag}_mass_fits.root", "READ")
        d1window_ws = d1window_file.Get(f"d_{d1_flag}_mass_fits")

        d1_mstart = d1window_ws.var(f"mean").getValV()
        d1_a_std = d1window_ws.var(f"width_a").getValV()
        d1_b_std = d1window_ws.var(f"width_b").getValV()

        if (d1_a_std > d1_b_std):
            d1_std = d1_a_std
        else:
            d1_std = d1_b_std
        return (d1_mstart, d1_std)

    if rflag == "apply":
        d1window_file = ROOT.TFile(f"d_window_root_files/d_{d1_flag}_mass_fits.root", "READ")
        d1window_ws = d1window_file.Get(f"d_{d1_flag}_mass_fits")
        d2window_file = ROOT.TFile(f"d_window_root_files/d_{d2_flag}_mass_fits.root", "READ")
        d2window_ws = d2window_file.Get(f"d_{d2_flag}_mass_fits")

        d1_mstart = d1window_ws.var(f"mean").getValV()
        d1_a_std = d1window_ws.var(f"width_a").getValV()
        d1_b_std = d1window_ws.var(f"width_b").getValV()

        # d1_y = d1window_ws.var(f"n_D_signal").getValV()
        d1_y_frac = d1window_ws.var(f"D1_a_frac").getValV()
        # d1_a = d1_y*d1_y_frac
        # d1_b = d1_y*(1 - d1_y_frac)

        d2_mstart = d2window_ws.var(f"mean").getValV()
        d2_a_std = d2window_ws.var(f"width_a").getValV()
        d2_b_std = d2window_ws.var(f"width_b").getValV()

        # d2_y = d2window_ws.var(f"n_D_signal").getValV()
        d2_y_frac = d2window_ws.var(f"D1_a_frac").getValV()
        # d2_a = d2_y*d2_y_frac
        # d2_b = d2_y*(1 - d2_y_frac)

        # d1_std = d1_a_std*d1_y_frac + d1_b_std*(1-d1_y_frac)
        # d2_std = d2_a_std*d2_y_frac + d2_b_std*(1-d2_y_frac)

        if (d1_a_std > d1_b_std):
            d1_std = d1_a_std
        else:
            d1_std = d1_b_std

        if (d2_a_std > d2_b_std):
            d2_std = d2_a_std
        else:
            d2_std = d2_b_std
        # if d1_flag == "mp" and d2_flag == "d0k3pi":
        #     d1_mstart = d1window_ws.var(f"mean").getValV()
        #     d1_a_std = d1window_ws.var(f"width_a").getValV()
        #
        #     d2_mstart = d2window_ws.var(f"mean").getValV()
        #     d2_a_std = d2window_ws.var(f"width_a").getValV()
        #     d2_b_std = d2window_ws.var(f"width_b").getValV()
        #     if (d2_a_std > d2_b_std):
        #         d2_std = d2_a_std
        #     else:
        #         d2_std = d2_b_std
        #     d1_std = d1_a_std

        if dst_flag:
            dstwindow_file = ROOT.TFile(f"d_window_root_files/d_dst_mass_fits.root", "READ")
            dstwindow_ws = dstwindow_file.Get(f"d_dst_mass_fits")
            dst_mstart = dstwindow_ws.var(f"mean").getValV()
            dst_a_std = dstwindow_ws.var(f"width_a").getValV()
            dst_b_std = dstwindow_ws.var(f"width_b").getValV()
            if (dst_a_std > dst_b_std):
                dst_std = dst_a_std
            else:
                dst_std = dst_b_std
            dstwindow = dst_std * 2

        d1window = d1_std * 2
        d2window = d2_std * 2

        if spec == "P_z_pst" or spec == "Z_m_pst":
            mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window} && abs(D2stmD_M - {dst_mstart}) < {dstwindow}) "
        else:
            mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window})"

        return mass_cut



def get_shapes_bkg(spec, flag, ws):
    if flag == "Exponential":
        ws.factory(f"Exponential:{spec}_spectrum_bkg(B_DTF_M, c0_{spec}[0, -0.01, 0.01])")
    if flag == "Chebychev":
        ws.factory(
            f"Chebychev:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[0.,-3,3],c1_{spec}[0.,-3,3]}})"
        )
    if flag == "Bernstein":
        ws.factory(
            f"Bernstein:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[1,0,10], c1_{spec}[1,0,10], c2_{spec}[1,0,10], c3_{spec}[1,0,10]}})"
        )

def get_free_shapes(ws, spec, fit_strat, shape_flag, mean_start, mean_window):

    #slist = [Harris ID, "Shape", "B Mean Guess", "B_Mean_Guess_Window", "B Fit Window for MC"]
    # if spec != s_list[0]:
    #     spec = f"{spec}_{s_list[0]}"

    b_low = mean_start - mean_window
    b_high = mean_start + mean_window

    if shape_flag == "G":
        ws.factory(
            f"Gaussian::{spec}_fit(B_DTF_M, mean_{spec}[{mean_start},{b_low},{b_high}], width_{spec}[5.0,0.1,30.0])"
        )
    if shape_flag == "DG":
        ws.factory(
            f"Gaussian::{spec}_a(B_DTF_M, mean_{spec}[{mean_start},{b_low},{b_high}], width_a_{spec}[10, 0.1, 40.0])"
        )
        ws.factory(
            f"Gaussian::{spec}_b(B_DTF_M, mean_{spec}, width_b_{spec}[20.0, 0.1, 40.0])"
        )
        ws.factory(f"SUM::{spec}_fit({spec}_a_frac[0.8,0,1]*{spec}_a, {spec}_b)")
    if shape_flag == "GEP":
        ws.factory(
            f"RooGaussExp::{spec}_fit(B_DTF_M,mean_{spec}[{mean_start},{b_low},{b_high}],width_{spec}[10,4.0,30.0],alpha_{spec}[3,0.05,7.0])"
        )
    if shape_flag == "BGEP":
        ws.factory(
            f"RooBifurGaussExp::{spec}_fit(B_DTF_M,mean_{spec}[{mean_start},{b_low},{b_high}],width_L_{spec}[10,1.0,30.0],width_R_{spec}[5,1.0,30.0],alpha_1_{spec}[2,0.00,10.0],alpha_2_{spec}[8,0.00,10.0])"
        )
    if shape_flag == "BG":
        ws.factory(
            f"BifurGauss::{spec}_fit(B_DTF_M,mean_{spec}[{mean_start},{b_low},{b_high}],width_1_{spec}[20,0.01,50.0],width_2_{spec}[5,0.01,50.0])"
        )
    if shape_flag == "GAddBG":
        ws.factory(
            f"BifurGauss::{spec}_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}],width_1_{spec}_a[20,0.01,50.0],width_2_{spec}_a[5,0.01,50.0])"
        )
        ws.factory(
            f"Gaussian:{spec}_fit_b(B_DTF_M, mean_{spec}_a, width_{spec}_b[10,1.0,30.0])"
        )
        ws.factory(f"SUM:{spec}_fit({spec}_a_frac[0.647,0.5,0.75]*{spec}_fit_a, {spec}_fit_b)")
    if shape_flag == "GEPAddBG":
        ws.factory(
            f"BifurGauss::{spec}_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}],width_1_{spec}_a[20,0.01,50.0],width_2_{spec}_a[5,0.01,50.0])"
        )
        ws.factory(
            f"RooGaussExp:{spec}_fit_b(B_DTF_M, mean_{spec}_a, width_{spec}_b[10,1.0,30.0],alpha_{spec}_b[3,0.05,7.0])"
        )
        ws.factory(f"SUM:{spec}_fit({spec}_a_frac[0.647,0.5,0.75]*{spec}_fit_a, {spec}_fit_b)")
    if shape_flag == "GAddBGEP":
        ws.factory(
            f"RooBifurGaussExp:{spec}_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}],width_L_{spec}_a[10,0.01,30.0],width_R_{spec}_a[20,0.01,30.0],alpha_1_{spec}_a[2,0.001,10.0],alpha_2_{spec}_a[3,0.001,10.0])"
        )
        ws.factory(
            f"Gaussian::{spec}_fit_b(B_DTF_M, mean_{spec}_b[{mean_start},{b_low},{b_high}], width_{spec}_b[5.0,0.1,30.0])"
        )
        ws.factory(f"SUM:{spec}_fit({spec}_a_frac[0.647,0.1,0.9]*{spec}_fit_a, {spec}_fit_b)")
    if shape_flag == "GEPAddBGEP":
        ws.factory(
            f"RooBifurGaussExp:{spec}_fit_a(B_DTF_M, mean_{spec}_a[{mean_start},{b_low},{b_high}],width_L_{spec}_a[10,1.0,30.0],width_R_{spec}_a[20,4.0,30.0],alpha_1_{spec}_a[2,0.05,7.0],alpha_2_{spec}_a[3,0.05,7.0])"
        )
        ws.factory(
            f"RooGaussExp:{spec}_fit_b(B_DTF_M, mean_{spec}_a,width_{spec}_b[10,1.0,30.0],alpha_{spec}_b[3,0.05,7.0])"
        )
        ws.factory(f"SUM:{spec}_fit({spec}_a_frac[0.647,0.5,0.75]*{spec}_fit_a, {spec}_fit_b)")




    # if shape_flag == "cb1R":
    #     ws.factory(
    #         f"CBShape::{spec}_fit(B_DTF_M,mean_{spec}[{mean_start},{b_low},{b_high}],width_{spec}[10,0,50],alpha_{spec}[-5,-10,-0.0], n_{spec}[5,0,10])"
    #     )
    # if shape_flag == "cb1L":
    #     ws.factory(
    #         f"CBShape::{spec}_fit(B_DTF_M,mean_{spec}[{mean_start},{b_low},{b_high}],width_{spec}[10,0,50], alpha_{spec}[5,0.0,10.0], n_{spec}[5,0,10])"
    #     )
    # if shape_flag == "cb2":
    #     ws.factory(
    #         f"CBShape::{spec}_a_fit(B_DTF_M,mean_{spec}[{mean_start},{b_low},{b_high}],width_{spec}[10,0,100],alpha_l_{spec}[-1,-20,-0.0],  n1_{spec}[20,0,1000])"
    #     )
    #     ws.factory(
    #         f"CBShape::{spec}_b_fit(B_DTF_M,mean_{spec},                                                  width_{spec},         alpha_r_{spec}[1,0.0,20.0],   n2_{spec}[20,0,1000])"
    #     )
    #     ws.factory(
    #     )
    #         f"SUM::{spec}_fit({spec}_a_frac[0.5,0,1]*{spec}_a_fit, {spec}_b_fit)"
# mc_spec_list = [
#     # "01_Z_m_p_11198006",
#     # "02_Z_m_p_11198400",
#     # "02_P_z_p_11198005",
#     # "04_Z_m_p_11198401",
#     # "04_P_z_p_11198410",
#     "04_Z_z_z_11198023",
#     # # "04_P_z_pst_11198023",
#     # # "05_P_z_p_12197023",
#     # # "06_P_z_p_12197410",
#     # # "07_P_z_p_12197400",
#     "07_Z_z_z_12197045",
#     # # "07_P_z_pst_12197045",
#     # # "08_P_z_p_12197401",
#     "08_Z_z_z_12197423",
#     # "08_P_z_pst_12197423",
#     # "09_Z_z_z_11196019",
#     # "10_Z_z_z_11196413",
#     # "12_Z_z_z_11196414",
#     # "13_Zs_sm_p_13198040",
#     # "14_Zs_sm_p_13198200",
#     # "15_Zs_sm_p_13198400",
#     # "16_Zs_sm_p_13198600",
#     # "norm7_norm7_12197008",
#     # "norm8_norm8_11198007",
# ]

# def make_folders():
#     for mc_name in mc_spec_list:
#         for year in ["2016", "2017", "2018"]:
#             new_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_mc/{mc_name}/{year}"
#             if os.path.isdir(new_path) == False:
#                 os.makedirs(new_path)
# #
# def rename_folders():
#
#     df = pandas.read_excel(
#         "/mnt/c/Users/Harris/Desktop/jobs2.xlsx"
#         )
#     print(df)
#     for name in mc_spec_list:
#         id = name.rsplit("_")[-1]
#         print(id)
#         if "P_z_pst" not in id:
#             df_6 = df.loc[df['SPEC'] == int(id)]
#             print(df_6)
#             jobs_n = df_6["JOB"].values
#             years_n = df_6["YEAR"].values
#             for job_name, year_name in zip(jobs_n, years_n):
#
#                 old_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_mc/raw/{job_name}"
#                 new_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_base_mc/{name}/{year_name}/"
#
#                 if os.path.isdir(old_path) == True:
#
#                     print(f"move {old_path} to {new_path}")
#                     shutil.move(old_path, new_path)
