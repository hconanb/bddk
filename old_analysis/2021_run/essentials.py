import ROOT, os, sys
import datetime

from uncertainties import ufloat, covariance_matrix,correlation_matrix
from uncertainties import unumpy

import glob
import xlrd

from pandas import *
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

sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run')
sys.path.append('/home/hbernste/lhcb-analysis/analysisbase/python')
ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis/rootclasses/lib/librootclasses.so")

def arg_as_list(s):
    v = " ".join(str(x) for x in s)
    return v

#fit strategies, bkg, zero missing charm peak, 1 missing charm peak, 2 missing charm peak
name_data_ws = "data_all_just_z"
name_fit_ws = "fit_data_all_just_z"

#determines if fit uses gaussian constratits for branching fractions and efficiecies
gc_onflag = 0

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
data_basepath = "/mnt/c/Users/Harris/Desktop/rootfiles/data_run2/"
analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run"

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

# def old_grab_file_list(type, event_list):
#     file_list = []
#     for event in event_list:
#         file_list = file_list + glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}/{event}_*.root")
#     return file_list

def get_shapes_bkg(spec, flag, ws):
    if flag == "Exponential":
        ws.factory(f"Exponential:{spec}_spectrum_bkg(B_DTF_M, c0_{spec}[0, -2, 2])")
    if flag == "Chebychev":
        ws.factory(
            f"Chebychev:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[0.,-3,3],c1_{spec}[0.,-3,3]}})"
        )
    if flag == "Bernstein":
        ws.factory(
            f"Bernstein:{spec}_spectrum_bkg(B_DTF_M,{{c0_{spec}[1,0,10], c1_{spec}[1,0,10], c2_{spec}[1,0,10], c3_{spec}[1,0,10]}})"
        )

def get_free_shapes(ws, spec, fit_strat, s_list):

    #slist = [Harris ID, "Shape", "B Mean Guess", "B_Mean_Guess_Window", "B Fit Window for MC"]
    if spec != s_list[0]:
        spec = f"{spec}_{s_list[0]}"
    shape_flag = s_list[1]
    mean_start = s_list[2]
    window = s_list[3]

    print(spec)
    b_low = mean_start - window
    b_high = mean_start + window

    if shape_flag == "G":
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_fit(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}], width_{spec}_{shape_flag}[5.0,0.0,20.0])"
        )
    if shape_flag == "DG":
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_a(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}], width_{spec}_{shape_flag}[12, 0.0, 50.0])"
        )
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_b(B_DTF_M, mean_{spec}_{shape_flag}, width_b_{spec}_{shape_flag}[21.0, 0.0, 50.0])"
        )
        ws.factory(f"SUM::{spec}_{shape_flag}_fit({spec}_{shape_flag}_a_frac[0.8,0,1]*{spec}_{shape_flag}_a, {spec}_{shape_flag}_b)")
    if shape_flag == "GEP":
        ws.factory(
            f"RooGaussExp::{spec}_{shape_flag}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,4.0,30.0],alpha_{spec}_{shape_flag}[3,0.05,7.0])"
        )
    if shape_flag == "BGEP":
        ws.factory(
            f"RooBifurGaussExp::{spec}_{shape_flag}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_L_{spec}_{shape_flag}[10,1.0,30.0],width_R_{spec}_{shape_flag}[5,1.0,30.0],alpha_1_{spec}_{shape_flag}[2,0.00,10.0],alpha_2_{spec}_{shape_flag}[8,0.00,10.0])"
        )
    if shape_flag == "BG":
        ws.factory(
            f"BifurGauss::{spec}_{shape_flag}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_1_{spec}_{shape_flag}[20,0.00,50.0],width_2_{spec}_{shape_flag}[5,0.00,50.0])"
        )
    if shape_flag == "GAddBGEP":
        ws.factory(
            f"BifurGauss::{spec}_{shape_flag}_fit_a(B_DTF_M, mean_{spec}_{shape_flag}_a[{mean_start},{b_low},{b_high}],width_1_{spec}_{shape_flag}_a[20,0.00,50.0],width_2_{spec}_{shape_flag}_a[5,0.00,50.0])"
        )
        # ws.factory(
        #     f"RooBifurGaussExp:{spec}_{shape_flag}_fit_a(B_DTF_M, mean_{spec}_{shape_flag}_a[{mean_start},{b_low},{b_high}],width_L_{spec}_{shape_flag}_a[10,1.0,30.0],width_R_{spec}_{shape_flag}_a[20,4.0,30.0],alpha_1_{spec}_{shape_flag}_a[2,0.05,7.0],alpha_2_{spec}_{shape_flag}_a[3,0.05,7.0])"
        # )
        ws.factory(
            f"Gaussian::{spec}_{shape_flag}_fit_b(B_DTF_M, mean_{spec}_{shape_flag}_b[{mean_start},{b_low},{b_high}], width_{spec}_{shape_flag}_b[5.0,0.0,20.0])"
        )
        # ws.factory(
        #     f"RooGaussExp:{spec}_{shape_flag}_fit_b(B_DTF_M, mean_{spec}_{shape_flag}_a,width_{spec}_{shape_flag}_a[10,1.0,30.0],alpha_{spec}_{shape_flag}_a[3,0.05,7.0])"
        # )
        # ws.factory(
        #     f"RooGaussExp:{spec}_{shape_flag}_fit_b(B_DTF_M, mean_{spec}_{shape_flag}_a,width_{spec}_{shape_flag}_b[10,4.0,30.0],alpha_{spec}_{shape_flag}_b[3,0.05,7.0])"
        # )

        ws.factory(f"SUM:{spec}_{shape_flag}_fit({spec}_{shape_flag}_a_frac[0.647]*{spec}_{shape_flag}_fit_a, {spec}_{shape_flag}_fit_b)")

    # if shape_flag == "cb1R":
    #     ws.factory(
    #         f"CBShape::{spec}_{shape_flag}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,0,50],alpha_{spec}_{shape_flag}[-5,-10,-0.0], n_{spec}_{shape_flag}[5,0,10])"
    #     )
    # if shape_flag == "cb1L":
    #     ws.factory(
    #         f"CBShape::{spec}_{shape_flag}_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,0,50], alpha_{spec}_{shape_flag}[5,0.0,10.0], n_{spec}_{shape_flag}[5,0,10])"
    #     )
    # if shape_flag == "cb2":
    #     ws.factory(
    #         f"CBShape::{spec}_a_fit(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{b_low},{b_high}],width_{spec}_{shape_flag}[10,0,100],alpha_l_{spec}_{shape_flag}[-1,-20,-0.0],  n1_{spec}_{shape_flag}[20,0,1000])"
    #     )
    #     ws.factory(
    #         f"CBShape::{spec}_b_fit(B_DTF_M,mean_{spec}_{shape_flag},                                                  width_{spec}_{shape_flag},         alpha_r_{spec}_{shape_flag}[1,0.0,20.0],   n2_{spec}_{shape_flag}[20,0,1000])"
    #     )
    #     ws.factory(
    #         f"SUM::{spec}_{shape_flag}_fit({spec}_a_frac[0.5,0,1]*{spec}_a_fit, {spec}_b_fit)"
    #     )
    if spec != s_list[0]:
        temp_pdf = ws.pdf(f"{spec}_{shape_flag}_fit")
        temp_pdf.SetName(f"{spec}_fit")
        temp_pdf.Print()
