import sys
import os
import pickle

import ROOT
import uproot
import pandas as pd
import argparse

import numpy as np
import pandas as pd

import argparse
import glob as glob

import math
import itertools

from uncertainties import ufloat, covariance_matrix, correlation_matrix
from uncertainties import unumpy

ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")

b0mass = 5279.61
bpmass = 5279.29
bsmass = 5366.79
dpmass = 1869.62
d0mass = 1864.84
dsmass = 1968.47
dstpmass = 2010.27
pi0mass = 134.9766

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

ntuple_col_list = [
"B_DTF_M",
"B_M",

]

id_to_spec_dict = {
     "Z_m_p":           "Z_m_p",
     "Z_z_z":           "Z_z_z",
     "P_z_p":           "P_z_p",
     "M_m_z":           "M_m_z",
     "P_z_pst":         "P_z_pst",
     "norm7":           "norm7",
     "norm8":           "norm8",
     "01_Z_m_p":        "01_Z_m_p_11198006" ,
     "02_Z_m_p":        "02_Z_m_p_11198400" ,
     "02_P_z_p":        "02_P_z_p_11198005" ,
     "04_Z_m_p":        "04_Z_m_p_11198401" ,
     "04_P_z_p":        "04_P_z_p_11198410" ,
     "04_Z_z_z":        "04_Z_z_z_11198023" ,
     "04_P_z_pst":      "04_P_z_pst_11198023" ,
     "05_P_z_p":        "05_P_z_p_12197023" ,
     "06_P_z_p":        "06_P_z_p_12197410",
     "07_P_z_p":        "07_P_z_p_12197400",
     "07_Z_z_z":        "07_Z_z_z_12197045",
     "07_P_z_pst":      "07_P_z_pst_12197045",
     "08_P_z_p":        "08_P_z_p_12197401",
     "08_Z_z_z":        "08_Z_z_z_12197423",
     "08_P_z_pst":      "08_P_z_pst_12197423",
     "09_Z_z_z":        "09_Z_z_z_11196019",
     "10_Z_z_z":        "10_Z_z_z_11196413",
     "12_Z_z_z":        "12_Z_z_z_11196414",
     "norm7_norm7":     "norm7_norm7_12197008",
     "norm8_norm8":     "norm8_norm8_11198007",
}

id_to_b_values_dict = {
     #mc spec : mean start, b window

     "Z_m_p_0":            [5280, 50],
     "Z_m_p_1":            [5130, 60],
     "Z_m_p_2":            [4980, 80],

     "Z_z_z_0":            [5280, 50],
     "Z_z_z_1":            [5130, 70],
     "Z_z_z_2":            [4980, 90],

     "P_z_p_0":            [5280, 50],
     "P_z_p_1":            [5130, 70],
     "P_z_p_2":            [4980, 90],

     "M_m_z_1":            [5130, 70],
     "M_m_z_2":            [4980, 90],

     "P_z_pst_0":            [5280, 50],
     "P_z_pst_1":            [5130, 70],

     "norm7_0":         [5280, 50],
     "norm8_0":         [5280, 50],

     "01_Z_m_p":            [5280, 25, 0],
     "02_Z_m_p":            [5130, 40, 1],
     "02_P_z_p":            [5130, 50, 1],
     "04_Z_m_p":            [4985, 60, 2],
     "04_P_z_p":            [4980, 65, 2],
     "04_Z_z_z":            [4980, 60, 2],
     "04_P_z_pst":          [5130, 50, 1],
     "05_P_z_p":            [5280, 30, 0],
     "06_P_z_p":            [5130, 80, 1],
     "07_P_z_p":            [5130, 60, 1],
     "07_Z_z_z":            [5130, 60, 1],
     "07_P_z_pst":          [5280, 30, 0],
     "08_P_z_p":            [4985, 80, 2],
     "08_Z_z_z":            [4985, 80, 2],
     "08_P_z_pst":          [5130, 60, 1],
     "09_Z_z_z":            [5280, 30, 0],
     "10_Z_z_z":            [5130, 75, 1],
     "12_Z_z_z":            [4980, 90, 2],
     "norm7_norm7":            [5280, 25, 0],
     "norm8_norm8":            [5280, 25, 0],
}

ids_to_bestfit_dict = {
     "Z_m_p_0":             "DG",
     "Z_m_p_1":             "BGEP",
     "Z_m_p_2":             "G",

     "Z_z_z_0":             "G",
     "Z_z_z_1":             "BG",
     "Z_z_z_2":             "BG",

     "P_z_p_0":             "G",
     "P_z_p_1":             "BGEP",
     "P_z_p_2":             "BG",

     "M_m_z_1":             "BG",
     "M_m_z_2":             "BG",

     "P_z_pst_0":            "BGEP",
     "P_z_pst_1":            "BGEP",

     "norm7_0":             "G",
     "norm8_0":             "G",

     "norm7_norm7":         "DG",
     "norm8_norm8":         "DG",
     "01_Z_m_p":            "DG",
     "02_Z_m_p":            "BGEP",
     "02_P_z_p":            "BGEP",
     "04_Z_m_p":            "BGEP",
     "04_P_z_p":            "BGEP",
     "04_Z_z_z":            "BGEP",
     "04_P_z_pst":          "BGEP",
     "05_P_z_p":            "BGEP",
     "06_P_z_p":            "GEPAddBGEP_fr",
     "07_P_z_p":            "GAddBGEP_fr",
     "07_Z_z_z":            "BGEP",
     "07_P_z_pst":          "BGEP",
     "08_P_z_p":            "GEPAddBGEP_fr",
     "08_Z_z_z":            "GEPAddBGEP_fr",
     "08_P_z_pst":          "BGEP",
     "09_Z_z_z":            "BGEP",
     "10_Z_z_z":            "GEPAddBGEP_fr",
     "12_Z_z_z":            "GAddBG_fr",
}

id_to_genstring_dict = {
     "norm7_norm7" : [
                    "D0b_0_D0_0_M2",
                    "D0b_0_Kp_0_M2",
                    "D0_0_Kp_0_M2",
                    "Km_0_pip_0_M2",
                    "Km_0_pip_1_M2",
                    "Km_0_pim_1_M2",
                    "pip_0_pip_1_M2",
                    "pip_0_pim_1_M2",
                    "pip_1_pim_1_M2",
                    "Km_0_pip_0_pip_1_M2",
                    "Km_0_pip_0_pim_1_M2",
                    "Km_0_pip_1_pim_1_M2",
                    "pip_0_pip_1_pim_1_M2"],

     "norm8_norm8" : ["Dm_0_D0_0_M2","Dm_0_Kp_0_M2","D0_0_Kp_0_M2","Km_0_pip_0_M2","Km_0_pip_1_M2","Km_0_pim_2_M2","pip_0_pip_1_M2","pip_0_pim_2_M2","pip_1_pim_2_M2","Km_0_pip_0_pip_1_M2","Km_0_pip_0_pim_2_M2","Km_0_pip_1_pim_2_M2","pip_0_pip_1_pim_2_M2"],
     "01_Z_m_p":  ["Dm_0_Dp_0_M2","Dm_0_Kst0_0_M2","Dp_0_Kst0_0_M2"],
     "02_Z_m_p":  ["Dstm_0_Dp_0_M2","Dstm_0_Kst0_0_M2","Dp_0_Kst0_0_M2"],
     "04_Z_m_p":  ["Dstm_0_Dstp_0_M2","Dstm_0_Kst0_0_M2","Dstp_0_Kst0_0_M2"],
}

id_to_color_values_dict = {
     "norm7_pdf" :  "#d43d51",
     "norm7_signal" : "#00876c",
     "norm7_signal_a" : "#f4e07f",
     "norm7_signal_b" : "#ef9250",
     "norm7_bkg" : "#84b76e",

     "norm8_pdf" :  "#d43d51",
     "norm8_signal" : "#00876c",
     "norm8_signal_a" : "#f4e07f",
     "norm8_signal_b" : "#ef9250",
     "norm8_bkg" : "#84b76e",

     "Z_m_p_pdf":       "#d43d51",
     "Z_m_p_bkg":       "#84b76e",
     "01_Z_m_p":        "#00876c",
     "02_Z_m_p":        "#f4e07f",
     "04_Z_m_p":        "#ef9250",

     "Z_z_z_pdf":       "#d43d51",
     "Z_z_z_bkg":       "#84b76e",
     "01_Z_z_z":        "#00876c",
     "02_Z_z_z":        "#f4e07f",
     "04_Z_z_z":        "#ef9250",

     "P_z_p_pdf":       "#d43d51",
     "P_z_p_bkg":       "#84b76e",
     "01_P_z_p":        "#00876c",
     "02_P_z_p":        "#f4e07f",
     "04_P_z_p":        "#ef9250",

     "M_m_z_pdf":       "#d43d51",
     "M_m_z_bkg":       "#84b76e",
     "01_M_m_z":        "#00876c",
     "02_M_m_z":        "#f4e07f",
     "04_P_z_p":        "#ef9250",

     "P_z_pst_pdf":       "#d43d51",
     "P_z_pst_bkg":       "#84b76e",
     "01_P_z_p":        "#00876c",
     "02_P_z_p":        "#f4e07f",
     "04_P_z_p":        "#ef9250",

}

id_to_meson_string_dict = {
     #spec : , [B name, D1 name, D2 name, K*0 or K+]
     "norm7":           ["B^{+}", "#bar{D}^{0}", "D^{0} #rightarrow K+#pi-#pi-#pi+", "K^{+}"],
     "norm8":           ["B^{0}", "D^{-}", "D^{0} #rightarrow K+#pi-#pi-#pi+", "K^{+}"],

     "norm7_norm7":           ["B^{+}", "#bar{D}^{0}", "D^{0} #rightarrow K+#pi-#pi-#pi+", "K^{+}"],
     "norm8_norm8":           ["B^{0}", "D^{-}", "D^{0} #rightarrow K+#pi-#pi-#pi+", "K^{+}"],

     "Z_m_p":           ["B^{0}", "D^{-}", "D^{+}", "K^{*0}"],
     "Z_z_z":           ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
     "P_z_p":           ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "M_m_z":           ["B^{-}", "D^{-}", "D^{0}", "K^{*0}"],
     "P_z_pst":         ["B^{+}", "#bar{D}^{0}",  "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],

     "01_Z_m_p":        ["B^{0}", "D^{-}", "D^{+}", "K^{*0}"],
     "02_Z_m_p":        ["B^{0}", "D^{-}", "D^{+}", "K^{*0}"],
     "02_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "04_Z_m_p":        ["B^{0}", "D^{-}", "D^{+}", "K^{*0}"],
     "04_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "04_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
     "04_P_z_pst":      ["B^{+}", "#bar{D}^{0}", "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],
     "05_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "06_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "07_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "07_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
     "07_P_z_pst":      ["B^{+}", "#bar{D}^{0}", "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],
     "08_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "08_P_z_pst":      ["B^{+}", "#bar{D}^{0}", "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],
     "08_Z_z_z":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "09_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
     "10_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
     "12_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
}

id_to_pid_dict = {

     "01_Z_m_p":        [B0_ID, Dp_ID, Dp_ID, kst0_ID],
     "02_Z_m_p":        [B0_ID, Dp_ID, Dp_ID, kst0_ID],
     "02_P_z_p":        [B0_ID, Dp_ID, Dp_ID, kst0_ID],
     # "04_Z_m_p":        [B0_ID, Dp_ID, Dp_ID, kst0_ID],
     # "04_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     # "04_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
     # "04_P_z_pst":      ["B^{+}", "#bar{D}^{0}", "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],
     "05_P_z_p":        [Bp_ID, D0_ID, Dp_ID, kst0_ID],
     "06_P_z_p":        [Bp_ID, D0_ID, Dp_ID, kst0_ID],
     "07_P_z_p":        [Bp_ID, D0_ID, Dp_ID, kst0_ID],
     "07_Z_z_z":        [Bp_ID, D0_ID, D0_ID, kst0_ID],
     # "07_P_z_pst":      ["B^{+}", "#bar{D}^{0}", "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],
     # "08_P_z_p":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     # "08_P_z_pst":      ["B^{+}", "#bar{D}^{0}", "(D^{*+} #rightarrow D^{0}#pi^{+})", "K^{*0}"],
     # "08_Z_z_z":        ["B^{+}", "#bar{D}^{0}", "D^{+}", "K^{*0}"],
     "09_Z_z_z":        [B0_ID, D0_ID, D0_ID, kst0_ID],
     "10_Z_z_z":        [B0_ID, D0_ID, D0_ID, kst0_ID],

     "norm7_norm7": [Bp_ID, D0_ID, D0_ID, k_ID],
     "norm8_norm8": [B0_ID, Dp_ID, D0_ID, k_ID],

     # "12_Z_z_z":        ["B^{0}", "#bar{D}^{0}", "D^{0}", "K^{*0}"],
}

id_to_scheme_dict = {
    "Z_m_p":           "Z",
    "Z_z_z":           "ZZ",
    "P_z_p":           "P",
    "M_m_z":           "M",
    "P_z_pst":         "ST",
    "norm7":           "N7",
    "norm8":           "N8",
    "01_Z_m_p_11198006" : "1",
    "02_Z_m_p_11198400" : "2a,3a",
    "02_P_z_p_11198005" : "2b,3b",
    "04_Z_m_p_11198401" : "4a",
    "04_P_z_p_11198410" : "4b",
    "04_Z_z_z_11198023" : "4c",
    "04_P_z_pst_11198023" : "4d",
    "05_P_z_p_12197023" : "5",
    "06_P_z_p_12197410" : "6",
    "07_P_z_p_12197400" : "7a",
    "07_Z_z_z_12197045": "7b",
    "07_P_z_pst_12197045": "7c",
    "08_P_z_p_12197401" : "8a",
    "08_Z_z_z_12197423" : "8b",
    "08_P_z_pst_12197423" : "8c",
    "09_Z_z_z_11196019" : "9",
    "10_Z_z_z_11196413" : "10",
    "12_Z_z_z_11196414" : "12",
    "13_Zs_sm_p_13198040" : "13",
    "14_Zs_sm_p_13198200" : "14",
    "15_Zs_sm_p_13198400" : "15",
    "16_Zs_sm_p_13198600" : "16",
    "norm7_norm7_12197008": "norm7",
    "norm8_norm8_11198007" : "norm8"
}

id_to_legend_loc = {
    "Z_m_p":   [0.20, 0.575, 0.40, 0.905],
    "Z_z_z":   [0.20, 0.575, 0.40, 0.905],
    "P_z_p":   [0.20, 0.575, 0.40, 0.905],
    "M_m_z":   [0.20, 0.575, 0.40, 0.905],
    "P_z_pst":   [0.20, 0.575, 0.40, 0.905],
    "norm7":   [0.175, 0.575, 0.375, 0.905],
    "norm8":   [0.175, 0.575, 0.375, 0.905],
    "rc":  [0.60, 0.575, 0.80, 0.905],
    "rc2":  [0.70, 0.575, 0.90, 0.905],
    "lc":  [0.15, 0.575, 0.35, 0.955],
    #         legend = ROOT.TLegend(0.60, 0.55, 0.85, 0.85)
    #         name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield","Missing 2 Particle Yield"]
    #
    #     if spec == "P_z_p":
    #         legend = ROOT.TLegend(0.65, 0.625, 0.90, 0.925)
    #         name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield","Missing 2 Particle Yield"]
    #
    #     if spec == "M_m_z":
    #         legend = ROOT.TLegend(0.60, 0.55, 0.85, 0.85)
    #         name_list =  ["Missing 1 Particle Yield","Missing 2 Particle Yield"]
    #
    #     if spec == "P_z_pst":
    #         legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.875)
    #         name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield"]
    #
    #     if spec == "Zs_sm_p":
    #         legend = ROOT.TLegend(0.20, 0.575, 0.40, 0.875)
    #         name_list =  ["Fully Reconstruced Yield", "Missing 1 Particle Yield", "Missing 1 Particle Yield", "Missing 2 Particle Yield"]
    #

}

data_to_mc_dict = {
    "Z_m_p" : ["01_Z_m_p", "02_Z_m_p","04_Z_m_p"]
}
def uncor_eff(pre, post):
    fail = pre - post
    next = post / (post + fail)
    return (next)

class rdf_numbers():
    def __init__(self, rdf_in, name):
        self.rdf = rdf_in
        self.name = name

        base_count = self.rdf.Count().GetValue()
        base_count_err = math.sqrt(base_count)
        self.base_count_ufloat = ufloat(base_count, base_count_err)

        npy_temp = self.rdf.AsNumpy(columns=["eventNumber","runNumber"])
        df_temp = pd.DataFrame(npy_temp)
        unique_events = df_temp[~df_temp.duplicated()].value_counts()
        n_unique_events = unique_events.size
        n_unique_events_err = math.sqrt(n_unique_events)
        self.n_unique_events_ufloat = ufloat(n_unique_events, n_unique_events_err)

        nodups = df_temp[~df_temp.duplicated(keep=False)].value_counts()
        dups = df_temp[df_temp.duplicated(keep=False)].value_counts()

        # both = nodups.append(dups)
        both = pd.concat([nodups, dups], axis = 0, join = 'outer')

        can_per_event = both.mean()
        can_per_event_err = both.std()/math.sqrt(both.size)
        self.can_per_event_ufloat = ufloat(can_per_event, can_per_event_err)

    def calc_can_eff(self):
        self.fail_count_ufloat = self.old_stuff.base_count_ufloat - self.base_count_ufloat
        next_eff = (self.base_count_ufloat/(self.base_count_ufloat + self.fail_count_ufloat))
        return next_eff

    def calc_event_eff(self):
        base_e_count_ufloat = self.old_stuff.n_unique_events_ufloat
        pass_e_count_ufloat = self.n_unique_events_ufloat
        fail_e_count_ufloat = base_e_count_ufloat - pass_e_count_ufloat
        next_eff = (pass_e_count_ufloat/(pass_e_count_ufloat + fail_e_count_ufloat))
        return next_eff

    def apply_filter(self, filter, filter_name):
        next_rdf = self.rdf.Filter(filter, filter_name)
        next_rdf_stuff = rdf_numbers(next_rdf, filter_name)
        next_rdf_stuff.old_stuff = self
        return next_rdf_stuff

def get_decay_string(m_list, type):

    b = m_list[0]
    d1 = m_list[1]
    d2 = m_list[2]
    k = m_list[3]

    if type == "Full":
        decay_string = f"{b} #rightarrow {d1} {d2} {k}"
    if type == "Rec":
        decay_string = f"{d1} {d2} {k}"
    return decay_string

def grab_files_and_chain(file_path, DecayTree_List):
    file_list = glob.glob(file_path)
    tchain = ROOT.TChain()
    print(file_list)
    for file_name in file_list:
        for DecayTree in DecayTree_List:
            tchain.Add(f"{file_name}/{DecayTree}")
    return tchain

def get_roofit_pdf(ws, spec, shape_flag, varoi):

    if "BKG:" not in shape_flag:
        mean_start =  id_to_b_values_dict[spec][0]
        mean_low = mean_start - 20
        mean_high = mean_start + 20

    if shape_flag == "G":
        ws.factory(f"Gaussian::{spec}_fit_{shape_flag}({varoi}, mean_{spec}_{shape_flag}[{mean_start},{mean_low}, {mean_high}], width_{spec}_{shape_flag}[5.0,0.1,30.0])")
    if shape_flag == "DG":
        ws.factory(
            f"Gaussian::{spec}_fit_{shape_flag}_a({varoi}, mean_{spec}_{shape_flag}[{mean_start},{mean_low}, {mean_high}], width_a_{spec}_{shape_flag}[10, 0.1, 30.0])"
        )
        ws.factory(
            f"Gaussian::{spec}_fit_{shape_flag}_b({varoi}, mean_{spec}_{shape_flag}, width_b_{spec}_{shape_flag}[15.0, 0.1, 20.0])"
        )
        ws.factory(f"SUM::{spec}_fit_{shape_flag}(a_frac_{spec}_{shape_flag}[0.5, 0.1, 0.9]*{spec}_fit_{shape_flag}_a, {spec}_fit_{shape_flag}_b)")

    if shape_flag == "GEP":
        ws.factory(
            f"RooGaussExp::{spec}_fit_{shape_flag}({varoi},mean_{spec}_{shape_flag}[{mean_start},{mean_low}, {mean_high}],width_{spec}_{shape_flag}[10,4.0,30.0],alpha_{spec}_{shape_flag}[3,0.05,7.0])"
        )

    if shape_flag == "BGEP":
        ws.factory(
            f"RooBifurGaussExp::{spec}_fit_{shape_flag}({varoi},mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}],width_L_{spec}_{shape_flag}[10,1.0,30.0],width_R_{spec}_{shape_flag}[5,1.0,30.0],alpha_1_{spec}_{shape_flag}[2,0.01,3.0],alpha_2_{spec}_{shape_flag}[2,0.01,3.0])"
        )

    if shape_flag == "BG":
        ws.factory(
            f"BifurGauss::{spec}_fit_{shape_flag}({varoi},mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}],width_1_{spec}_{shape_flag}[20,0.01,50.0],width_2_{spec}_{shape_flag}[5,0.01,50.0])"
        )
    if shape_flag == "cb1R":
        ws.factory(
            f"CBShape::{spec}_fit_{shape_flag}(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}],width_{spec}_{shape_flag}[10,0.01,20],alpha_{spec}_{shape_flag}[-3,-5,-0.00001], n_{spec}_{shape_flag}[5,0,50])"
        )
    if shape_flag == "cb1L":
        ws.factory(
            f"CBShape::{spec}_fit_{shape_flag}(B_DTF_M,mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}],width_{spec}_{shape_flag}[10,0.01,20], alpha_{spec}_{shape_flag}[2,0.01,4.0], n_{spec}_{shape_flag}[5,0,50])"
        )
    #Gaussian + BifurGauss Shared Mean Floating Ratio
    if shape_flag == "GAddBG_fr":
        ws.factory(
            f"Gaussian:{spec}_fit_{shape_flag}_a(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}], width_{spec}_{shape_flag}_a[10, 1.0, 50.0])"
        )
        ws.factory(
            f"BifurGauss::{spec}_fit_{shape_flag}_b(B_DTF_M, mean_{spec}_{shape_flag}, width_1_{spec}_{shape_flag}_b[20,0.01,20.0],width_2_{spec}_{shape_flag}_b[5,0.01,25.0])"
        )
        ws.factory(f"SUM:{spec}_fit_{shape_flag}(a_frac_{spec}_{shape_flag}[0.647, 0.1, 0.95]*{spec}_fit_{shape_flag}_a, {spec}_fit_{shape_flag}_b)")
    # Gaussian w exp tail + BifurGauss Shared Mean Floating Ratio
    if shape_flag == "GEPAddBG_fr":
        ws.factory(
            f"RooGaussExp:{spec}_fit_{shape_flag}_a(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}], width_{spec}_{shape_flag}_a[10,1.0,50.0],alpha_{spec}_{shape_flag}_a[3,0.01,5.0])"
        )
        ws.factory(
            f"BifurGauss::{spec}_fit_{shape_flag}_b(B_DTF_M, mean_{spec}_{shape_flag},width_1_{spec}_{shape_flag}_b[20,0.01,40.0],width_2_{spec}_{shape_flag}_b[5,0.01,40.0])"
        )
        ws.factory(f"SUM:{spec}_fit_{shape_flag}(a_frac_{spec}_{shape_flag}[0.647, 0.1, 0.95]*{spec}_fit_{shape_flag}_a, {spec}_fit_{shape_flag}_b)")
    # Gaussian + BifurGauss w exp tail Shared Mean Floating Ratio

    if shape_flag == "GAddBGEP_fr":
        ws.factory(
            f"Gaussian::{spec}_fit_{shape_flag}_a(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}], width_{spec}_{shape_flag}_a[25.0,10.0,50.0])"
        )
        ws.factory(
            f"RooBifurGaussExp:{spec}_fit_{shape_flag}_b(B_DTF_M, mean_{spec}_{shape_flag},width_L_{spec}_{shape_flag}_b[10,0.01,35.0],width_R_{spec}_{shape_flag}_b[15,0.01,35.0],alpha_1_{spec}_{shape_flag}_b[2,0.001,5],alpha_2_{spec}_{shape_flag}_b[3,0.001,4])"
        )
        ws.factory(f"SUM:{spec}_fit_{shape_flag}(a_frac_{spec}_{shape_flag}[0.647, 0.1, 0.95]*{spec}_fit_{shape_flag}_a, {spec}_fit_{shape_flag}_b)")
    # #Gaussian w exp tail+ BifurGauss w exp tail Shared Mean Floating Ratio

    if shape_flag == "GEPAddBGEP_fr":
        ws.factory(
            f"RooGaussExp:{spec}_fit_{shape_flag}_a(B_DTF_M, mean_{spec}_{shape_flag}[{mean_start},{mean_low},{mean_high}], width_{spec}_{shape_flag}_a[25,3.0,50.0],alpha_{spec}_{shape_flag}_a[3,0.01,5])"
        )
        ws.factory(
            f"RooBifurGaussExp:{spec}_fit_{shape_flag}_b(B_DTF_M, mean_{spec}_{shape_flag},width_L_{spec}_{shape_flag}_b[10,1.0,40.0],width_R_{spec}_{shape_flag}_b[15,1.0,40.0],alpha_1_{spec}_{shape_flag}_b[2,0.05,5],alpha_2_{spec}_{shape_flag}_b[3,0.05,5])"
        )
        ws.factory(f"SUM:{spec}_fit_{shape_flag}(a_frac_{spec}_{shape_flag}[0.647, 0.1, 0.95]*{spec}_fit_{shape_flag}_a, {spec}_fit_{shape_flag}_b)")


    if shape_flag == "BKG: Exponential":
        ws.factory(f"Exponential:{spec}_fit_bkg({varoi}, c0_{spec}[0, -0.01, 0.01])")

    if shape_flag == "BKG: Chebychev":
        ws.factory(
            f"Chebychev:{spec}_fit_bkg({varoi},{{c0_{spec}[0.,-3,3],c1_{spec}[0.,-3,3], c2_{spec}[0.,-3,3]}})"
        )

    if shape_flag == "BKG: Bernstein":
        ws.factory(
            f"Bernstein:{spec}_fit_bkg({varoi},{{c0_{spec}[1,0,10], c1_{spec}[1,0,10], c2_{spec}[1,0,10], c3_{spec}[1,0,10]}})"
        )

def save_pdf(thing, folder, name, rpflag = 0):
    import datetime
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}/{folder}/'):
        os.makedirs(f'plots/{now.month}_{now.day}/{folder}/')
    if rpflag == 1:
        thing.save(f"plots/{now.month}_{now.day}/{folder}/{name}.pdf")
    else:
        thing.SaveAs(f"plots/{now.month}_{now.day}/{folder}/{name}.pdf")
    print(f"saved: {name}.pdf in plots/{now.month}_{now.day}/{folder}/")

def save_png(thing, folder, name, rpflag = 0):
    import datetime
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}/{folder}/'):
        os.makedirs(f'plots/{now.month}_{now.day}/{folder}/')
    if rpflag == 1:
        thing.save(f"plots/{now.month}_{now.day}/{folder}/{name}.png")
    else:
        thing.SaveAs(f"plots/{now.month}_{now.day}/{folder}/{name}.png")
    print(f"saved: {name}.png in plots/{now.month}_{now.day}/{folder}/")

def get_dwindow_values(spec, rflag):

    dst_flag = False
    if "Z_m_p" in spec:
        d1_flag = "mp"
        d2_flag = "mp"
    if "Z_z_z" in spec:
        d1_flag = "z"
        d2_flag = "z"
    if "P_z_p" in spec and "P_z_pst" not in spec:
        d1_flag = "z"
        d2_flag = "mp"
    if "M_m_z" in spec:
        d1_flag = "mp"
        d2_flag = "z"
    if "P_z_pst" in spec:
        d1_flag = "z"
        d2_flag = "z"
        dst_flag = True
    if "norm7" in spec:
        d1_flag = "z"
        d2_flag = "d0k3pi"
    if "norm8" in spec:
        d1_flag = "mp"
        d2_flag = "d0k3pi"
    if spec in ["z","mp","d0k3pi","dst"]:
        d1_flag = spec
        d2_flag = spec
        if spec == "dst":
            d1_flag = "z"
            d2_flag = "z"
            dst_flag = True

    d1window_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/d_window_files/d_{d1_flag}_mass_fits.root", "READ")
    d1window_ws = d1window_file.Get(f"d_{d1_flag}_mass_fits")
    d2window_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/d_window_files/d_{d2_flag}_mass_fits.root", "READ")
    d2window_ws = d2window_file.Get(f"d_{d2_flag}_mass_fits")

    d1_mstart = d1window_ws.var(f"mean").getValV()
    d2_mstart = d2window_ws.var(f"mean").getValV()

    if d1_flag == "d0k3pi":
        d1_window = 17.5
    if d2_flag == "d0k3pi":
        d2_window = 17.5
    if d1_flag == "mp":
        d1_window = 17
    if d2_flag == "mp":
        d2_window = 17
    if d1_flag == "z":
        d1_window = 17
    if d2_flag == "z":
        d2_window = 17

    d1_std = d1_window/2
    d2_std = d2_window/2

    d1_sig_max = 2*d1_std
    d2_sig_max = 2*d2_std

    d1_sb_min = 3*d1_std
    d2_sb_min = 3*d2_std

    d1_sb_max = 5*d1_std
    d2_sb_max = 5*d2_std

    d1_sb_min_line = f"(abs({d1_mstart} - D1_M) > {d1_sb_min})"
    d2_sb_min_line = f"(abs({d2_mstart} - D2_M) > {d2_sb_min})"

    d1_sb_max_line = f"(abs({d1_mstart} - D1_M) < {d1_sb_max})"
    d2_sb_max_line = f"(abs({d2_mstart} - D2_M) < {d2_sb_max})"

    d1_sb_line = f"({d1_sb_min_line} && {d1_sb_max_line})"
    d2_sb_line = f"({d2_sb_min_line} && {d2_sb_max_line})"

    d1_sig_line = f"(abs({d1_mstart} - D1_M) < {d1_sig_max})"
    d2_sig_line = f"(abs({d2_mstart} - D2_M) < {d2_sig_max})"

    sig_sig_line = f"({d1_sig_line} && {d2_sig_line})"
    sig_sb_line = f"({d1_sig_line} && {d2_sb_line})"
    sb_sig_line = f"({d1_sb_line} && {d2_sig_line})"
    sb2_line = f"({sig_sb_line} || {sb_sig_line})"

    if dst_flag:
        dstwindow_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/build_rootfiles/d_window_root_files/d_dst_mass_fits.root", "READ")
        dstwindow_ws = dstwindow_file.Get(f"d_dst_mass_fits")
        dst_mstart = dstwindow_ws.var(f"mean").getValV()
        dstwindow = 2.0
        sig_sig_line = f"({sig_sig_line}) && (abs(D2stmD_M - {dst_mstart}) < {dstwindow})"

    if rflag == "sig":
        return sig_sig_line
    if rflag == "sb":
        return sb2_line
    if rflag == "table":
        if spec != "dst":
            return d1_mstart, d1_std
        if spec == "dst":
            return dst_mstart, dstwindow/2
    if rflag == "charmless_start":
        return d1_mstart, d2_mstart
    if rflag == "charmless_2d":
        return d1_sb_min, d2_sb_min, d1_sb_max, d2_sb_max

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

    # ROOT.RooMsgService.instance().Print()
    # rme = ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.INFO)
    # ROOT.RooMsgService.instance().addStream(ROOT.RooFit.DEBUG, ROOT.RooFit.Topic(ROOT.RooFit.Contents), ROOT.RooFit.OutputFile("rfsmear_debug.log"))
    #projDeps = ROOT.RooArgSet(), useWeights = False, copyDataSet = True, newName = "test", fitToarg5 =  ROOT.RooFit.PrintLevel(3)

    sData = ROOT.RooStats.SPlot("sData","An SPlot", data, model, yields)
    print("npds")

    # print("Check SWeights:")
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

def build_truth_base(spec):

    b_spec = id_to_pid_dict[spec][0]
    c1_spec = id_to_pid_dict[spec][1]
    c2_spec = id_to_pid_dict[spec][2]

    if id_to_pid_dict[spec][3] == k_ID:
        norm_flag = 1
    else:
        norm_flag = 0

    truth_MC = f"abs(B_TRUEID) == {b_spec}"
    truth_MC = f"{truth_MC} && abs(D1H1_TRUEID) == {k_ID}"
    truth_MC = f"{truth_MC} && abs(D1H2_TRUEID) == {pi_ID}"
    if c1_spec == Dp_ID and b_spec != B0_ID:
        truth_MC = f"{truth_MC} && abs(D1H3_TRUEID) == {pi_ID}"
    truth_MC = f"{truth_MC} && abs(D2H1_TRUEID) == {k_ID}"
    truth_MC = f"{truth_MC} && abs(D2H2_TRUEID) == {pi_ID}"
    if c2_spec == Dp_ID:
        truth_MC = f"{truth_MC} && abs(D2H3_TRUEID) == {pi_ID}"
    if norm_flag == 1:
        truth_MC = f"{truth_MC} && abs(D2H4_TRUEID) == {pi_ID}"
        truth_MC = f"{truth_MC} && abs(K_TRUEID) == {k_ID}"
    if norm_flag == 0:
        truth_MC = f"{truth_MC} && abs(KSTH1_TRUEID) == {k_ID}"
        truth_MC = f"{truth_MC} && abs(KSTH2_TRUEID) == {pi_ID}"
    print(truth_MC)
    return truth_MC
