import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
import os
import ROOT
import pandas as pd
import uproot
import glob
import datetime
import uproot

class spectrum_class:

  def __init__(self, name, d1_string, d2_string, d1_mass, d2_mass, b_string):
    self.spec = name

    self.d1_string = d1_string
    self.d2_string = d2_string

    self.d1_mass = d1_mass
    self.d2_mass = d2_mass

    if self.spec == "st":
        self.d3_string = "(D^{*+} #rightarrow D^{0} #pi+)"
        self.d3_mass = 150

    self.b_string = b_string

z_c = spectrum_class("Z_m_p", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
zz_c = spectrum_class("Z_z_z", "#barD^{0}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D^{0} K*0")
p_c = spectrum_class("P_z_p", "#bar{D^{0}}", "D^{+}", d0mass, dpmass, "B #rightarrow #bar{D^{0}} D+ K*0")
m_c = spectrum_class("M_m_z", "D^{-}", "D^{0}", dpmass, d0mass,  "B #rightarrow D- D^{0} K*0")
st_c = spectrum_class("P_z_pst","#bar{D^{0}}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D*+ K*0")
# # s_c = spectrum_class("s","D^{-}_{s}", "D^{+}" ,dsmass, dpmass, "e_g")
n7_c = spectrum_class("norm7","#bar{D^{0}}", "D^{0} #rightarrow k#pi#pi#pi", d0mass, d0mass,"B #rightarrow D-#bar{D^{0}}K+")
n8_c = spectrum_class("norm8","D^{-}", "D^{0} #rightarrow k#pi#pi#pi", dpmass, d0mass, "B #rightarrow D-#bar{D^{0}}K+")


# z_c_01 = spectrum_class("01_Z_m_p_11198006", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
# z_c_01 = spectrum_class("01_Z_m_p_11198006", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
# z_c_02 = spectrum_class("02_Z_m_p_11198400", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
# z_c_04 = spectrum_class("04_Z_m_p_11198401", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")

def getBestCand(type, year, spec):
    print(f"Getting Best candidates for {spec}, {year}")
    if type == "mc":
        dsig_file = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/sig_d/{spec}.root"
    if type == "data":
        dsig_file = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_d/{spec}.root"

    dsb_file = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/sb_d/{spec}.root"

    with uproot.open(f"{dsig_file}:DecayTreeTuple") as sigin, uproot.open(f"{dsb_file}:DecayTreeTuple") as sbin:

        ###Base Files
        d_sig_base = sigin.arrays(sigin.keys(), library="pd")[0]
        d_sb_base = sbin.arrays(sbin.keys(), library="pd")[0]
        ##Only read in branches for mult candidate selection
        d_sig_muc = sigin.arrays(["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2"], library="pd")
        d_sb_muc = sbin.arrays(["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2"], library="pd")
        cans_sb = len(d_sb_muc.index)

        dsig_bestCands = d_sig_muc.sort_values('B_dtf_chi2', ascending=False).drop_duplicates(subset=['runNumber' ,'eventNumber']).sort_index().reset_index(level=0, drop=True)
        dsb_bestCands = d_sb_muc.sample(frac=1).drop_duplicates(subset=['runNumber' ,'eventNumber']).sort_index().reset_index(level=0, drop=True)

        og_can_sig = len(d_sig_muc.index)
        og_can_sb = len(d_sb_muc.index)
        cc_dtf = len(d_sig_muc.index) - len(dsig_bestCands.index)
        cc_r = len(d_sb_muc.index) - len(dsb_bestCands.index)

        print(f"candidates cut with dtf selection {len(d_sig_muc.index) - len(dsig_bestCands.index)}")
        print(f"candidates cut with random selection {len(d_sb_muc.index) - len(dsb_bestCands.index)}")

        ###Fine the intersecting candidates between sigbc and sbbc
        df_merge = pd.merge(dsig_bestCands, dsb_bestCands, on=["runNumber", "eventNumber"], how='inner')
        print(f"Overlaping Events between sig and sb {len(df_merge.index)}")
        # delete columns we dont need
        del df_merge["nCandidate_x"]
        del df_merge["B_dtf_chi2_x"]
        # rename columns for later merge
        eolap = len(df_merge.index)

        df_merge = df_merge.rename(columns={"nCandidate_y": "nCandidate", "B_dtf_chi2_y": "B_dtf_chi2"})
        # print(df_merge)
        ###Append those to sb
        d_sb_olap = dsb_bestCands.append(df_merge).copy()
        # keep=False marks all duplicated row with a True
        d_sb_olap['Duplicated'] = d_sb_olap.duplicated(subset=['runNumber' ,'eventNumber'], keep=False)
        # print(d_sb_olap)
        # # # selects only rows which are not duplicated.
        d_sb_muc_final = d_sb_olap[~d_sb_olap['Duplicated']]
        # print(d_sb_muc_final)
        # delete the indicator column
        del d_sb_muc_final['Duplicated']

        # # print(f"Total candidates left in sb {len(d_sb_muc_final.index)}")
        # # print(f"Total cans eliminated check {cans_sb - len(d_sb_muc_final.index)}")

        df_sig_output = pd.merge(d_sig_base, dsig_bestCands, on=["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2"], how='inner')
        df_sb_output = pd.merge(d_sb_base, d_sb_muc_final, on=["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2"], how='inner')

        fileOut = uproot.recreate(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/final_sample/{spec}.root")
        fileOut["DecayTreeTuple_SIG"] = df_sig_output
        fileOut["DecayTreeTuple_SB"] = df_sb_output

        og_can_sig = len(d_sig_muc.index)
        og_can_sb = len(d_sb_muc.index)
        cc_dtf = len(d_sig_muc.index) - len(dsig_bestCands.index)
        cc_r = len(d_sb_muc.index) - len(dsb_bestCands.index)

        dict_temp = {
                    'Spectrum' : spec,
                    'Year': year,
                    'Type': type,
                    'Canidates/Unique Event In Signal': og_can_sig,
                    'Canidates/Unique Event In Sideband': og_can_sb,
                    'Canidates Cut w dtfx2 selection in Signal': cc_dtf,
                    'Canidates Cut w random selection in Sideband': cc_r,
                    'Event Overlap: Signal and Sideband': eolap,
                    'Final Events in Signal' : len(df_sig_output.index),
                    'Final Events in Sideband': len(df_sb_output.index),
                    }

        print(dict_temp)
        text_df = pandas.DataFrame(dict_temp , index=[0])
        text_df.to_csv(f"charmless_text/{spec}_{year}_{type}_cc.txt")

for year in [2016, 2017, 2018]:
    for cspec in [z_c, zz_c, p_c, m_c, st_c, n7_c, n8_c]
        getBestCand("data", year, cspec.spec)



    # getBestCand("mc", year, z_c_02.spec)
    # getBestCand("mc", year, z_c_04.spec)
