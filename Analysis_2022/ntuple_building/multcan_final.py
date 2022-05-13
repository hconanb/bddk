import sys
import os

basedir = os.getcwd().split('ntuple_building')[0]
sys.path.append(basedir)

from essential_functions import *

def getBestCand_old(spec, type, year, snap_flag, txt_flag):

    print(f"Getting Best candidates for {spec} in {year} as {type}")


    file = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_veto/{spec}.root"


    with uproot.open(f"{file}:DecayTreeTuple_SIG") as sigin, uproot.open(f"{file}:DecayTreeTuple_SB") as sbin:

        ###Base Files
        d_sig_base = sigin.arrays(sigin.keys(), library="pd")[0]
        d_sb_base = sbin.arrays(sbin.keys(), library="pd")[0]
        ##Only read in branches for mult candidate selection
        d_sig_muc = sigin.arrays(["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2", "B_M","D1_M","D2_M"], library="pd")
        d_sb_muc = sbin.arrays(["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2", "B_M","D1_M","D2_M"], library="pd")
        cans_sb = len(d_sb_muc.index)

        dsig_bestCands = d_sig_muc.sort_values('B_dtf_chi2', ascending=False).drop_duplicates(subset=['runNumber' ,'eventNumber']).sort_index().reset_index(level=0, drop=True)
        dsb_bestCands = d_sb_muc.sample(frac=1).drop_duplicates(subset=['runNumber' ,'eventNumber']).sort_index().reset_index(level=0, drop=True)

        print(f"Candidates cut in signal region and sb region: {len(d_sig_base.index)}, {len(d_sb_base)}")

        print(f"Candidates cut with dtf selection in signal region: {len(d_sig_muc.index) - len(dsig_bestCands.index)}")
        print(f"Candidates cut with random selection in Sideband region {len(d_sb_muc.index) - len(dsb_bestCands.index)}")

        ###Fine the intersecting candidates between sigbc and sbbc
        df_merge = pd.merge(dsig_bestCands, dsb_bestCands, on=["runNumber", "eventNumber"], how='inner')
        print(f"Overlaping Events between Signal and Sideband region after selection {len(df_merge.index)}")
        # delete columns we dont need
        # del df_merge["nCandidate_x"]
        # del df_merge["B_dtf_chi2_x"]
        # del df_merge["B_M_x"]

        # rename columns for later merge
        eolap = len(df_merge.index)
        # df_merge = df_merge.rename(columns={"nCandidate_y": "nCandidate", "B_dtf_chi2_y": "B_dtf_chi2", "B_M_y" : "B_M"})
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

        df_sig_output = pd.merge(d_sig_base, dsig_bestCands, on=["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2", 'B_M', 'D1_M', 'D2_M'], how='inner')
        df_sb_output = pd.merge(d_sb_base, d_sb_muc_final, on=["runNumber", "eventNumber", "nCandidate", "B_dtf_chi2", 'B_M','D1_M', 'D2_M'], how='inner')

        if snap_flag:

            output_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/final_sample/{spec}.root"
            print(f"Recreating {output_path}")
            fileOut = uproot.recreate(output_path)

            fileOut["DecayTreeTuple_SIG"] = df_sig_output
            fileOut["DecayTreeTuple_SB"] = df_sb_output

            fileOut["DTT_m"] = df_merge

            df_T = df_sig_output[(df_sig_output['B_L0HadronDecision_TOS'] == 1) | (df_sig_output['B_L0MuonDecision_TOS'] == 1) | (df_sig_output['B_L0PhotonDecision_TOS'] == 1) | (df_sig_output['B_L0ElectronDecision_TOS'] == 1)]
            df_nTaT = df_sig_output[((df_sig_output['B_L0HadronDecision_TOS'] == 0) & (df_sig_output['B_L0MuonDecision_TOS'] == 0) & (df_sig_output['B_L0PhotonDecision_TOS'] == 0) & (df_sig_output['B_L0ElectronDecision_TOS'] == 0)) & ((df_sig_output['B_L0HadronDecision_TIS'] == 1) | (df_sig_output['B_L0MuonDecision_TIS'] == 1) | (df_sig_output['B_L0PhotonDecision_TIS'] == 1) | (df_sig_output['B_L0ElectronDecision_TIS'] == 1))]

            print(f"Unique Events in Signal {len(df_sig_output.index)} \n")
            print(f"Unique Events in Sb {len(df_sb_output.index)} \n")
            print(f"Unique TOS Events in Signal: {len(df_T.index)} \n")
            print(f"Unique TIS Events in Signal: {len(df_nTaT.index)}\n")

            fileOut["DecayTreeTuple_T"] = df_T
            fileOut["DecayTreeTuple_nTaT"] = df_nTaT

        if txt_flag:

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

            text_df = pd.DataFrame(dict_temp , index=[0])
            text_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/mult_can/{spec}_{year}.txt")

def getBestCand(spec, type, year, snap_flag, txt_flag):

    print(f"Getting Best candidates for {spec} in {year} as {type}")

    file = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/post_veto/{spec}.root"

    with uproot.open(f"{file}:DecayTreeTuple_SIG") as sigin, uproot.open(f"{file}:DecayTreeTuple_SB") as sbin:

        ###Base Files
        d_sig_base = sigin.arrays(sigin.keys(), library="pd")[0]
        d_sb_base = sbin.arrays(sbin.keys(), library="pd")[0]

        df_all = pd.concat([d_sig_base, d_sb_base], keys=["SIG", "SB"])

        df_no_mc = df_all.sample(frac=1).drop_duplicates(subset=['runNumber' ,'eventNumber']).sort_index()
        #.reset_index(level=0, drop=True)

        og_sig_can = len(d_sig_base.index)
        og_sb_can = len(d_sb_base.index)
        og_all_can = len(df_all.index)
        events_all = len(df_no_mc.index)

        print(f"Candidates in signal region and sb region: {og_sig_can}, {og_sb_can}")
        print(f"Candidates in all: {og_all_can}")
        print(f"Events after random selection in all: {events_all}")

        df_sig_final = df_no_mc.loc[["SIG"]].reset_index(drop=True)
        df_sb_final = df_no_mc.loc[["SB"]].reset_index(drop=True)
        #
        sig_all = len(df_sig_final.index)
        sb_all = len(df_sb_final.index)

        print(f"Events in Sig after Random Selection: {sig_all}")
        print(f"Events in SB after Random Selection: {sb_all}")

        print(df_sig_final)

        if snap_flag:

            output_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_{type}/{year}/final_sample/{spec}.root"

            print(f"Recreating {output_path}")

            fileOut = uproot.recreate(output_path)

            fileOut["DecayTreeTuple_SIG"] = df_sig_final
            fileOut["DecayTreeTuple_SB"] = df_sb_final

            trigger_TOS_ALL_L0 = "(B_L0HadronDecision_TOS == 1 | B_L0MuonDecision_TOS == 1 |  B_L0ElectronDecision_TOS == 1 | B_L0PhotonDecision_TOS == 1)"
            trigger_TIS_ALL_L0 = "(B_L0HadronDecision_TIS == 1 | B_L0MuonDecision_TIS == 1 |  B_L0ElectronDecision_TIS == 1 | B_L0PhotonDecision_TIS == 1)"

            df_T = df_sig_final.query(trigger_TOS_ALL_L0)
            df_nTaT = df_sig_final.query(f"~{trigger_TOS_ALL_L0} & {trigger_TIS_ALL_L0}")

            fileOut["DecayTreeTuple_T"] = df_T
            fileOut["DecayTreeTuple_nTaT"] = df_nTaT

        if txt_flag:

            dict_temp = {
                        'Spectrum' : spec,
                        'Year': year,
                        'Type': type,
                        'Canidates in Signal OG': og_sig_can,
                        'Canidates in Sideband OG': og_sb_can,
                        'Canidates Cut w random selection over All': og_all_can - events_all,
                        'Events in Sig after Random Selection': sig_all,
                        'Events in SB after Random Selection': sb_all
                        }

            text_df = pd.DataFrame(dict_temp , index=[0])
            text_df.to_csv(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/eff_txt_files/mult_can/{spec}_{year}.txt")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Apply Multiple Canidate Selection")
    parser.add_argument("--Spec_List", choices = id_to_spec_dict.keys(),  nargs="+", help = 'Spec')
    parser.add_argument('--Snap_Flag', action='store_true')
    parser.add_argument('--Txt_Flag', action='store_true')

    args = parser.parse_args()
    spec_list = args.Spec_List
    snap_flag = args.Snap_Flag
    txt_flag = args.Txt_Flag

    for spec in spec_list:
        if id_to_spec_dict[spec] == spec:
            type = "data"
        else:
            type = "mc"
        for year in ["2016","2017","2018"]:
            getBestCand(id_to_spec_dict[spec], type, year, snap_flag, txt_flag)
