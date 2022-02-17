import pandas
import glob as glob
import ROOT
import numpy as np
from uncertainties import ufloat
from uncertainties import covariance_matrix, correlation_matrix
from uncertainties import unumpy
from uncertainties import ufloat_fromstr
import datetime
import os

RDF = ROOT.ROOT.RDataFrame

mc_spec_dict = {
    # "01_Z_m_p_11198006" : ("1", "B^{0} #rightarrow D^{-} D^{+} K^{*0}"),
    # "02_Z_m_p_11198400" : ("2a,3a", "B^{0} #rightarrow D^{*-} D^{+} K^{*0} / D^{-} D^{*+} K^{*0}"),
    # # "02_P_z_p_11198005" : ("2b,3b","B^{0} #rightarrow  (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) D^{+} K^{*0} "),
    # "04_Z_m_p_11198401" : ("4a","B^{0} #rightarrow D^{*-} D^{*+} K^{*0}"),
    # "04_P_z_p_11198410" : ("4b","B^{0} #rightarrow (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) (D^{*+} #rightarrow D^{+} #pi^{0})"),
    # "04_Z_z_z_11198023" : ("4c","B^{0} #rightarrow (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) (D^{*+} #rightarrow D^{0} #pi^{+})"),
    # "04_P_z_pst_11198023" : ("4d","B^{0} #rightarrow (D^{*-} #rightarrow #bar{D^{0}} #pi^{-}) (D^{*+} #rightarrow D^{0} #pi^{+})"),
    # # "05_P_z_p_12197023" : ("5","B^{+} #rightarrow #bar{D^{0}} D^{+} K^{*0}"),
    # # "06_P_z_p_12197410" : ("6","B^{+} #rightarrow #bar{D^{*0}} D^{+} K^{*0}"),
    # "07_P_z_p_12197400" : ("7a","B^{+} #rightarrow #bar{D^{0}} D^{*+} K^{*0}"),
    # "07_Z_z_z_12197045": ("7b","B^{+} #rightarrow #bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi^{+})"),
    # "07_P_z_pst_12197045": ("7c","B^{+} #rightarrow #bar{D^{*0}} D^{*+} K^{*0}"),
    # "08_P_z_p_12197401" : ("8a","B^{+} #rightarrow #bar{D^{*0}} D^{*+} K^{*0}"),
    # "08_Z_z_z_12197423" : ("8b","B^{+} #rightarrow (#bar{D^{0}} #rightarrow #bar{D^{0}} #bar{#pi^{0}}) (D^{*+} #rightarrow D^{0} #pi^{+})"),
    # "08_P_z_pst_12197423" : ("8c","B^{+} #rightarrow #bar{D^{0}} D^{*+} K^{*0}"),
    # "09_Z_z_z_11196019" : ("9", "B^{0} #rightarrow #bar{D^{0}} D^{0} K^{*0}"),
    # "10_Z_z_z_11196413" : ("10", "B^{0} #rightarrow #bar{D^{*0}} D^{0} K^{*0}"),
    # "12_Z_z_z_11196414" : ("12", "B^{0} #rightarrow #bar{D^{*0}} D^{*0} K^{*0}"),
    # "13_Zs_sm_p_13198040" : ("13", "B^{0}_{s} #rightarrow #bar{D^{-}_{s}} D^{+} K^{*0}"),
    # "14_Zs_sm_p_13198200" : ("14", "B^{0}_{s} #rightarrow #bar{D^{*-}_{s}} D^{+} K^{*0}"),
    # "15_Zs_sm_p_13198400" : ("15", "B^{0}_{s} #rightarrow #bar{D^{-}_{s}} D^{*+} K^{*0}"),
    # "16_Zs_sm_p_13198600" : ("16", "B^{0}_{s} #rightarrow #bar{D^{-}_{*s}} D^{*+} K^{*0}"),
    "norm7_norm7_12197008": ("norm7","B^{+} #rightarrow #bar{D^{0}} D^{*+} K^{*0}"),
    "norm8_norm8_11198007" : ("norm8", "B^{0} #rightarrow D^{-} (D^{0} -> #K-#pi+#pi+#pi-) K+")
}

for key in mc_spec_dict:
    for year in ["2016","2017","2018"]:
        file_id = key

        mc_file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/{year}/final_sample/{file_id}.root"
        mcid = key.split("_")[:-1]
        mcid = "_".join(mcid)

        r_file = ROOT.TFile(mc_file_path)
        t_tree = r_file.Get("DecayTreeTuple_T")
        ntat_tree = r_file.Get("DecayTreeTuple_nTaT")

        rdf_rec_T = RDF(t_tree)
        rdf_rec_nTaT = RDF(ntat_tree)

        t_count = rdf_rec_T.Count().GetValue()
        ntat_count = rdf_rec_nTaT.Count().GetValue()

        bootrec_t = 0
        bootrec_ntat = 0

        for boot_tree, trigger in zip((t_tree, ntat_tree), ("TOS","TIS")):

            bootrec_count = np.zeros(1000)
            start_flag = 0

            if trigger == "TOS":
                rec_count = t_count

            if trigger == "TIS":
                rec_count = ntat_count

            for ev in boot_tree:
                if start_flag == 0:
                    e_value = ev.RD_org_eventNumber
                    r_value = ev.RD_org_runNumber
                    count = 0
                    start_flag = 1
                if e_value != ev.RD_org_eventNumber or r_value != ev.RD_org_runNumber:
                    e_value = ev.RD_org_eventNumber
                    r_value = ev.RD_org_runNumber
                    bootrec_count += count*np.random.poisson(1,1000)
                    count = 0
                count += 1

            bootrec_count += count*np.random.poisson(1,1000)

            bootrec_all = ufloat(bootrec_count.mean(), bootrec_count.std())
            rec_all = ufloat(rec_count, np.sqrt(rec_count))

            if trigger == "TOS":
                bootrec_t = bootrec_all
            if trigger == "TIS":
                bootrec_ntat = bootrec_all

        dict_temp = {
                    'Scheme ID' : mc_spec_dict[file_id][0],
                    'Decay Description' :mc_spec_dict[file_id][1],
                    'Year': year,
                    f'Final Bootstraped TOS Events': f"{bootrec_t:.3f}",
                    f'Final Bootstraped TIS Events': f"{bootrec_ntat:.3f}",
                    }

        print(dict_temp)

        # now = datetime.datetime.now()
        # if not os.path.exists(f'boot_txt_files/{now.month}_{now.day}/'):
        #     os.makedirs(f'boot_txt_files/{now.month}_{now.day}/')

        df_dict = pandas.DataFrame(dict_temp, index=[0])
        outefftxt = f"boot_effs/{file_id}_{year}"
        df_dict.to_csv(outefftxt)
