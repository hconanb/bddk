import numpy
import root_numpy
import pandas
from hep_ml import reweight
from matplotlib import pyplot as plt
import glob as glob
import ROOT as ROOT
from Analysis_2021.essentials import save_png
import numpy
import os
RDF = ROOT.ROOT.RDataFrame

PE_d1d2 = "(pow(D1_PE + D2_PE,2))"
PX_d1d2 = "(pow(D1_PX + D2_PX,2))"
PY_d1d2 = "(pow(D1_PY + D2_PY,2))"
PZ_d1d2 = "(pow(D1_PZ + D2_PZ,2))"
M2_d1d2 = f"{PE_d1d2} - ({PX_d1d2} + {PY_d1d2} + {PZ_d1d2})"
M_d1d2 = f"sqrt({M2_d1d2})"

PE_d1kst = "(pow(D1_PE + KST_PE,2))"
PX_d1kst = "(pow(D1_PX + KST_PX,2))"
PY_d1kst = "(pow(D1_PY + KST_PY,2))"
PZ_d1kst = "(pow(D1_PZ + KST_PZ,2))"
M2_d1kst = f"{PE_d1kst} - ({PX_d1kst} + {PY_d1kst} + {PZ_d1kst})"
M_d1kst = f"sqrt({M2_d1kst})"

PE_d2kst = "(pow(D2_PE + KST_PE,2))"
PX_d2kst = "(pow(D2_PX + KST_PX,2))"
PY_d2kst = "(pow(D2_PY + KST_PY,2))"
PZ_d2kst = "(pow(D2_PZ + KST_PZ,2))"
M2_d2kst = f"{PE_d2kst} - ({PX_d2kst} + {PY_d2kst} + {PZ_d2kst})"
M_d2kst = f"sqrt({M2_d2kst})"

#########################################MCDecayTreeTuple
PE_d1d2_mcdtt = "(pow(C1_TRUEP_E + C2_TRUEP_E,2))"
PX_d1d2_mcdtt = "(pow(C1_TRUEP_X + C2_TRUEP_X,2))"
PY_d1d2_mcdtt = "(pow(C1_TRUEP_Y + C2_TRUEP_Y,2))"
PZ_d1d2_mcdtt = "(pow(C1_TRUEP_Z + C2_TRUEP_Z,2))"
M2_d1d2_mcdtt = f"{PE_d1d2_mcdtt} - ({PX_d1d2_mcdtt} + {PY_d1d2_mcdtt} + {PZ_d1d2_mcdtt})"
M_d1d2_mcdtt = f"sqrt({M2_d1d2_mcdtt})"

PE_d1kst_mcdtt = "(pow(C1_TRUEP_E + Kst_TRUEP_E,2))"
PX_d1kst_mcdtt = "(pow(C1_TRUEP_X + Kst_TRUEP_X,2))"
PY_d1kst_mcdtt = "(pow(C1_TRUEP_Y + Kst_TRUEP_Y,2))"
PZ_d1kst_mcdtt = "(pow(C1_TRUEP_Z + Kst_TRUEP_Z,2))"
M2_d1kst_mcdtt = f"{PE_d1kst_mcdtt} - ({PX_d1kst_mcdtt} + {PY_d1kst_mcdtt} + {PZ_d1kst_mcdtt})"
M_d1kst_mcdtt = f"sqrt({M2_d1kst_mcdtt})"

PE_d2kst_mcdtt = "(pow(C2_TRUEP_E + Kst_TRUEP_E,2))"
PX_d2kst_mcdtt = "(pow(C2_TRUEP_X + Kst_TRUEP_X,2))"
PY_d2kst_mcdtt = "(pow(C2_TRUEP_Y + Kst_TRUEP_Y,2))"
PZ_d2kst_mcdtt = "(pow(C2_TRUEP_Z + Kst_TRUEP_Z,2))"
M2_d2kst_mcdtt = f"{PE_d2kst_mcdtt} - ({PX_d2kst_mcdtt} + {PY_d2kst_mcdtt} + {PZ_d2kst_mcdtt})"
M_d2kst_mcdtt = f"sqrt({M2_d2kst_mcdtt})"

M2_d1d2_norm = M2_d1d2.replace("KST","K")
M2_d1kst_norm = M2_d1kst.replace("KST","K")
M2_d2kst_norm = M2_d2kst.replace("KST","K")

dalitz_columns = ["M2_d1d2", "M2_d2k", "M2_d1k"]
dalitz_columns_dict = {
"M2_d1d2" : "D^{-}D^{+}",
"M2_d2k" : "D^{+}K^{*0}",
"M2_d1k" : "D^{-}K^{*0}"
}
import ROOT

from hep_ml.metrics_utils import ks_2samp_weighted
hist_settings = {'bins': 20, 'density': True, 'alpha': 0.7}
# from Analysis_2021.essentials import save_png

ROOT.gROOT.ProcessLine(".L lhcbStyle.C")
ROOT.gSystem.Load("/home/hbernste/lhcb-analysis-master/rootclasses/lib/librootclasses.so")
ROOT.gStyle.SetPalette(ROOT.kBird)


def grab_files_and_chain(file_path, chain_name, type):
    file_list = glob.glob(file_path)
    print(file_list)
    if type != "MC":
        tchain = ROOT.TChain(chain_name)
        for file_name in file_list:
            tchain.Add(file_name)
        return tchain
    if type == "MC":
        tc_T = ROOT.TChain(f"DecayTreeTuple_T")
        tc_nTaT = ROOT.TChain(f"DecayTreeTuple_nTaT")
        for file_name in file_list:
            tc_T.Add(file_name)
            tc_nTaT.Add(file_name)
        return tc_T, tc_nTaT
def draw_sw_test(DATA, sw_np, MC, mc_weights, spec, mcid, year, savename):
    c1 = ROOT.TCanvas("c1","c1")
    for column in dalitz_columns:
        bins = 20

        data_np = DATA[column].to_numpy()
        mc_np = MC[column].to_numpy()

        data_df = ROOT.RDF.MakeNumpyDataFrame({column:data_np, "data_sweights": sw_np})
        mc_df = ROOT.RDF.MakeNumpyDataFrame({column:mc_np, "mc_weights":  mc_weights})

        xmin = data_df.Min(column).GetValue()
        xmax = data_df.Max(column).GetValue()

        print(xmin, xmax)

        hist_data_sw = data_df.Histo1D((f"data_{column}_sw", f"data_{column}_sw", bins, xmin, xmax), column, "data_sweights")
        # hist_data = data_df.Histo1D((f"data_{column}", f"data_{column}", bins, xmin, xmax), column)

        hist_mc = mc_df.Histo1D((f"mc_{column}", f"mc_{column}", bins, xmin, xmax), column, "mc_weights")

        # temp_data = hist_data.DrawCopy("E")
        temp_data_sw = hist_data_sw.DrawCopy("E")
        temp_mc = hist_mc.DrawCopy("E SAME")

        # temp_data.SetMarkerColor(ROOT.kBlue)
        # temp_data.Scale(1/temp_data.Integral())

        temp_data_sw.SetMarkerColor(ROOT.kCyan)
        temp_data_sw.Scale(1/temp_data_sw.Integral())

        temp_mc.SetMarkerColor(ROOT.kRed)
        temp_mc.Scale(1/temp_mc.Integral())

        plot_min = (temp_data_sw.GetXaxis().GetXmin()/1000000)
        plot_max = (temp_data_sw.GetXaxis().GetXmax()/1000000)

        temp_data_sw.GetXaxis().SetTitle(f"m^{{2}}_{{({dalitz_columns_dict[column]})}} [GeV^{{2}}]")
        # temp_data.GetXaxis().SetLimits(plot_min, plot_max)
        temp_data_sw.GetXaxis().SetLimits(plot_min, plot_max)
        temp_mc.GetXaxis().SetLimits(plot_min, plot_max)

        legend = ROOT.TLegend(0.70, 0.60, 0.95, 0.95)
        legend.AddEntry(temp_data_sw, f"Data_SWeighted","lp")
        # legend.AddEntry(temp_data, f"Data","lp")
        legend.AddEntry(temp_mc, f"MC","lp")
        legend.Draw()

        save_png(c1, f"{savename}", f"{spec}_{mcid}_{column}_{year}", rpflag = 0)



    # data_np = DATA["B_DTF_M"].to_numpy()
    # data_df = ROOT.RDF.MakeNumpyDataFrame({"B_DTF_M":data_np, "data_sweights_pds": sweights_np_pds, "data_sweights_npds": sweights_np_npds})
    # c2 = ROOT.TCanvas("c2","c2")
    # hist_data_pds_b = data_df.Histo1D((f"B_DTF_M_pds", f"B_DTF_M_pds", 100, 4800, 5600), "B_DTF_M", "data_sweights_pds")
    # hist_data_npds_b = data_df.Histo1D((f"B_DTF_M_npds", f"B_DTF_M_npds", 100, 4800, 5600), "B_DTF_M", "data_sweights_npds")
    #
    # temp_data_pds = hist_data_pds_b.DrawCopy("E")
    # temp_data_npds = hist_data_npds_b.DrawCopy("SAME")
    #
    # temp_data_pds.SetMarkerColor(ROOT.kRed)
    # temp_data_npds.SetMarkerColor(ROOT.kBlue)
    #
    # temp_data_pds.GetXaxis().SetTitle(f"B_DTF_M [MeV]")
    # # temp_data_pds.GetXaxis().SetLimits(temp_data_pds.GetXaxis().GetXmin()/1000000,temp_data_pds.GetXaxis().GetXmax()/1000000)
    # # temp_data_npds.GetXaxis().SetLimits(temp_data_npds.GetXaxis().GetXmin()/1000000,temp_data_npds.GetXaxis().GetXmax()/1000000)
    #
    # legend = ROOT.TLegend(0.70, 0.60, 0.95, 0.95)
    # legend.AddEntry(temp_data_pds, f"Data_SWeighted_pds_2","lp")
    # legend.AddEntry(temp_data_npds, f"Data_SWeighted_pds_3","lp")
    # legend.Draw()
    #
    # save_png(c2, f"{savename}", f"{spec}_{mcid}_bm", rpflag = 0)
def bdt_test(run_name, mc_list, strat):
    for tuple in mc_list:

        spec = tuple[0]
        mcid = tuple[1]

        data_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/final_sample/{spec}.root", "DecayTreeTuple", "DATA")
        sw_data_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/sw_files/sw_{spec}_test_2_10.root", "SW_tree", "DATA")

        # if strat == "split":
        #     mc_tchain, mc_tchain_nTaT = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/{year}/final_sample/{mcid}_{spec}*.root", "DecayTreeTuple", "MC")
        #     for year in ["2016","2017","2018"]:
        #         spec = tuple[0]
        #         mcid = tuple[1]
        #         if mcid == "02":
        #             yield_name = f"{spec}_0203"
        #         else:
        #             yield_name = f"{spec}_{mcid}"
        #
        #
        #         data_tchain.AddFriend(sw_data_tchain, "SW_TREE")
        #
        #         data_rdf = RDF(data_tchain)
        #         mc_T_rdf = RDF(mc_tchain_T)
        #         mc_nTaT_rdf = RDF(mc_tchain_nTaT)
        #         # sw_df = RDF(sw_data_tchain)
        #
        #         data_rdf = data_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
        #                          .Define("M2_d1k", f"({M2_d1kst})") \
        #                          .Define("M2_d2k", f"({M2_d2kst})")
        #                          # .Define(f"{yield_name}_yield_sw", f"SW_TREE.{yield_name}_yield_sw")
        #
        #         mc_T_rdf = mc_T_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
        #                          .Define("M2_d1k", f"({M2_d1kst})") \
        #                          .Define("M2_d2k", f"({M2_d2kst})")
        #
        #         mc_nTaT_rdf = mc_nTaT_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
        #                                .Define("M2_d1k", f"({M2_d1kst})") \
        #                                .Define("M2_d2k", f"({M2_d2kst})")
        #
        #         if mcid == "04":
        #
        #             d1kcut  = "M2_d1k  < 10000000"
        #             d2kcut  = "M2_d2k  < 10000000"
        #             d1d2cut = "M2_d1d2 < 17000000"
        #
        #             data_rdf = data_rdf.Filter(d1kcut, "test")
        #             data_rdf = data_rdf.Filter(d1d2cut, "test2")
        #             data_rdf = data_rdf.Filter(d2kcut, "test3")
        #
        #         # mc_T_rdf = mc_T_rdf.Filter(d1kcut, "test")
        #         # mc_T_rdf = mc_T_rdf.Filter(d1d2cut, "test2")
        #         # mc_T_rdf = mc_T_rdf.Filter(d2kcut, "test2")
        #         #
        #         # mc_nTaT_rdf =  mc_nTaT_rdf.Filter(d1kcut, "test")
        #         # mc_nTaT_rdf =  mc_nTaT_rdf.Filter(d1d2cut, "test2")
        #         # mc_nTaT_rdf = mc_nTaT_rdf.Filter(d2kcut, "test2")
        #
        #         data_allCutsReport = data_rdf.Report()
        #         mc_allCutsReport = mc_T_rdf.Report()
        #         data_allCutsReport.Print()
        #         mc_allCutsReport.Print()
        #
        #         data_np = data_rdf.AsNumpy(columns = dalitz_columns)
        #         data_df = pandas.DataFrame(data_np)
        #
        #         mc_T_np = mc_T_rdf.AsNumpy(columns = dalitz_columns)
        #         mc_T_df = pandas.DataFrame(mc_T_np)
        #         mc_T_weights_og = numpy.ones(len(mc_T_df))
        #
        #         mc_nTaT_np = mc_nTaT_rdf.AsNumpy(columns = dalitz_columns)
        #         mc_nTaT_df = pandas.DataFrame(mc_nTaT_np)
        #         mc_nTaT_weights_og = numpy.ones(len(mc_nTaT_df))
        #
        #         sw_np = data_rdf.AsNumpy(columns=[f"{yield_name}_yield_sw"])
        #         sweights_df = pandas.DataFrame(sw_np)
        #         sweights_np = sweights_df[f"{yield_name}_yield_sw"].to_numpy()
        #
        #         draw_sw_test(data_df, sweights_np, mc_T_df, mc_nTaT_df, mc_T_weights_og, mc_nTaT_weights_og, spec, mcid, year, "rw_no_sw")
        #
        #         data_df_rw = data_df[dalitz_columns]
        #         mc_T_df_rw = mc_T_df[dalitz_columns]
        #         mc_nTaT_df_rw = mc_nTaT_df[dalitz_columns]
        #
        #         reweighter_base_T = reweight.GBReweighter(max_depth=3, gb_args={'subsample': 0.5})
        #         reweighter_T = reweight.FoldingReweighter(reweighter_base_T, n_folds = 2)
        #
        #         reweighter_base_nTaT = reweight.GBReweighter(max_depth=3, gb_args={'subsample': 0.5})
        #         reweighter_nTaT = reweight.FoldingReweighter(reweighter_base_nTaT, n_folds= 2)
        #
        #         reweighter_T = reweighter_T.fit(mc_T_df_rw, data_df_rw, target_weight = sweights_np)
        #         reweighter_nTaT = reweighter_nTaT.fit(mc_nTaT_df_rw, data_df_rw, target_weight = sweights_np)
        #
        #         # test_data_weights = reweighter_T.predict_weights(data_df, sweights_np)
        #
        #         weights_T = reweighter_T.predict_weights(mc_T_df_rw)
        #         weights_nTaT = reweighter_nTaT.predict_weights(mc_nTaT_df_rw)
        #         print(mc_T_weights_og.sum(), weights_T.sum())
        #         print(mc_nTaT_weights_og.sum(), weights_nTaT.sum())
        #         draw_sw_test(data_df, sweights_np, mc_T_df, mc_nTaT_df, weights_T, weights_nTaT, spec, mcid,  year, "rw_sw")
        #
        #         mc_T_df = ROOT.RDF.MakeNumpyDataFrame({"rw_weights":  weights_T})
        #         mc_nTaT_df = ROOT.RDF.MakeNumpyDataFrame({"rw_weights":   weights_nTaT})
        #
        #         opts = ROOT.ROOT.RDF.RSnapshotOptions()
        #         opts.fMode = "UPDATE"
        #         t_clist = mc_T_df.GetColumnNames()
        #         ntat_clist = mc_nTaT_df.GetColumnNames()
        #
        #         outputfile = f"Analysis_2021/dalitz_reweighting/rw_mc_files/{mcid}_{spec}_{year}.root"
        #         if os.path.exists(outputfile):
        #             os.remove(outputfile)
        #             print("First deleting old ", outputfile)
        #         else:
        #             print("making ", outputfile)
        #         rdfsnap = mc_T_df.Snapshot(f"rwTuple_T", outputfile, t_clist, opts)
        #         rdfsnap_pre_d = mc_nTaT_df.Snapshot("rwTuple_nTaT", outputfile, ntat_clist, opts)
        #         print(f"Done")
        #
        #
        #         # # # reweighter
        #         # # # # # predict method provides unbiased weights prediction for the whole sample
        #         # # # # # folding reweighter contains two reweighters, each is trained on one half of samples
        #         # # # # # # during predictions each reweighter predicts another half of samples not used in training
        #         # #
        #         # # # draw_distributions(original, target, folding_weights, target_weights, spec, mcid, "2Fold")
        #         # #
        #         # # # for i in mc_columns:
        #         # # #     # col_original = original_test.eval(i, engine='python')
        #         # # #     # col_target = target_test.eval(i, engine='python')
        #         # # #     col_og_f = original.eval(i, engine='python')
        #         # # #     col_tar_f = target.eval(i, engine='python')
        #         # # #
        #         # # #     # w_target = target_weights
        #         # # #     print('No reweight   KS:', ks_2samp_weighted(col_original, col_target,
        #         # # #                                                  weights1=original_weights_test, weights2=target_weights_test))
        #         # # #     # print('Bins reweight KS:', ks_2samp_weighted(col_original, col_target,
        #         # # #     #                                              weights1=bins_weights_test, weights2=w_target))
        #         # # #     # print('GB Reweight   KS:', ks_2samp_weighted(col_original, col_target,
        #         # # #     #                                              weights1=gb_weights_test, weights2=target_weights_test))
        #         # # #     print('Fold Reweight   KS:', ks_2samp_weighted(col_og_f , col_tar_f,
        #         # # #                                                     weights1=folding_weights , weights2=target_weights))
        #         # #
        #         # #
        #         # # # from sklearn.ensemble import GradientBoostingClassifier
        #         # # # from sklearn.model_selection import train_test_split
        #         # # # from sklearn.metrics import roc_auc_score
        #         # # #
        #         # # # data = numpy.concatenate([original_test, target_test])
        #         # # # labels = numpy.array([0] * len(original_test) + [1] * len(target_test))
        #         # #
        #         # # # weights = {}
        #         # # # weights['original'] = original_weights_test
        #         # # # # weights['bins'] = bins_weights_test
        #         # # # weights['gb_weights'] = gb_weights_test
        #         # # #
        #         # # # for name, new_weights in weights.items():
        #         # # #     W = numpy.concatenate([new_weights / new_weights.sum() * len(target_test), target_weights_test / target_weights_test.sum() * len(target_test)])
        #         # # #     Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(data, labels, W, random_state=42, train_size=0.50)
        #         # # #     clf = GradientBoostingClassifier(subsample=0.5, n_estimators=50).fit(Xtr, Ytr, sample_weight=Wtr)
        #         # # #
        #         # # #     print(name, roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts))
        #         # # #
        #         # # # #
        #         # # # data = numpy.concatenate([original, target])
        #         # # # labels = numpy.array([0] * len(original) + [1] * len(target))
        #         # # #
        #         # # # weights = {}
        #         # # # weights['original'] = original_weights
        #         # # # weights['2-folding'] = folding_weights
        #         # # #
        #         # # # for name, new_weights in weights.items():
        #         # # #     W = numpy.concatenate([new_weights / new_weights.sum() * len(target), target_weights / target_weights.sum() * len(target)])
        #         # # #     Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(data, labels, W, random_state=42, train_size=0.50)
        #         # # #     clf = GradientBoostingClassifier(subsample=0.5, n_estimators=50).fit(Xtr, Ytr, sample_weight=Wtr)
        #         # # #     print(name, roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts))
        if strat == "all":

            if spec == "Z_m_p" and mcid == "02":
                yield_name = f"{spec}_0203"
            elif spec == "Z_z_z" and (mcid == "07" or mcid == "10"):
                yield_name = f"{spec}_0710"
            elif spec == "Z_z_z" and (mcid == "04" or mcid == "08" or mcid == "12"):
                yield_name = f"{spec}_040812"
            else:
                yield_name = f"{spec}_{mcid}"

            mc_tchain_T, mc_tchain_nTaT = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/*/final_sample/{mcid}_{spec}*.root", "DecayTreeTuple", "MC")

            data_tchain.AddFriend(sw_data_tchain, "SW_TREE")

            data_rdf = RDF(data_tchain)
            mc_T_rdf = RDF(mc_tchain_T)
            mc_nTaT_rdf = RDF(mc_tchain_nTaT)

            data_rdf = data_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
                             .Define("M2_d1k", f"({M2_d1kst})") \
                             .Define("M2_d2k", f"({M2_d2kst})")

            mc_T_rdf = mc_T_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
                               .Define("M2_d1k", f"({M2_d1kst})") \
                               .Define("M2_d2k", f"({M2_d2kst})")

            mc_nTaT_rdf = mc_nTaT_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
                                     .Define("M2_d1k", f"({M2_d1kst})") \
                                     .Define("M2_d2k", f"({M2_d2kst})")

            d1kcut  = f"M2_d1d2  < {tuple[2][0]}"
            d2kcut  = f"M2_d1k  < {tuple[2][1]}"
            d1d2cut = f"M2_d2k < {tuple[2][2]}"


            for name in ["no_cut", "cut"]:
                if name == "cut":
                    data_rdf_next = data_rdf.Filter(d1kcut, "test") \
                                            .Filter(d1d2cut, "test2") \
                                            .Filter(d2kcut, "test3")

                    mc_T_rdf_next = mc_T_rdf.Filter(d1kcut, "test") \
                                            .Filter(d1d2cut, "test2") \
                                            .Filter(d2kcut, "test3")

                    mc_nTaT_rdf_next = mc_nTaT_rdf.Filter(d1kcut, "test") \
                                                  .Filter(d1d2cut, "test2") \
                                                  .Filter(d2kcut, "test3")
                else:
                    data_rdf_next = data_rdf
                    mc_T_rdf_next = mc_T_rdf
                    mc_nTaT_rdf_next = mc_nTaT_rdf

                data_np = data_rdf_next.AsNumpy(columns = dalitz_columns)
                data_df = pandas.DataFrame(data_np)

                mc_T_np = mc_T_rdf_next.AsNumpy(columns = dalitz_columns)
                mc_T_df = pandas.DataFrame(mc_T_np)

                mc_nTaT_np = mc_nTaT_rdf.AsNumpy(columns = dalitz_columns)
                mc_nTaT_df = pandas.DataFrame(mc_nTaT_np)

                frames = [mc_T_df, mc_nTaT_df]
                mc_df = pandas.concat(frames, keys = ["T, nTaT"])
                mc_weights_og = numpy.ones(len(mc_df))

                sw_np = data_rdf_next.AsNumpy(columns=[f"{yield_name}_yield_sw"])
                sweights_df = pandas.DataFrame(sw_np)
                sweights_np = sweights_df[f"{yield_name}_yield_sw"].to_numpy()
                year = "all"

                draw_sw_test(data_df, sweights_np, mc_df, mc_weights_og, spec, mcid, year, f"rw_sw_{name}")

            data_df_rw = data_df[dalitz_columns]
            mc_df_rw = mc_df[dalitz_columns]

            reweighter_base = reweight.GBReweighter(max_depth=3, gb_args={'subsample': 0.5})
            reweighter = reweight.FoldingReweighter(reweighter_base, n_folds = 2)
            reweighter = reweighter.fit(mc_df_rw, data_df_rw, target_weight = sweights_np)

            # test_data_weights = reweighter_T.predict_weights(data_df, sweights_np)

            weights = reweighter.predict_weights(mc_df_rw)
            print(mc_weights_og.sum(),weights.sum())
            print(mc_weights_og.sum()/weights.sum())

            draw_sw_test(data_df, sweights_np, mc_df, weights, spec, mcid, year, f"rw_sw")

            with open(f"Analysis_2021/dalitz_reweighting/sys_txt_files/{mcid}_{spec}_{year}.txt", 'w') as f:
                f.write(f"{weights.sum()/mc_weights_og.sum()}")

            # mc_df = ROOT.RDF.MakeNumpyDataFrame({"rw_weights":  weights})

            # opts = ROOT.ROOT.RDF.RSnapshotOptions()
            # opts.fMode = "UPDATE"
            # clist = mc_df.GetColumnNames()
            #
            # outputfile = f"Analysis_2021/dalitz_reweighting/rw_mc_files/{mcid}_{spec}_{year}.root"
            # if os.path.exists(outputfile):
            #     os.remove(outputfile)
            #     print("First deleting old ", outputfile)
            # else:
            #     print("making ", outputfile)
            # rdfsnap = mc_df.Snapshot(f"rwTuple_T", outputfile, clist)
            # print(f"Done")
