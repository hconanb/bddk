import numpy
import root_numpy
import pandas
from hep_ml import reweight
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score

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

PE_d2h1d2h2 = "(pow(D2H1_PE + D2H2_PE,2))"
PX_d1d2 = "(pow(D1_PX + D2_PX,2))"
PY_d1d2 = "(pow(D1_PY + D2_PY,2))"
PZ_d1d2 = "(pow(D1_PZ + D2_PZ,2))"
M2_d1d2 = f"{PE_d1d2} - ({PX_d1d2} + {PY_d1d2} + {PZ_d1d2})"
M_d1d2 = f"sqrt({M2_d1d2})"

#########################################MCDecayTreeTuple
# PE_d1d2_mcdtt = "(pow(C1_TRUEP_E + C2_TRUEP_E,2))"
# PX_d1d2_mcdtt = "(pow(C1_TRUEP_X + C2_TRUEP_X,2))"
# PY_d1d2_mcdtt = "(pow(C1_TRUEP_Y + C2_TRUEP_Y,2))"
# PZ_d1d2_mcdtt = "(pow(C1_TRUEP_Z + C2_TRUEP_Z,2))"
# M2_d1d2_mcdtt = f"{PE_d1d2_mcdtt} - ({PX_d1d2_mcdtt} + {PY_d1d2_mcdtt} + {PZ_d1d2_mcdtt})"
# M_d1d2_mcdtt = f"sqrt({M2_d1d2_mcdtt})"
#
# PE_d1kst_mcdtt = "(pow(C1_TRUEP_E + Kst_TRUEP_E,2))"
# PX_d1kst_mcdtt = "(pow(C1_TRUEP_X + Kst_TRUEP_X,2))"
# PY_d1kst_mcdtt = "(pow(C1_TRUEP_Y + Kst_TRUEP_Y,2))"
# PZ_d1kst_mcdtt = "(pow(C1_TRUEP_Z + Kst_TRUEP_Z,2))"
# M2_d1kst_mcdtt = f"{PE_d1kst_mcdtt} - ({PX_d1kst_mcdtt} + {PY_d1kst_mcdtt} + {PZ_d1kst_mcdtt})"
# M_d1kst_mcdtt = f"sqrt({M2_d1kst_mcdtt})"
#
# PE_d2kst_mcdtt = "(pow(C2_TRUEP_E + Kst_TRUEP_E,2))"
# PX_d2kst_mcdtt = "(pow(C2_TRUEP_X + Kst_TRUEP_X,2))"
# PY_d2kst_mcdtt = "(pow(C2_TRUEP_Y + Kst_TRUEP_Y,2))"
# PZ_d2kst_mcdtt = "(pow(C2_TRUEP_Z + Kst_TRUEP_Z,2))"
# M2_d2kst_mcdtt = f"{PE_d2kst_mcdtt} - ({PX_d2kst_mcdtt} + {PY_d2kst_mcdtt} + {PZ_d2kst_mcdtt})"
# M_d2kst_mcdtt = f"sqrt({M2_d2kst_mcdtt})"

M2_d1d2_norm = M2_d1d2.replace("KST","K")
M2_d1kst_norm = M2_d1kst.replace("KST","K")
M2_d2kst_norm = M2_d2kst.replace("KST","K")

dalitz_columns = ["M2_d1d2", "M2_d2k", "M2_d1k"]
# dalitz_columns = ["M2_d1d2", "M2_d1k"]

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
    if type == "DATA":
        tchain = ROOT.TChain(chain_name)
        for file_name in file_list:
            tchain.Add(file_name)
        return tchain
    if type == "MC" or type == "norm":
        tchain = ROOT.TChain()
        for file_name in file_list:
            tchain.Add(f"{file_name}/DecayTreeTuple_T")
            tchain.Add(f"{file_name}/DecayTreeTuple_nTaT")
        return tchain
    if type == "Gen":
        tchain = ROOT.TChain()
        for file_name in file_list:
            tchain.Add(f"{file_name}/DecayTree")
        return tchain

def draw_sw_test(dt, rect, gent, spec, mcid, year, savename):
    c1 = ROOT.TCanvas("c1","c1")
    for column in dalitz_columns:
        bins = 50

        DATA = dt[0]
        data_weights = dt[1]
        REC = rect[0]
        mc_weights = rect[1]
        GEN = gent[0]
        gen_weights = gent[1]

        data_np = DATA[column].to_numpy()
        rec_np = REC[column].to_numpy()
        gen_np = GEN[column].to_numpy()

        data_df = ROOT.RDF.MakeNumpyDataFrame({column:data_np, "weights": data_weights})
        rec_df = ROOT.RDF.MakeNumpyDataFrame({column:rec_np, "weights":  mc_weights})
        gen_df = ROOT.RDF.MakeNumpyDataFrame({column:gen_np, "weights": gen_weights})

        xmin = data_df.Min(column).GetValue()
        xmax = data_df.Max(column).GetValue()

        hist_data = data_df.Histo1D((f"data_{column}_sw", f"data_{column}_sw", bins, xmin, xmax), column, "weights")
        hist_rec = rec_df.Histo1D((f"rec_{column}", f"rec_{column}", bins, xmin, xmax), column, "weights")
        hist_gen = gen_df.Histo1D((f"gen_{column}", f"gen_{column}", bins, xmin, xmax), column, "weights")

        temp_data = hist_data.DrawCopy("E")
        temp_rec = hist_rec.DrawCopy("E SAME")
        temp_gen = hist_gen.DrawCopy("E SAME")

        temp_data.SetMarkerColor(ROOT.kRed)
        temp_data.Scale(1/temp_data.Integral())

        temp_rec.SetMarkerColor(ROOT.kBlue)
        temp_rec.Scale(1/temp_rec.Integral())

        temp_gen.SetMarkerColor(ROOT.kGreen)
        temp_gen.Scale(1/temp_gen.Integral())

        plot_min = (temp_data.GetXaxis().GetXmin()/1000000)
        plot_max = (temp_data.GetXaxis().GetXmax()/1000000)

        temp_data.GetXaxis().SetTitle(f"m^{{2}}_{{({dalitz_columns_dict[column]})}} [GeV^{{2}}]")
        # temp_data.GetXaxis().SetLimits(plot_min, plot_max)
        temp_data.GetXaxis().SetLimits(plot_min, plot_max)
        temp_rec.GetXaxis().SetLimits(plot_min, plot_max)
        temp_gen.GetXaxis().SetLimits(plot_min, plot_max)
        leg_f = "tl"
        print(column)
        if column == "M2_d1d2":
            leg_f = "tl"
        if column == "M2_d1k":
            leg_f = "tr"
        if column == "M2_d2k":
            leg_f = "tr"


        if leg_f == "tl":
            legend = ROOT.TLegend(0.15, 0.65, 0.40, 0.95)
        if leg_f == "tr":
            legend = ROOT.TLegend(0.70, 0.65, 0.95, 0.90)

        legend.AddEntry(temp_data, f"Data_SWeighted","lp")
        legend.AddEntry(temp_rec, f"MC: Final Sample","lp")
        legend.AddEntry(temp_gen, f"MC: Generator Level ","lp")

        legend.SetTextSize(0.0300)
        legend.SetFillStyle(1001)

        legend.Draw()

        save_png(c1, f"{savename}", f"{spec}_{mcid}_{column}_{year}_{savename}", rpflag = 0)

def draw_dalitz(tup, tuname, spec, mcid, max_dd, max_dk, year, savename):

    events_np = tup[0].to_numpy()
    weights = tup[1]

    data_rdf = ROOT.RDF.MakeNumpyDataFrame({dalitz_columns[0]:events_np[:,0], dalitz_columns[1]:events_np[:,1], dalitz_columns[2]:events_np[:,2], "weights": weights})

    dbins = 20
    min_dk = 6.5e6
    min_dd = 13.5e6

    data_hdalitz_1k_12 = data_rdf.Histo2D(("", f"{tuname}_M2_d1k_M2_d1d2", dbins, min_dk, max_dk, dbins, min_dd, max_dd), "M2_d1k", "M2_d1d2", "weights")
    data_hdalitz_2k_12 = data_rdf.Histo2D(("", f"{tuname}_M2_d2k_M2_d1d2", dbins, min_dk, max_dk, dbins, min_dd, max_dd), "M2_d2k", "M2_d1d2", "weights")
    data_hdalitz_1k_2k = data_rdf.Histo2D(("", f"{tuname}_M2_d1k_M2_d2k", dbins, min_dk, max_dk, dbins, min_dk, max_dk), "M2_d1k", "M2_d2k", "weights")

    for histo in [data_hdalitz_1k_12, data_hdalitz_2k_12, data_hdalitz_1k_2k]:

        tltxt = ROOT.TLatex()
        pname = histo.GetTitle()

        type = pname.split("_")[0]
        xnt = pname.split("_")[2]
        ynt = pname.split("_")[4]

        if spec == "Z_m_p":
            if xnt == "d1k":
                xname = "D^{-}K^{*0}"
            if xnt == "d2k":
                xname = "D^{+}K^{*0}"
            if xnt == "d1d2":
                xname = "D^{-}D^{+}"
            if ynt == "d1k":
                yname = "D^{-}K^{*0}"
            if ynt == "d2k":
                yname = "D^{+}K^{*0}"
            if ynt == "d1d2":
                yname = "D^{-}D^{+}"

        histo.GetXaxis().SetTitle(f"m^{{2}}({xname}) [GeV^{{2}}]")
        histo.GetXaxis().SetLimits(histo.GetXaxis().GetXmin()/1000000,histo.GetXaxis().GetXmax()/1000000)
        histo.GetYaxis().SetTitle(f"m^{{2}}({yname}) [GeV^{{2}}]")
        histo.GetYaxis().SetLimits(histo.GetYaxis().GetXmin()/1000000,histo.GetYaxis().GetXmax()/1000000)
        # histo.GetZaxis().SetTitle(f"m^{{2}}({yname}) [GeV^{{2}}]")

        c1 = ROOT.TCanvas("c1","c1")
        c1.SetRightMargin(0.15)
        c1.SetLeftMargin(0.9)
        c1.SetTopMargin(0.1)
        c1.SetBottomMargin(0.9)
        histo.Draw("COLZ1")
        tltxt.DrawLatexNDC(.3,.925, f"Dalitz Distribution: {tuname}")
        dfa = savename.split("_")[1]
        save_png(c1, f"{savename}", f"{pname}_{type}_{spec}_{mcid}_{dfa}", rpflag = 0)

def bdt_test(run_name, mc_list):

    for tuple in mc_list:

        spec = tuple[0]
        mcid = tuple[1]

        if spec == "norm7" or spec == "norm8":
            data_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{spec}.root", "DecayTreeTuple", "norm")
            sw_data_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/normalization_fit/sw_files/sw_{spec}.root", "SW_tree", "DATA")
        else:
            data_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{spec}.root", "DecayTreeTuple_SIG", "DATA")
            sw_data_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/sw_files/sw_{spec}_{run_name}.root", "SW_tree", "DATA")

        if spec == "Z_m_p" and mcid == "02":
            yield_name = f"{spec}_0203_yield_sw"
        # elif spec == "Z_z_z" and (mcid == "07" or mcid == "10"):
        #     yield_name = f"{spec}_0710_yield_sw"
        # elif spec == "Z_z_z" and (mcid == "04" or mcid == "08" or mcid == "12"):
        #     yield_name = f"{spec}_040812_yield_sw"
        # elif spec == "norm7" or spec == "norm8":
        #     yield_name = f"n_{spec}_signal_sw"
        else:
            yield_name = f"{spec}_{mcid}_yield_sw"

        mc_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/*/final_sample/{mcid}_{spec}*.root", "DecayTreeTuple", "MC")
        data_tchain.AddFriend(sw_data_tchain, "SW_TREE")
        # gen_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/B2DplusDminusKstarzero_tree_BIG.root", "DecayTree", "Gen")
        gen_tchain = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_redecay/{mcid}_{spec}*.root", "DecayTree", "Gen")

        data_rdf = RDF(data_tchain)
        mc_rdf = RDF(mc_tchain)
        gen_rdf = RDF(gen_tchain)

        if "norm" not in spec:
            data_rdf = data_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
                             .Define("M2_d1k", f"({M2_d1kst})") \
                             .Define("M2_d2k", f"({M2_d2kst})")

            mc_rdf = mc_rdf.Define("M2_d1d2", f"({M2_d1d2})") \
                               .Define("M2_d1k", f"({M2_d1kst})") \
                               .Define("M2_d2k", f"({M2_d2kst})")
            if mcid == "02":
                gen_rdf = gen_rdf.Define("M2_d1d2", "Dstp_0_Dm_0_M2*1000000") \
                                 .Define("M2_d1k", "Dm_0_Kst0_0_M2*1000000") \
                                 .Define("M2_d2k", "Dstp_0_Kst0_0_M2*1000000")
            if mcid == "01":
                gen_rdf = gen_rdf.Define("M2_d1d2", "Dp_0_Dm_0_M2*1000000") \
                                 .Define("M2_d1k", "Dm_0_Kst0_0_M2*1000000") \
                                 .Define("M2_d2k", "Dp_0_Kst0_0_M2*1000000")

        if "norm" in spec:
            data_rdf = data_rdf.Define("M2_d1d2", f"({M2_d1d2_norm})") \
                             .Define("M2_d1k", f"({M2_d1kst_norm})") \
                             .Define("M2_d2k", f"({M2_d2kst_norm})")

            mc_rdf = mc_rdf.Define("M2_d1d2", f"({M2_d1d2_norm})") \
                               .Define("M2_d1k", f"({M2_d1kst_norm})") \
                               .Define("M2_d2k", f"({M2_d2kst_norm})")

        d1kcut  = f"M2_d1d2  < {tuple[2][0]}"
        d2kcut  = f"M2_d1k  < {tuple[2][1]}"
        d1d2cut = f"M2_d2k < {tuple[2][2]}"

        data_np_base = data_rdf.AsNumpy(columns = dalitz_columns)
        data_df_base = pandas.DataFrame(data_np_base)
        data_weights_no_sw = numpy.ones(len(data_df_base))

        sw_np_base = data_rdf.AsNumpy(columns = [yield_name])
        sw_df_base = pandas.DataFrame(sw_np_base)
        data_weights_sw = sw_df_base.to_numpy()

        mc_np_base = mc_rdf.AsNumpy(columns = dalitz_columns)
        mc_df_base = pandas.DataFrame(mc_np_base)
        mc_weights_og = numpy.ones(len(mc_df_base))

        gen_np_base = gen_rdf.AsNumpy(columns = dalitz_columns)
        gen_df_base = pandas.DataFrame(gen_np_base)
        gen_weights_og = numpy.ones(len(gen_df_base))

        data_t = (data_df_base, data_weights_sw)
        mc_rec_t = (mc_df_base, mc_weights_og)
        mc_gen_t = (gen_df_base, gen_weights_og)

        # draw_sw_test(data_t, mc_rec_t, mc_gen_t, spec, mcid, "all", f"prerw_all_nocut")
        # draw_dalitz(data_df_base, sw_df_base, mc_df_base, mc_weights_og, spec, mcid, tuple[2][0], tuple[2][1], "all", f"dalitz_prerw_all_nocut")

        data_rdf_next = data_rdf.Filter(d1kcut, "test") \
                                .Filter(d1d2cut, "test2") \
                                .Filter(d2kcut, "test3")

        mc_rdf_next = mc_rdf.Filter(d1kcut, "test") \
                            .Filter(d1d2cut, "test2") \
                            .Filter(d2kcut, "test3")

        gen_rdf_next = gen_rdf.Filter(d1kcut, "test") \
                            .Filter(d1d2cut, "test2") \
                            .Filter(d2kcut, "test3") \
                            .Filter("abs(Kst0_0_M - 0.895) < 0.050", "Kst_Cut")

        #
        data_np = data_rdf_next.AsNumpy(columns = dalitz_columns)
        data_df = pandas.DataFrame(data_np)

        sw_np = data_rdf_next.AsNumpy(columns = [yield_name])
        sw_df = pandas.DataFrame(sw_np)
        data_weights = sw_df[f"{yield_name}"].to_numpy()

        mc_np = mc_rdf_next.AsNumpy(columns = dalitz_columns)
        mc_df = pandas.DataFrame(mc_np)
        mc_weights = numpy.ones(len(mc_df))

        gen_np = gen_rdf_next.AsNumpy(columns = dalitz_columns)
        gen_df = pandas.DataFrame(gen_np)
        gen_weights = numpy.ones(len(gen_df))

        data_t2 = (data_df, data_weights)
        mc_rec_t2 = (mc_df, mc_weights)
        mc_gen_t2 = (gen_df, gen_weights)

        draw_sw_test(data_t2, mc_rec_t2, mc_gen_t2, spec, mcid, "all", f"prerw_wcut")
        for tu, tuname in zip([data_t2, mc_rec_t2, mc_gen_t2], ["Data", "Rec", "Gen"]):
            draw_dalitz(tu, tuname, spec, mcid, tuple[2][0], tuple[2][1], "all", f"dalitz_prerw_wcut")


        reweighter_base_kfold = reweight.GBReweighter(max_depth = 2, gb_args={'subsample': 0.5})
        reweighter_kfold = reweight.FoldingReweighter(reweighter_base_kfold, n_folds = 2, random_state = 42)
        reweighter_kfold = reweighter_kfold.fit(original=mc_df, target=data_df[dalitz_columns], target_weight = data_weights)

        # # genweighter_base_kfold = reweight.GBReweighter(max_depth = 2, gb_args={'subsample': 0.5})
        # genweighter_kfold = reweight.FoldingReweighter(genweighter_base_kfold, n_folds = 2)
        # gen_reweighter_kfold = genweighter_kfold.fit(original=gen_df, target=data_df[dalitz_columns], target_weight = sweights_np)
        #
        mc_weights_kfold = reweighter_kfold.predict_weights(mc_df)
        gen_weights_kfold = reweighter_kfold.predict_weights(gen_df)

        mc_rec_t3 = (mc_df, mc_weights_kfold)
        mc_gen_t3 = (gen_df, gen_weights_kfold)

        draw_sw_test(data_t2, mc_rec_t3, mc_gen_t3, spec, mcid, "all", f"rw_wcut")
        for tu, tuname in zip([mc_rec_t3, mc_gen_t3], ["Rec", "Gen"]):
            draw_dalitz(tu, tuname, spec, mcid, tuple[2][0], tuple[2][1], "all", f"dalitz_rw_wcut")

        print(f"OG {mc_weights.sum()} NEW {mc_weights_kfold.sum()} Change {(mc_weights.sum() - mc_weights_kfold.sum())/mc_weights.sum()}")
        print(f"OG {gen_weights.sum()} NEW {gen_weights_kfold.sum()} Change {(gen_weights.sum() - gen_weights_kfold.sum())/gen_weights.sum()}")

        oe = mc_weights.sum()/gen_weights.sum()
        ne = mc_weights_kfold.sum()/gen_weights_kfold.sum()

        print(f"OG Eff {oe} New Eff {ne} Change {(oe - ne)/oe}")

        with open(f"Analysis_2021/dalitz_reweighting/sys_txt_files/{mcid}_{spec}_2022.txt", 'w') as f:
           f.write(f"OG {mc_weights.sum()} NEW {mc_weights_kfold.sum()} Change {(mc_weights.sum() - mc_weights_kfold.sum())/mc_weights.sum()}\n")
           f.write(f"OG {gen_weights.sum()} NEW {gen_weights_kfold.sum()} Change {(gen_weights.sum() - gen_weights_kfold.sum())/gen_weights.sum()}\n")
           f.write(f"OG Eff {oe} New Eff {ne} Change {(oe - ne)/oe}")

        # total = numpy.concatenate([mc_df, data_df[dalitz_columns]])
        # labels = numpy.array([0] * len(mc_df) + [1] * len(data_df[dalitz_columns]))

        # weights = {}
        # weights['original'] = mc_weights
        # weights['2-folding'] = mc_weights_kfold
        # og_w_auc = 0
        # new_w_auc = 0
        #
        # for name, new_weights in weights.items():
        #     # / new_weights.sum() * len(data_df[dalitz_columns]),
        #     # # /sweights_np.sum() * len(data_df[dalitz_columns])
        #     # print(name)
        #     # print(new_weights)
        #     # print(new_weights.sum())
        #     # print(new_weights/new_weights.sum())
        #     # print(len(data_df[dalitz_columns]))
        #     # print(new_weights/new_weights.sum()*len(data_df[dalitz_columns]))
        #
        #     # W = numpy.concatenate([(new_weights/new_weights.sum())*len(data_df[dalitz_columns]), [1]*len(data_df[dalitz_columns])])
        #
        #     W = numpy.concatenate([(new_weights/new_weights.sum())*sweights_np.sum(), sweights_np])
        #     Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(total, labels, W, random_state=42, train_size=0.51)
        #     clf = GradientBoostingClassifier(subsample=0.3, n_estimators=30).fit(Xtr, Ytr, sample_weight=Wtr)
        #
        #     if name == 'original':
        #         og_w_auc = roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts)
        #     if name == '2-folding':
        #         new_w_auc = roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts)
        #
        #     print(name, roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts))
        #

        # mc_train, mc_test = train_test_split(mc_df)
        # data_train_w_sw, data_test_w_sw = train_test_split(data_df)

        # data_train = data_train_w_sw[dalitz_columns]
        # data_test = data_test_w_sw[dalitz_columns]
        #
        # sw_train = data_train_w_sw[yield_name]
        # sw_test = data_test_w_sw[yield_name]
        #
        # mc_weights_train = numpy.ones(len(mc_train))
        # mc_weights_test = numpy.ones(len(mc_test))
        # mc_weights_og = numpy.ones(len(mc_df))

        # draw_sw_test(data_test, sw_test, mc_test, mc_weights_test, spec, mcid, "all", f"prerw_test_wcut")

        # reweighter = reweight.GBReweighter(n_estimators=50, learning_rate=0.1, max_depth=3, min_samples_leaf=1000,
        #                                    gb_args={'subsample': 0.5})
        #
        # reweighter.fit(mc_train, data_train, target_weight=sw_train)
        #
        # mc_gbweights_test = reweighter.predict_weights(mc_test)

        # draw_sw_test(data_test, sw_test, mc_test, mc_gbweights_test, spec, mcid, "all", f"rw_gb_wcut")

def bdt_comp(run_name):

        spec = "Z_m_p"
        mcid = "01"

        yield_name = f"{spec}_{mcid}_yield_sw"
        run_name = "sw_test_2_18"

        data_tchain_all = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/final_sample/{spec}.root", "DecayTreeTuple", "DATA")
        sw_data_tchain_all = grab_files_and_chain(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/sw_files/sw_{spec}_{run_name}.root", "SW_tree", "DATA")

        data_tchain_r = grab_files_and_chain(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/final_sample/{spec}.root", "DecayTreeTuple", "DATA")
        sw_data_tchain_r = grab_files_and_chain(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/sw_files/sw_{spec}_swcomp.root", "SW_tree", "DATA")

        data_tchain_all.AddFriend(sw_data_tchain_all, "SW_TREE")
        data_tchain_r.AddFriend(sw_data_tchain_r, "SW_TREE")

        data_rdf_all= RDF(data_tchain_all)
        data_rdf_r = RDF(data_tchain_r)

        data_rdf_all = data_rdf_all.Define("M2_d1d2", f"({M2_d1d2})") \
                         .Define("M2_d1k", f"({M2_d1kst})") \
                         .Define("M2_d2k", f"({M2_d2kst})")

        # data_rdf_r = data_rdf_r.Filter("B_DTF_M < 5350 && B_DTF_M > 5200")
        #
        # data_rdf_r = data_rdf_r.Define("M2_d1d2", f"({M2_d1d2})") \
        #                        .Define("M2_d1k", f"({M2_d1kst})") \
        #                        .Define("M2_d2k", f"({M2_d2kst})")


        # hprof2d = data_rdf_all.Profile1D(("hprof1d", "Profile of pz versus px", 100, -0.5, 1.5), f"{yield_name}", "B_DTF_M")
        # hprof2d_swap = data_rdf_all.Profile1D(("hprof1ds", "Profile of pz versus pxs", 100, 4800,5600),"B_DTF_M",f"{yield_name}")

        # hprof2d.SetMinimum(4800)
        # hprof2d.SetMaximum(5600)
        #
        # hprof2d_swap.SetMinimum(-0.5)
        # hprof2d_swap.SetMaximum(1.5)
        #
        # c1 = ROOT.TCanvas("c1","c1")
        # hprof2d.Draw()
        # c1.SaveAs("test.png")
        # c2 = ROOT.TCanvas("c2","c2")
        # hprof2d_swap.Draw()
        # c2.SaveAs("test2.png")
        # data_np_base_all = data_rdf_all.AsNumpy(columns = dalitz_columns)
        # data_df_base_all = pandas.DataFrame(data_np_base_all)

        # sw_np_all = data_rdf_all.AsNumpy(columns=[f"{yield_name}"])
        # sweights_df_all = pandas.DataFrame(sw_np_all)
        # sweights_np_all = sweights_df_all[f"{yield_name}"].to_numpy()
        #
        # data_np_base_r = data_rdf_r.AsNumpy(columns = dalitz_columns)
        # data_df_base_r = pandas.DataFrame(data_np_base_r)
        #
        # sw_np_r = data_rdf_r.AsNumpy(columns=[f"{yield_name}"])
        # sweights_df_r = pandas.DataFrame(sw_np_r)
        # sweights_np_r = sweights_df_r[f"{yield_name}"].to_numpy()


        # draw_sw_test(data_df_base_all, sweights_np_all, data_df_base_r,sweights_np_r, spec, mcid, "all", f"prerw_comp")
