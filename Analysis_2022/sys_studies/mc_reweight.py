import sys
import os

basedir = os.getcwd().split('sys_studies')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

from hep_ml import reweight
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from hep_ml.metrics_utils import ks_2samp_weighted
import matplotlib.pyplot as plt
import warnings



seed = 7
np.random.seed(seed)

def get_dalitz_string(plist):
    """arguments:
    plist -- list of particles, e.g., ["K1", "K2", ...]
    """
    sflag = 0
    for p in plist:
        if sflag == 0:
            e = f"{p}_PE"
            x = f"{p}_PX"
            y = f"{p}_PY"
            z = f"{p}_PZ"
        if sflag == 1:
            e = f"{e}+{p}_PE"
            x = f"{x}+{p}_PX"
            y = f"{y}+{p}_PY"
            z = f"{z}+{p}_PZ"
        sflag = 1
    e2 = f"pow({e},2)"
    x2 = f"pow({x},2)"
    y2 = f"pow({y},2)"
    z2 = f"pow({z},2)"
    s = f"({e2}-{x2}-{y2}-{z2})"
    return s
def draw_1d_dalitz(data_rdf, mc_rdf, gen_rdf, data_spec, mc_spec, dalitz_columns_dict, w_flag, save_name):

    c1 = ROOT.TCanvas("c1","c1")

    for column in dalitz_columns_dict:

        bins = 50

        xmin = data_rdf.Min(column).GetValue()
        xmax = data_rdf.Max(column).GetValue()

        print(xmin, xmax)

        if w_flag == "og":

            hist_data = data_rdf.Histo1D((f"data_{column}_sw", f"data_{column}_sw", bins, xmin, xmax), column, f"{sw_id}")
            hist_rec = mc_rdf.Histo1D((f"rec_{column}", f"rec_{column}", bins, xmin, xmax), column, "og_weights")
            hist_gen = gen_rdf.Histo1D((f"gen_{column}", f"gen_{column}", bins, xmin, xmax), column, "og_weights")

        if w_flag == "rw":
            hist_data = data_rdf.Histo1D((f"data_{column}_sw", f"data_{column}_sw", bins, xmin, xmax), column, f"{sw_id}")
            hist_rec = mc_rdf.Histo1D((f"rec_{column}", f"rec_{column}", bins, xmin, xmax), column, "new_weights")
            hist_gen = gen_rdf.Histo1D((f"gen_{column}", f"gen_{column}", bins, xmin, xmax), column, "new_weights")

        temp_data = hist_data.DrawCopy("E")
        temp_rec = hist_rec.DrawCopy("E SAME")
        temp_gen = hist_gen.DrawCopy("E SAME")
        #
        temp_data.SetMarkerColor(ROOT.kRed)
        temp_data.Scale(1/temp_data.Integral())

        temp_rec.SetMarkerColor(ROOT.kBlue)
        temp_rec.Scale(1/temp_rec.Integral())
        #
        temp_gen.SetMarkerColor(ROOT.kGreen)
        temp_gen.Scale(1/temp_gen.Integral())
        #
        plot_min = (temp_data.GetXaxis().GetXmin()/1000000)
        plot_max = (temp_data.GetXaxis().GetXmax()/1000000)

        temp_data.GetXaxis().SetTitle(f"m^{{2}}_{{({dalitz_columns_dict[column]})}} [GeV^{{2}}]")
        temp_data.GetXaxis().SetLimits(plot_min, plot_max)
        temp_rec.GetXaxis().SetLimits(plot_min, plot_max)
        temp_gen.GetXaxis().SetLimits(plot_min, plot_max)

        leg_f = "tl"

        legend = ROOT.TLegend(0.70, 0.65, 0.95, 0.90)

        legend.AddEntry(temp_data,f"Data_SWeighted","lp")
        legend.AddEntry(temp_rec, f"MC: Final Sample","lp")
        legend.AddEntry(temp_gen, f"MC: RapidSim ","lp")

        legend.SetTextSize(0.0300)
        legend.SetFillStyle(1001)

        legend.Draw()

        save_pdf(c1, f"{save_name}", f"{mc_spec}_{column}_1d")
def draw_2d_dalitz(rdf, data_spec, mc_spec, dalitz_columns_dict, save_name):
    dbins = 20
    min_dk = 6.5e6
    min_dd = 13.5e6
    max_dk = 13e6
    max_dd = 21e6

    if "data" in save_name:
        wname = sw_id
        dfname = data_spec
    else:
        wname = "og_weights"
        dfname = mc_spec

    hdalitz_1k_12 = rdf.Histo2D(("", f"M2_d1k_M2_d1d2", dbins, min_dk, max_dk, dbins, min_dd, max_dd), "M2_d1k", "M2_d1d2", f"{wname}")
    hdalitz_2k_12 = rdf.Histo2D(("", f"M2_d2k_M2_d1d2", dbins, min_dk, max_dk, dbins, min_dd, max_dd), "M2_d2k", "M2_d1d2", f"{wname}")
    hdalitz_1k_2k = rdf.Histo2D(("", f"M2_d1k_M2_d2k", dbins, min_dk, max_dk, dbins, min_dk, max_dk), "M2_d1k", "M2_d2k", f"{wname}")

    for histo in [hdalitz_1k_12, hdalitz_2k_12, hdalitz_1k_2k]:

        tltxt = ROOT.TLatex()
        pname = histo.GetTitle()

        xnt = pname.split("_")[1]
        ynt = pname.split("_")[3]

        if data_spec == "Z_m_p":
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
        tltxt.DrawLatexNDC(.3,.925, f"Dalitz Distribution: {dfname}")

        save_pdf(c1, f"{save_name}", f"{dfname}_{pname}")

def optimise(data_rdf, mc_rdf, gen_rdf, dalitz_columns_dict, data_spec, mc_spec):


    data_np = data_rdf.AsNumpy(columns = dalitz_columns_dict.keys())
    data_df = pd.DataFrame(data_np)

    sw_np = data_rdf.AsNumpy(columns = [sw_id])
    sw_df = pd.DataFrame(sw_np)
    data_weights = sw_df[f"{sw_id}"].to_numpy()

    mc_np = mc_rdf.AsNumpy(columns = dalitz_columns_dict.keys())
    mc_df = pd.DataFrame(mc_np)
    mc_og_sum = mc_rdf.Count().GetValue()

    gen_np = gen_rdf.AsNumpy(columns = dalitz_columns_dict.keys())
    gen_df = pd.DataFrame(gen_np)
    gen_og_sum = gen_rdf.Count().GetValue()

    aucs = []
    fprs = []
    tprs = []
    params = []

    n_estimators = [65, 75, 85, 95, 105]
    depths = [2, 3, 4, 5, 6]

    for n_est in n_estimators:
        for depth in depths:
            print(f'Estimators: {n_est}, Depth: {depth}')

            reweighter_base = reweight.GBReweighter(n_estimators=n_est, max_depth=depth, gb_args={'subsample': 0.5})
            reweighter = reweight.FoldingReweighter(reweighter_base, n_folds=10)
            reweighter.fit(original=mc_df[dalitz_columns_dict.keys()], target=data_df[dalitz_columns_dict.keys()], target_weight=data_weights)

            MC_weights = reweighter.predict_weights(mc_df[dalitz_columns_dict.keys()], vote_function =  lambda x: np.mean(x, axis=0))


            # checks = [True if x < 1 else False for x in MC_weights]

            # if True in checks:
            #     print checks
            #     print 'Not all 1s'
            #     print MC_weights
            # else:
            #     print 'All 1s'

            data_full = np.concatenate([mc_df, data_df])
            labels = np.array([0] * len(mc_df) + [1] * len(data_df))

            W = np.concatenate([MC_weights / MC_weights.sum() * len(mc_df), data_weights / data_weights.sum() * len(data_df)])
            Xtr, Xts, Ytr, Yts, Wtr, Wts = train_test_split(data_full, labels, W, train_size=0.51, shuffle = True, random_state = seed)
            clf = GradientBoostingClassifier(subsample = 0.5, n_estimators = n_est, random_state = seed).fit(Xtr, Ytr, sample_weight=Wtr)

            auc = roc_auc_score(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts)
            fpr, tpr, __ = roc_curve(Yts, clf.predict_proba(Xts)[:, 1], sample_weight=Wts)

            aucs.append(auc)
            fprs.append(fpr)
            tprs.append(tpr)
            params.append((n_est, depth))
            print('ROC AUC = {:.3f}'.format(auc))
            print(MC_weights.sum(), mc_og_sum)

    closest = abs(np.subtract(np.array(aucs), 0.5))
    top5 = sorted(range(len(aucs)), key=lambda i: closest[i])[:5]

    fig = plt.figure()

    for i in top5:
        plt.plot(fprs[i], tprs[i], label = 'Estimators: {}, Depth: {}, AUC: {:.3f}'.format(params[i][0], params[i][1], aucs[i]))

    plt.plot([0,1], [0,1], 'k--')
    plt.legend()
    plt.savefig(f'plots/GBRweighting_OPT_ROC.png')
def mc_rw(data_rdf, mc_rdf, gen_rdf, dalitz_columns_dict, data_spec, mc_spec):

    data_np = data_rdf.AsNumpy(columns = dalitz_columns_dict.keys())
    data_df = pd.DataFrame(data_np)

    sw_np = data_rdf.AsNumpy(columns = [sw_id])
    sw_df = pd.DataFrame(sw_np)
    data_weights = sw_df[f"{sw_id}"].to_numpy()

    mc_np = mc_rdf.AsNumpy(columns = dalitz_columns_dict.keys())
    mc_df = pd.DataFrame(mc_np)
    mc_og_sum = mc_rdf.Count().GetValue()

    gen_np = gen_rdf.AsNumpy(columns = dalitz_columns_dict.keys())
    gen_df = pd.DataFrame(gen_np)
    gen_og_sum = gen_rdf.Count().GetValue()

    aucs = []
    fprs = []
    tprs = []
    params = []

    if "norm" not in data_spec:
        n_est = 75
        depth = 5
    if "norm" in data_spec:
        n_est = 75
        depth = 3

    print(f'Estimators: {n_est}, Depth: {depth}')

    reweighter_base = reweight.GBReweighter(n_estimators=n_est, max_depth=depth, gb_args={'subsample': 0.5})
    reweighter = reweight.FoldingReweighter(reweighter_base, n_folds=10)
    reweighter.fit(original=mc_df[dalitz_columns_dict.keys()], target=data_df[dalitz_columns_dict.keys()], target_weight=data_weights)

    MC_weights = reweighter.predict_weights(mc_df[dalitz_columns_dict.keys()], vote_function =  lambda x: np.mean(x, axis=0))
    GEN_weights = reweighter.predict_weights(gen_df[dalitz_columns_dict.keys()], vote_function =  lambda x: np.mean(x, axis=0))

    old_mc_count = mc_rdf.Count().GetValue()
    old_gen_count = gen_rdf.Count().GetValue()

    new_mc_count = MC_weights.sum()
    new_gen_count = GEN_weights.sum()

    oe = old_mc_count/old_gen_count
    ne = new_mc_count/new_gen_count

    print(f"OG Eff {oe} New Eff {ne} Change {(oe - ne)/oe}")

    dict_temp = {
                    'Data Spectrum' : id_to_scheme_dict[data_spec],
                    'MC Scheme' : id_to_scheme_dict[id_to_spec_dict[mc_spec]],
                    'Old Overall Eff' : f'{oe:.2E}',
                    'New Overall Eff' : f'{ne:.2E}',
                    'Relative Change (Systematic)' : f'{(oe - ne)/oe:.2E}'
                }

    df = pd.DataFrame(dict_temp, index=[0])
    txt_temp = f"txt_files/{mc_spec}_rw.txt"
    df.to_csv(txt_temp)

    mc_rdf_new = ROOT.RDF.MakeNumpyDataFrame({"new_weights": MC_weights})
    gen_rdf_new = ROOT.RDF.MakeNumpyDataFrame({"new_weights": GEN_weights})

    mc_rdf_new.Snapshot("DecayTree", f"weight_files/{mc_spec}_mc.root")
    gen_rdf_new.Snapshot("DecayTree", f"weight_files/{mc_spec}_gen.root")

    print(f"saved weights for {mc_spec}")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Reweighting")
    parser.add_argument('--Data_Spec', choices=["Z_m_p","Z_z_z","P_z_p","M_m_z", "P_z_pst", "norm7","norm8"], help = 'Data Spec')
    parser.add_argument('--MC_Spec', choices=["01_Z_m_p", "02_Z_m_p", "04_Z_m_p", "norm7_norm7","norm8_norm8"], help = 'MC Spec')
    parser.add_argument('--Rw', action = 'store_true')
    parser.add_argument('--Plot', action = 'store_true')
    parser.add_argument('--Opt', action = 'store_true')
    args = parser.parse_args()

    data_spec = args.Data_Spec
    mc_spec = args.MC_Spec
    rw = args.Rw
    opt = args.Opt
    plot = args.Plot

    mean = id_to_b_values_dict[mc_spec][0]
    window = id_to_b_values_dict[mc_spec][1]
    bmax = mean + window
    bmin = mean - window

    data_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/final_sample/{data_spec}.root"
    mc_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/final_sample/{id_to_spec_dict[mc_spec]}.root"
    gen_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_redecay/{mc_spec}_tree.root"


    peak = id_to_b_values_dict[mc_spec][2]
    sw_path = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/sw_files/sw_{data_spec}_{peak}.root"


    data_tchain = grab_files_and_chain(data_path, ["DecayTreeTuple_SIG"])
    mc_tchain = grab_files_and_chain(mc_path, ["DecayTreeTuple_SIG"])
    gen_tchain = grab_files_and_chain(gen_path, ["DecayTree"])

    sw_data_tchain = grab_files_and_chain(sw_path, ["SW_tree"])
    data_tchain.AddFriend(sw_data_tchain, "SW_tree")

    if plot:
        mc_w_tchain = grab_files_and_chain(f"weight_files/{mc_spec}_mc.root", ["DecayTree"])
        gen_w_tchain = grab_files_and_chain(f"weight_files/{mc_spec}_gen.root", ["DecayTree"])

        mc_tchain.AddFriend(mc_w_tchain,"DecayTree")
        gen_tchain.AddFriend(gen_w_tchain,"DecayTree")

    data_rdf = RDF(data_tchain)
    mc_rdf = RDF(mc_tchain)
    gen_rdf = RDF(gen_tchain)

    data_rdf = data_rdf.Filter(f"B_DTF_M < {bmax} && B_DTF_M > {bmin}")

    m = id_to_meson_string_dict[data_spec]

    sw_id = f"{data_spec}_{id_to_b_values_dict[mc_spec][2]}_yield_sw"

    if "norm" not in data_spec:
        dalitz_columns_dict = {
        "M2_d1d2" : f"{m[1]}{m[2]}",
        "M2_d2k" :  f"{m[2]}{m[3]}",
        "M2_d1k" :  f"{m[1]}{m[3]}",
        }
        gen_M2_d1d2 = id_to_genstring_dict[mc_spec][0]
        gen_M2_d1k = id_to_genstring_dict[mc_spec][1]
        gen_M2_d2k = id_to_genstring_dict[mc_spec][2]

        data_rdf = data_rdf.Define("M2_d1d2",  get_dalitz_string(["D1","D2"])) \
                            .Define("M2_d1k", get_dalitz_string(["D1","KST"])) \
                            .Define("M2_d2k", get_dalitz_string(["D2","KST"])) \

        mc_rdf = mc_rdf.Define("M2_d1d2",  get_dalitz_string(["D1","D2"])) \
                        .Define("M2_d1k", get_dalitz_string(["D1","KST"])) \
                        .Define("M2_d2k", get_dalitz_string(["D2","KST"])) \
                        .Define("og_weights", "1")

        gen_rdf =  gen_rdf.Define("M2_d1d2",  f"({gen_M2_d1d2})*1e6") \
                           .Define("M2_d1k",  f"({gen_M2_d1k})*1e6") \
                           .Define("M2_d2k",  f"({gen_M2_d2k})*1e6") \
                           .Define("og_weights", "1")

    if "norm" in data_spec:
        dalitz_columns_dict =  {
        "M2_d1d2" : f"{m[1]}{m[2]}",
        "M2_d2k" :  f"{m[2]}{m[3]}",
        "M2_d1k" :  f"{m[1]}{m[3]}",
        "M2_h1h2" : "K-#pi+",
        "M2_h1h3" : "K-#pi+",
        "M2_h1h4" : "K-#pi-",
        "M2_h2h3" : "#pi+#pi+",
        "M2_h2h4" : "#pi+#pi-",
        # "M2_h3h4" : "#pi+#pi-",
        # "M2_h1h2h3" : "K-#pi+#pi+",
        # "M2_h1h2h4" : "K-#pi+#pi-",
        # "M2_h1h3h4" : "K-#pi+#pi-",
        # "M2_h2h3h4" : "#pi+#pi+#pi-",
        }

        gen_M2_d1d2 = id_to_genstring_dict[mc_spec][0]
        gen_M2_d1k = id_to_genstring_dict[mc_spec][1]
        gen_M2_d2k = id_to_genstring_dict[mc_spec][2]

        gen_M2_h1h2 = id_to_genstring_dict[mc_spec][3]
        gen_M2_h1h3 = id_to_genstring_dict[mc_spec][4]
        gen_M2_h1h4 = id_to_genstring_dict[mc_spec][5]

        gen_M2_h2h3 = id_to_genstring_dict[mc_spec][6]
        gen_M2_h2h4 = id_to_genstring_dict[mc_spec][7]
        gen_M2_h3h4 = id_to_genstring_dict[mc_spec][8]

        gen_M2_h1h2h3 = id_to_genstring_dict[mc_spec][9]
        gen_M2_h1h2h4 = id_to_genstring_dict[mc_spec][10]
        gen_M2_h1h3h4 = id_to_genstring_dict[mc_spec][11]
        gen_M2_h2h3h4 = id_to_genstring_dict[mc_spec][12]

        data_rdf = data_rdf.Define("M2_d1d2", get_dalitz_string(["D1","D2"])) \
                            .Define("M2_d1k", get_dalitz_string(["D1","K"])) \
                            .Define("M2_d2k", get_dalitz_string(["D2","K"])) \
                            .Define("M2_h1h2",  get_dalitz_string(["D2H1","D2H2"])) \
                            .Define("M2_h1h3",  get_dalitz_string(["D2H1","D2H3"])) \
                            .Define("M2_h1h4",  get_dalitz_string(["D2H1","D2H4"])) \
                            .Define("M2_h2h3",  get_dalitz_string(["D2H2","D2H3"])) \
                            .Define("M2_h2h4",  get_dalitz_string(["D2H2","D2H4"])) \
                            .Define("M2_h3h4",  get_dalitz_string(["D2H3","D2H4"])) \
                            .Define("M2_h1h2h3", get_dalitz_string(["D2H1","D2H2","D2H3"])) \
                            .Define("M2_h1h2h4", get_dalitz_string(["D2H1","D2H2","D2H4"])) \
                            .Define("M2_h1h3h4", get_dalitz_string(["D2H1","D2H3","D2H4"])) \
                            .Define("M2_h2h3h4", get_dalitz_string(["D2H2","D2H3","D2H4"]))

        mc_rdf = mc_rdf.Define("M2_d1d2", get_dalitz_string(["D1","D2"])) \
                            .Define("M2_d1k", get_dalitz_string(["D1","K"])) \
                            .Define("M2_d2k", get_dalitz_string(["D2","K"])) \
                            .Define("M2_h1h2",   get_dalitz_string(["D2H1","D2H2"])) \
                            .Define("M2_h1h3",  get_dalitz_string(["D2H1","D2H3"])) \
                            .Define("M2_h1h4",  get_dalitz_string(["D2H1","D2H4"])) \
                            .Define("M2_h2h3",  get_dalitz_string(["D2H2","D2H3"])) \
                            .Define("M2_h2h4",  get_dalitz_string(["D2H2","D2H4"])) \
                            .Define("M2_h3h4",  get_dalitz_string(["D2H3","D2H4"])) \
                            .Define("M2_h1h2h3", get_dalitz_string(["D2H1","D2H2","D2H3"])) \
                            .Define("M2_h1h2h4", get_dalitz_string(["D2H1","D2H2","D2H4"])) \
                            .Define("M2_h1h3h4", get_dalitz_string(["D2H1","D2H3","D2H4"])) \
                            .Define("M2_h2h3h4", get_dalitz_string(["D2H2","D2H3","D2H4"])) \
                            .Define("og_weights", "1")

        gen_rdf = gen_rdf.Define("M2_d1d2",  f"({gen_M2_d1d2})*1e6") \
                            .Define("M2_d1k",  f"({gen_M2_d1k})*1e6") \
                            .Define("M2_d2k",  f"({gen_M2_d2k})*1e6") \
                            .Define("M2_h1h2", f"({gen_M2_h1h2})*1e6") \
                            .Define("M2_h1h3", f"({gen_M2_h1h3})*1e6") \
                            .Define("M2_h1h4", f"({gen_M2_h1h4})*1e6") \
                            .Define("M2_h2h3", f"({gen_M2_h2h3})*1e6") \
                            .Define("M2_h2h4", f"({gen_M2_h2h4})*1e6") \
                            .Define("M2_h3h4", f"({gen_M2_h3h4})*1e6") \
                            .Define("M2_h1h2h3", f"({gen_M2_h1h2h3})*1e6") \
                            .Define("M2_h1h2h4", f"({gen_M2_h1h2h4})*1e6") \
                            .Define("M2_h1h3h4", f"({gen_M2_h1h3h4})*1e6") \
                            .Define("M2_h2h3h4", f"({gen_M2_h2h3h4})*1e6") \
                            .Define("og_weights", "1")

    if opt:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            warnings.warn("deprecated", DeprecationWarning)
            optimise(data_rdf, mc_rdf, gen_rdf, dalitz_columns_dict, data_spec, mc_spec)
    if rw:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            warnings.warn("deprecated", DeprecationWarning)
            mc_rw(data_rdf, mc_rdf, gen_rdf, dalitz_columns_dict, data_spec, mc_spec)
    if plot:
        draw_1d_dalitz(data_rdf, mc_rdf, gen_rdf, data_spec, mc_spec, dalitz_columns_dict, "og", "pre_rw_1d_Dalitz")
        draw_1d_dalitz(data_rdf, mc_rdf, gen_rdf, data_spec, mc_spec, dalitz_columns_dict, "rw", "post_rw_1d_Dalitz")
