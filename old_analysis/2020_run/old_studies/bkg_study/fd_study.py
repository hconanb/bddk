from fit_essentials import *

# tla = ROOT.TLatex()
# tla.SetTextFont(32)
# tla.SetTextColor(1)
# tla.SetTextSize(0.03)
# tla.SetTextAlign(12)

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")

z_tree = z_data_file.Get("DecayTreeTuple")
zz_tree = zz_data_file.Get("DecayTreeTuple")
p_tree = p_data_file.Get("DecayTreeTuple")

peak = 0

b0_sig = f"abs(B_DTF_M - {B0mass}) < 50"
bp_sig = f"abs(B_DTF_M - {Bpmass}) < 50"

d1_z_sig = f"(abs(D1_M - {d0mass}) < {dwindow})"
d1_p_sig = f"(abs(D1_M - {dpmass}) < {dwindow})"

d2_z_sig = f"(abs(D2_M - {d0mass}) < {dwindow})"
d2_p_sig = f"(abs(D2_M - {dpmass}) < {dwindow})"

d1_z_sb = f"(abs(D1_M - {d0mass}) < {dmaxbkg} && abs(D1_M - {d0mass}) > {dminbkg})"
d1_p_sb = f"(abs(D1_M - {dpmass}) < {dmaxbkg} && abs(D1_M - {dpmass}) > {dminbkg})"

d2_z_sb = f"(abs(D2_M - {d0mass}) < {dmaxbkg} && abs(D2_M - {d0mass}) > {dminbkg})"
d2_p_sb = f"(abs(D2_M - {dpmass}) < {dmaxbkg} && abs(D2_M - {dpmass}) > {dminbkg})"

RDF = ROOT.ROOT.RDataFrame

tree_list = [zz_tree, z_tree, p_tree]
tag_list = ["zz","z","p"]

fd_rdf_array = []
plot_list = []

bbins = 100
bmin = 4800
bmax = 5600
dbins = 100
dmin = 1750
dmax = 2000

x2_cut_list = [500, 400, 300, 200, 150, 100, 90, 80, 70, 60, 50, 40, 30, 20, 10]
fd_cut_list = [0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20]

ROOT.gStyle.SetOptStat(11)

for tree,tag in zip(tree_list, tag_list):

    d1_pl = []
    d2_pl = []
    d1sb_pl = []
    d2sb_pl = []
    ch2_pl = []
    bm_pl = []
    bdtfm_pl = []
    bmsig_pl = []

    rdf_base = RDF(tree)

    for fd_cut in x2_cut_list:
        rdf_temp = rdf_base.Filter(f"(D1 {fd_cut} )")
        if tag == "z":
            rdf_ch1_d1sb = rdf_temp.Filter(d1_p_sb + "&&" + d2_p_sig)
            rdf_ch1_d2sb = rdf_temp.Filter(d1_p_sig + "&&" + d2_p_sb)
            rdf_ch2 = rdf_temp.Filter(d1_p_sb + "&&" + d2_p_sb)
            rdf_sig = rdf_temp.Filter(d1_p_sig + "&&" + d2_p_sig)
        if tag == "zz":
            rdf_ch1_d1sb = rdf_temp.Filter(d1_z_sb + "&&" + d2_z_sig)
            rdf_ch1_d2sb = rdf_temp.Filter(d1_z_sig + "&&" + d2_z_sb)
            rdf_ch2 = rdf_temp.Filter(d1_z_sb + "&&" + d2_z_sb)
            rdf_sig = rdf_temp.Filter(d1_z_sig + "&&" + d2_z_sig)
        if tag == "p":
            rdf_ch1_d1sb = rdf_temp.Filter(d1_z_sb + "&&" + d2_p_sig)
            rdf_ch1_d2sb = rdf_temp.Filter(d1_z_sig + "&&" + d2_p_sb)
            rdf_ch2 = rdf_temp.Filter(d1_z_sb + "&&" + d2_p_sb)
            rdf_sig = rdf_temp.Filter(d1_z_sig + "&&" + d2_p_sig)


        D1_M = rdf_temp.Histo1D((f"D1_M_x2<{fd_cut}", f"D1_M_{tag}_{fd_cut}", dbins, dmin, dmax), "D1_M")
        D2_M = rdf_temp.Histo1D((f"D2_M_x2<{fd_cut}", f"D2_M_{tag}_{fd_cut}", dbins, dmin, dmax), "D2_M")
        B_M = rdf_temp.Histo1D((f"B_M_x2<{fd_cut}", f"B_M_{tag}_{fd_cut}", bbins, bmin, bmax), "B_M")
        B_dtf_M = rdf_temp.Histo1D((f"B_dtf_M_x2<{fd_cut}", f"B_dtf_M_{tag}_{fd_cut}", bbins, bmin, bmax), "B_DTF_M")
        ch1_d1sb_hist = rdf_ch1_d1sb.Histo1D((f"B_d1sbM_x2<{fd_cut}", f"B_M_d1sb_{tag}_{fd_cut}", bbins, bmin, bmax), "B_M")
        ch1_d2sb_hist = rdf_ch1_d2sb.Histo1D((f"B_d2sbM_x2<{fd_cut}", f"B_M_d2sb_{tag}_{fd_cut}", bbins, bmin, bmax), "B_M")
        ch2_hist = rdf_ch2.Histo1D((f"B_ch2M_x2<{fd_cut}, B_dtf_chi2[0] < {fd_cut}", f"B_M_ch2_{tag}_{fd_cut}", bbins, bmin, bmax), "B_M")
        B_M_sig = rdf_sig.Histo1D((f"B_sig_x2<{fd_cut}, B_dtf_chi2[0] < {fd_cut}", f"B_M_ch2_{tag}_{fd_cut}", bbins, bmin, bmax), "B_M")
        # hist_list = [D1_M, D2_M, B_M, B_dtf_M, ch1_d1sb_hist, ch1_d2sb_hist, ch2_hist]
        # plot_list.extend(hist_list)
        d1_pl.append(D1_M)
        d2_pl.append(D2_M)
        ch2_pl.append(ch2_hist)
        bm_pl.append(B_M)
        bdtfm_pl.append(B_dtf_M)
        d1sb_pl.append(ch1_d1sb_hist)
        d2sb_pl.append(ch1_d2sb_hist)
        bmsig_pl.append(B_M_sig)

    plarray = [d1_pl, d2_pl, ch2_pl, bm_pl, bdtfm_pl, d1sb_pl, d2sb_pl, bmsig_pl]
    plnames = ["D1_M", "D2_M", "B_M_ch2", "B_M", "B_DTF_M", "B_M_D1sb", "B_M_D2sb", "B_M_sig"]

    for i,j in zip(plarray, plnames):
        canvas = ROOT.TCanvas("c1","c1")
        canvas.Divide(5,3 )
        for k in range(len(i)):
            print(k)
            canvas.cd(k+1)
            i[k].Draw()
            if not os.path.exists(f'{tag}/'):
                os.makedirs(f'{tag}/')
            canvas.SaveAs(f"{tag}/{j}.png")


# for i in plot_list:
#     canvas = ROOT.TCanvas("c1","c1")
#     hist = i.DrawCopy()
#     hist.SetStats(1)
#     canvas.Update()
#     name = hist.GetTitle()
#     saveplot(canvas, name)
