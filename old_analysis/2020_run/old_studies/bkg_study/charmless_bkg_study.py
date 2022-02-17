from bkg_essentials import *

ROOT.gStyle.SetOptStat(11)

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")

z_tree = z_data_file.Get("DecayTreeTuple")
zz_tree = zz_data_file.Get("DecayTreeTuple")
p_tree = p_data_file.Get("DecayTreeTuple")

peak = 0

RDF = ROOT.ROOT.RDataFrame

tree_list = [zz_tree, z_tree, p_tree]
tag_list = ["zz","z","p"]
dtype_list = [[d0mass, d0mass], [dpmass, dpmass], [d0mass, dpmass]]

fdx2_cut_list = ["0","0.5","1","1.5","2","2.5","3","3.5"]

bbins = 100
bmin = 4800
bmax = 5600

dwindow = 25.0
dmaxbkg = 85.0
dminbkg = 35.0

dbins = 100
dmin = 1750
dmax = 2000



fdx_l = []

fd_forward_cut = "D1_DIRA_ORIVX*D1_FDCHI2_ORIVX > 0 && D2_DIRA_ORIVX*D2_FDCHI2_ORIVX > 0"
fd_backward_cut = "D1_DIRA_ORIVX*D1_FDCHI2_ORIVX < 0 && D2_DIRA_ORIVX*D2_FDCHI2_ORIVX < 0"

for tree, tag, dtype in zip(tree_list, tag_list, dtype_list):

    rdf_base_0 = RDF(tree)
    rdf_base = rdf_base_0.Filter("D1_DIRA_ORIVX>0 && D2_DIRA_ORIVX>0")

    rdf_f_temp = rdf_base.Filter(f"D1_FDCHI2_ORIVX > 1 && D2_FDCHI2_ORIVX > 1")
    rdf_b_temp = rdf_base.Filter(f"D1_FDCHI2_ORIVX <= 1 && D2_FDCHI2_ORIVX <= 1")

    d1_sig = f"(abs(D1_M - {dtype[0]}) < {dwindow})"
    d1_usb =  f"((D1_M - {dtype[0]}) < {dmaxbkg} && (D1_M - {dtype[0]}) > {dminbkg})"
    d1_lsb =  f"(-(D1_M - {dtype[0]}) < {dmaxbkg} && -(D1_M - {dtype[0]}) > {dminbkg})"
    d1_asb = f"(abs(D1_M - {dtype[0]}) < {dmaxbkg} && abs(D1_M - {dtype[0]}) > {dminbkg})"

    d2_sig = f"(abs(D2_M - {dtype[1]}) < {dwindow})"
    d2_usb =  f"((D2_M - {dtype[1]}) < {dmaxbkg} && (D2_M - {dtype[1]}) > {dminbkg})"
    d2_lsb =  f"(-(D2_M - {dtype[1]}) < {dmaxbkg} && -(D2_M - {dtype[1]}) > {dminbkg})"
    d2_asb = f"(abs(D2_M - {dtype[1]}) < {dmaxbkg} && abs(D2_M - {dtype[1]}) > {dminbkg})"

    d1_rlist = [[d1_sig, "d1_sig"], [d1_asb, "d1_asb"]]
    d2_rlist = [[d2_sig, "d2_sig"], [d2_asb, "d2_asb"]]

    for d1_r in d1_rlist:
        for d2_r in d2_rlist:
            rdf_f = rdf_f_temp.Filter(d1_r[0] + "&&" + d2_r[0])
            rdf_b = rdf_b_temp.Filter(d1_r[0] + "&&" + d2_r[0])

            B_M_f = rdf_f.Histo1D((f"B_M_f_{d1_r[1]}_{d2_r[1]}_{tag}", f"B_M_f_{d1_r[1]}_{d2_r[1]}_{tag}", bbins, bmin, bmax), "B_M")
            B_M_b = rdf_b.Histo1D((f"B_M_b_{d1_r[1]}_{d2_r[1]}_{tag}", f"B_M_b_{d1_r[1]}_{d2_r[1]}_{tag}", bbins, bmin, bmax), "B_M")

            #B_dtf_M = rdf_next.Histo1D((f"B_dtf_M_{d1_r[1]}_{d2_r[1]}_{tag}", f"B_dtf_M_{d1_r[1]}_{d2_r[1]}_{tag}", bbins, bmin, bmax), "B_dtf_M")
            #D1_M = rdf.Histo1D((f"D1_M_{d1_r[1]}_{d2_r[1]}_{tag}", f"D1_M_{d1_r[1]}_{d2_r[1]}_{tag}", dbins, dmin, dmax), "D1_M")
            #D2_M = rdf.Histo1D((f"D2_M_{d1_r[1]}_{d2_r[1]}_{tag}", f"D2_M_{d1_r[1]}_{d2_r[1]}_{tag}", dbins, dmin, dmax), "D2_M")

            hist_list = [B_M_f, B_M_b]

            for hist in hist_list:
                canvas = ROOT.TCanvas("c1","c1")
                new_hist = hist.DrawCopy()
                name = new_hist.GetTitle()
                saveplot(canvas, name)
                
for tree, tag, dtype in zip(tree_list, tag_list, dtype_list):

    rdf_base_0 = RDF(tree)
    rdf_base = rdf_base_0.Filter("D1_DIRA_ORIVX>0 && D2_DIRA_ORIVX>0")

    bm_sig_sig_f = []
    bm_sig_asb_f = []
    bm_asb_sig_f = []
    bm_asb_asb_f = []

    bm_sig_asb_b = []
    bm_asb_sig_b = []
    bm_asb_asb_b = []
    bm_sig_sig_b = []

    for fdcut in fdx2_cut_list:
        rdf_f_temp = rdf_base.Filter(f"D1_FDCHI2_ORIVX > {fdcut} && D2_FDCHI2_ORIVX > {fdcut}")
        rdf_b_temp = rdf_base.Filter(f"D1_FDCHI2_ORIVX <= {fdcut} && D2_FDCHI2_ORIVX <= {fdcut}")

        d1_sig = f"(abs(D1_M - {dtype[0]}) < {dwindow})"
        d1_usb =  f"((D1_M - {dtype[0]}) < {dmaxbkg} && (D1_M - {dtype[0]}) > {dminbkg})"
        d1_lsb =  f"(-(D1_M - {dtype[0]}) < {dmaxbkg} && -(D1_M - {dtype[0]}) > {dminbkg})"
        d1_asb = f"(abs(D1_M - {dtype[0]}) < {dmaxbkg} && abs(D1_M - {dtype[0]}) > {dminbkg})"

        d2_sig = f"(abs(D2_M - {dtype[1]}) < {dwindow})"
        d2_usb =  f"((D2_M - {dtype[1]}) < {dmaxbkg} && (D2_M - {dtype[1]}) > {dminbkg})"
        d2_lsb =  f"(-(D2_M - {dtype[1]}) < {dmaxbkg} && -(D2_M - {dtype[1]}) > {dminbkg})"
        d2_asb = f"(abs(D2_M - {dtype[1]}) < {dmaxbkg} && abs(D2_M - {dtype[1]}) > {dminbkg})"

        d1_rlist = [[d1_sig, "d1_sig"], [d1_asb, "d1_asb"]]
        d2_rlist = [[d2_sig, "d2_sig"], [d2_asb, "d2_asb"]]

        for d1_r in d1_rlist:
            for d2_r in d2_rlist:

                rdf_f = rdf_f_temp.Filter(d1_r[0] + "&&" + d2_r[0])
                rdf_b = rdf_b_temp.Filter(d1_r[0] + "&&" + d2_r[0])

                B_M_f = rdf_f.Histo1D((f"B_M_f_{fdcut}_{d1_r[1]}_{d2_r[1]}_{tag}", f"B_M_f_{d1_r[1]}_{d2_r[1]}_{tag}", bbins, bmin, bmax), "B_M")
                B_M_b = rdf_b.Histo1D((f"B_M_b_{fdcut}_{d1_r[1]}_{d2_r[1]}_{tag}", f"B_M_b_{d1_r[1]}_{d2_r[1]}_{tag}", bbins, bmin, bmax), "B_M")

                if d1_r[1] == "d1_sig" and d2_r[1] == "d2_sig":
                    bm_sig_sig_f.append(B_M_f)
                    bm_sig_sig_b.append(B_M_b)
                if d1_r[1] == "d1_asb" and d2_r[1] == "d2_sig":
                    bm_asb_sig_f.append(B_M_f)
                    bm_asb_sig_b.append(B_M_b)
                if d1_r[1] == "d1_sig" and d2_r[1] == "d2_asb":
                    bm_sig_asb_f.append(B_M_f)
                    bm_sig_asb_b.append(B_M_b)
                if d1_r[1] == "d1_asb" and d2_r[1] == "d2_asb":
                    bm_asb_asb_f.append(B_M_f)
                    bm_asb_asb_b.append(B_M_b)



    savesplitplot(bm_sig_sig_f, 4, 2, "B_M_f_ss")
    savesplitplot(bm_sig_sig_b, 4, 2, "B_M_b_ss")

    savesplitplot(bm_asb_sig_f, 4, 2, "B_M_f_as")
    savesplitplot(bm_asb_sig_b, 4, 2, "B_M_b_as")

    savesplitplot(bm_sig_asb_f, 4, 2, "B_M_f_sa")
    savesplitplot(bm_sig_asb_b, 4, 2, "B_M_b_sa")

    savesplitplot(bm_asb_asb_f, 4, 2, "B_M_f_aa")
    savesplitplot(bm_asb_asb_b, 4, 2, "B_M_b_aa")
