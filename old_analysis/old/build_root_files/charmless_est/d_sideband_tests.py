import sys
# from createXFD import *
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/")
from essentials import *

class spectrum_class:

  def __init__(self, name, d1_string, d2_string, d1_mass, d2_mass, d_strat):
    self.spec = name

    self.d1_string = d1_string
    self.d2_string = d2_string

    self.d1_mass = d1_mass
    self.d2_mass = d2_mass

    if self.spec == "st":
        self.d3_string = "(D^{*+} #rightarrow D^{0} #pi+)"
        self.d3_mass = 150

    self.d_strat = d_strat

z_c = spectrum_class("z", "D^{-}", "D^{+}", dpmass, dpmass, "e_dg_b")
zz_c = spectrum_class("zz", "#barD^{0}", "D^{0}", d0mass, d0mass, "e_dg_b")
p_c = spectrum_class("p", "#bar{D^{0}}", "D^{+}", d0mass, dpmass, "e_g")
m_c = spectrum_class("m", "D^{-}", "D^{0}", dpmass, d0mass, "e_g")
st_c = spectrum_class("st","#bar{D^{0}}", "D^{0}", d0mass, d0mass, "e_dg_b")
s_c = spectrum_class("s","D^{-}_{s}", "D^{+}" ,dsmass, dpmass, "e_g")
n7_c = spectrum_class("norm7","#bar{D^{0}}", "D^{0} #rightarrow k#pi#pi#pi", d0mass, d0mass, "e_dg_b")
n8_c = spectrum_class("norm8","D^{-}", "D^{0} #rightarrow k#pi#pi#pi", dpmass, d0mass, "e_g")


def get_dwindow_values(dwindow_ws, tag):

    if tag == "z" or tag == "zz" or tag == "st" or tag == "norm7":

        d1_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
        d2_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
        d1_std = dwindow_ws.var(f"width_a_{tag}_D1").getValV()
        d2_std = dwindow_ws.var(f"width_a_{tag}_D1").getValV()

    if tag == "p" or tag == "m" or tag == "s" or tag == "norm8":

        d1_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
        d2_mstart = dwindow_ws.var(f"mean_{tag}_D2").getValV()

        d1_std = dwindow_ws.var(f"width_a_{tag}_D1").getValV()
        d2_std = dwindow_ws.var(f"width_a_{tag}_D2").getValV()

    if tag == "st":
        d3_mstart = dwindow_ws.var(f"mean_{tag}_D3").getValV()
        d3_std = dwindow_ws.var(f"width_a_{tag}_D3").getValV()
        d3window = d3_std * 2.5

    return (d1_mstart, d2_mstart, d1_std, d2_std)

def build_d_window_ws(sc):


    spec = sc.spec
    # file_list = [f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root_nod/DATA/{spec}_spectrum_filtered.root", f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/DATA/z_none_ToT_spectrum_filtered.root"]

    file = ROOT.TFile(
        f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root_nod/DATA/{spec}_spectrum_filtered.root"
    )
    # file = ROOT.TFile(
    #     f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/DATA/z_none_ToT_spectrum_filtered.root"
    # )

    # for filename,nn in zip(file_list,["n","og"]):
    #
    # file = ROOT.TFile(filename)

    tree = file.Get("DecayTreeTuple")
    rdf_data_base = RDF(tree)

    dwindow_file = ROOT.TFile(f"d_window_root_files/d_{spec}_mass_fits.root", "READ")
    dwindow_ws = dwindow_file.Get(f"d_{spec}_mass_fits")
    d1_mstart, d2_mstart, d1_std, d2_std = get_dwindow_values(dwindow_ws, spec)
    dwindow_file.Close()

    # print (d1_mstart, d2_mstart, d1_std, d2_std)

    # rdf_sideband_ld1 = rdf_data_base.Filter(f"-({d1_mstart} - D1_M) > 2*{d1_std}")
    # rdf_sideband_ld2 = rdf_data_base.Filter(f"-({d2_mstart} - D1_M) > 2*{d2_std}")

    # rdf_sideband_hd1 = rdf_data_base.Filter(f"({d1_mstart} - D1_M) > 2*{d1_std}")
    # rdf_sideband_hd2 = rdf_data_base.Filter(f"({d2_mstart} - D1_M) > 2*{d2_std}")
    #
    # rdf_sideband_ld1ld2 = rdf_data_base.Filter(f"-({d1_mstart} - D1_M) > 2*{d1_std} && -({d2_mstart} - D1_M) > 2*{d2_std}")
    # rdf_sideband_ld1hd2 = rdf_data_base.Filter(f"-({d1_mstart} - D1_M) > 2*{d1_std} && ({d2_mstart} - D1_M) > 2*{d2_std}")
    #
    # rdf_sideband_hd1ld2 = rdf_data_base.Filter(f"({d1_mstart} - D1_M) > 2*{d1_std} && -({d2_mstart} - D1_M) > 2*{d2_std}")
    # rdf_sideband_hd1hd2 = rdf_data_base.Filter(f"({d1_mstart} - D1_M) > 2*{d1_std} && ({d2_mstart} - D1_M) > 2*{d2_std}")
    #
    # rdf_sideband_ad1 = rdf_data_base.Filter(f"abs({d1_mstart} - D1_M) > 2*{d1_std}")
    # rdf_sideband_ad2 = rdf_data_base.Filter(f"abs({d2_mstart} - D1_M) > 2*{d2_std}")
    #
    # rdf_sideband_ad1ld2 = rdf_data_base.Filter(f"abs({d1_mstart} - D1_M) > 2*{d1_std} && -({d2_mstart} - D1_M) > 2*{d2_std}")
    # rdf_sideband_ad1hd2 = rdf_data_base.Filter(f"abs({d1_mstart} - D1_M) > 2*{d1_std} && ({d2_mstart} - D1_M) > 2*{d2_std}")
    #
    # rdf_sideband_ld1ad2 = rdf_data_base.Filter(f"-({d1_mstart} - D1_M) > 2*{d1_std} && abs({d2_mstart} - D1_M) > 2*{d2_std}")
    # rdf_sideband_hd1ad2 = rdf_data_base.Filter(f"({d1_mstart} - D1_M) > 2*{d1_std} && abs({d2_mstart} - D1_M) > 2*{d2_std}")

    # print (f"{nn} pre", rdf_data_base.Count().GetValue())
    rdf_sideband_sig = rdf_data_base.Filter(f"(abs({d1_mstart} - D1_M) < 2*{d1_std}) && (abs({d2_mstart} - D2_M) < 2*{d2_std})")
    rdf_sideband_sb1_sg2 = rdf_data_base.Filter(f"(abs({d1_mstart} - D1_M) > 2*{d1_std}) && (abs({d1_mstart} - D1_M) < 5*{d1_std}) && (abs({d2_mstart} - D2_M) < 2*{d2_std})")
    rdf_sideband_sg1_sb2 = rdf_data_base.Filter(f"(abs({d1_mstart} - D1_M) < 2*{d1_std}) && (abs({d2_mstart} - D2_M) > 2*{d2_std}) && (abs({d2_mstart} - D2_M) < 5*{d2_std})")
    rdf_sideband_sb1_sb2 = rdf_data_base.Filter(f"(abs({d1_mstart} - D1_M) > 2*{d1_std}) && (abs({d1_mstart} - D1_M) < 5*{d1_std}) && (abs({d2_mstart} - D2_M) > 2*{d2_std}) && (abs({d2_mstart} - D2_M) < 5*{d2_std})")

    # rdf_sig_sidebands = rdf_data_base.Filter(f"abs(B_M - 5280) < 50 &&  (abs({d1_mstart} - D1_M) < 4*{d1_std}) &&  (abs({d2_mstart} - D1_M) < 4*{d2_std})")
    # rdf_sig_sidebands2 = rdf_data_base.Filter(f"abs(B_M - 5125) < 50 &&  (abs({d1_mstart} - D1_M) < 4*{d1_std}) &&  (abs({d2_mstart} - D1_M) < 4*{d2_std})")

    rdf_sig_sidebands_sb = rdf_data_base.Filter(f"((abs({d1_mstart} - D1_M) > 2*{d1_std}) && (abs({d1_mstart} - D1_M) < 4*{d1_std}) && (abs({d2_mstart} - D2_M) < 2*{d2_std})) || ((abs({d1_mstart} - D1_M) < 2*{d1_std}) && (abs({d2_mstart} - D2_M) > 2*{d2_std}) && (abs({d2_mstart} - D2_M) < 4*{d2_std}))")
    # rdf_sig_sidebands_sb2 = rdf_sig_sidebands2.Filter(f"(abs({d1_mstart} - D1_M) > 2*{d1_std}) && (abs({d1_mstart} - D1_M) < 4*{d1_std}) && (abs({d2_mstart} - D2_M) < 2*{d2_std}) || (abs({d1_mstart} - D1_M) < 2*{d1_std}) && (abs({d2_mstart} - D2_M) > 2*{d2_std}) && (abs({d2_mstart} - D2_M) < 4*{d2_std})")


    # Histo1D((f"B_M_f_{d1_r[1]}_{d2_r[1]}_{tag}", f"B_M_f_{d1_r[1]}_{d2_r[1]}_{tag}", bbins, bmin, bmax), "B_M")

    # hist_sig = rdf_sideband_sig.Histo1D(("sig", "sig", 100, 4800, 5600), 'B_M')
    # hist_sb1sg2 = rdf_sideband_sb1_sg2.Histo1D(("sb1sg2", "sb1sg2", 100, 4800, 5600), 'B_M')
    # hist_sg1sb2 = rdf_sideband_sg1_sb2.Histo1D(("sg1sb2", "sg1sb2", 100, 4800, 5600), 'B_M')
    # hist_sb1sb2 = rdf_sideband_sb1_sb2.Histo1D(("sb1sb2", "sb1sb2", 100, 4800, 5600), 'B_M')

    ROOT.gStyle.SetOptStat("ne")
    ROOT.gStyle.SetPalette(ROOT.kBird)
   # hcontz->Draw("CONTZ");

    c1 = ROOT.TCanvas("c1","c1")
    hist_sss1d = rdf_sig_sidebands_sb.Histo1D(("D-D+K*0: Cross - Signal", "D-D+K*0: Cross - Signal", 100, 4800, 5600), 'B_M')
    # hist_sss = rdf_sig_sidebands.Histo2D(("B->D-D+K*0, 5280 #pm 50", f"B->D-D+K*0, 5280 #pm 50", 80, 1800, 1960, 80, 1800, 1960), "D1_M", "D2_M")
    hist_sss1d.SetTitle("B->D-D+K*0, 5280 #pm 50")
    hist_sss1d.GetXaxis().SetTitle("D-D+K*0 MeV")
    # hist_sss.GetYaxis().SetTitle("D+")
    hist_sss1d.Draw("")
    c1.SaveAs(f"test_sss1d.png")

    # c2 = ROOT.TCanvas("c2","c2")
    # # hist_sss2 = rdf_sig_sidebands2.Histo2D(("B->D-D+K*0, 5125 #pm 50", f"B->D-D+K*0, 5125 #pm 50", 80, 1800, 1960, 80, 1800, 1960), "D1_M", "D2_M")
    # # hist_sss2.SetTitle("B->D-D+K*0, 5125 #pm 50")
    # # hist_sss2.GetXaxis().SetTitle("D-")
    # # hist_sss2.GetYaxis().SetTitle("D+")
    # # hist_sss2.Draw("CONTZ")
    # c2.SaveAs(f"test_sss21d.png")

    # for rdf_t, nn in zip([rdf_sideband_sig, rdf_sideband_sb1_sg2, rdf_sideband_sg1_sb2, rdf_sideband_sb1_sb2],["sig", "sb1sg2", "sg1sb2", "sb1sb2"]):
    #
    #     c1 = ROOT.TCanvas("c1","c1")
    #     c1.Divide(1,2)
    #     hist_d1 = rdf_t.Histo1D((f"d1_{nn}", f"d1_{nn}", 100, 1800, 1960), 'D1_M')
    #     hist_d2 = rdf_t.Histo1D((f"d2_{nn}", f"d2_{nn}", 100, 1800, 1960), 'D2_M')
    #
    #     c1.cd(1)
    #     hist_d1.Draw()
    #
    #     c1.cd(2)
    #     hist_d2.Draw()
    #
    #     c1.SaveAs(f"test_d1d2_{nn}.png")
    #
    # for rdf_t, nn in zip([rdf_sideband_sig, rdf_sideband_sb1_sg2, rdf_sideband_sg1_sb2, rdf_sideband_sb1_sb2], ["sig", "sb1sg2", "sg1sb2", "sb1sb2"]):
    #
    #     c1 = ROOT.TCanvas("c1","c1")
    #     c1.Divide(1,2)
    #
    #     hist_bm = rdf_t.Histo1D((f"{nn}", f"{nn}", 100, 4800, 5600), 'B_M')
    #     hist_bdtfm = rdf_t.Histo1D((f"{nn}", f"{nn}", 100, 4800, 5600), 'B_DTF_M')
    #
    #     c1.cd(1)
    #     hist_bm.Draw()
    #
    #     c1.cd(2)
    #     hist_bdtfm.Draw()
    #
    #     c1.SaveAs(f"test_bm_bdtfm_{nn}.png")

build_d_window_ws(z_c)
