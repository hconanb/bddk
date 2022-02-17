import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
from uncertainties import correlated_values
from uncertainties.umath import sqrt
import root_pandas as rp
import itertools
RDF = ROOT.ROOT.RDataFrame


branch_dict = {
"D1H1D1H2D2H1D2H2KSTH1" : "m(K#piK#piK)",
"D1H1D1H2D2H1D2H2KSTH1KSTH2" : "m(K#piK#piK#pi)",
"D2H1D2H2KSTH2": "m(D0(K#pi)K^{*0}(pi))",
"D1H1D1H2KSTH2": "m(#bar{D0}(K#pi)K^{*0}(pi))",

"D1H1" : "Kp_D0bar",
"D1H2" : "pim_D0bar",
"D2H1" : "Km_D",
"D2H2" : "pip1_Dp",
"D2H3" : "pip2_Dp",
"KSTH1" : "Kp_Kst0",
"KSTH2" : "pim_Kst0",
}

def rdf_plot(rdf, column, nbins, min, max, title):
# # D1_M = rdf.Histo1D((f"D1_M", f"D1_M", dbins, dmin, dmax), "D1_M")
    hist = rdf.Histo1D((f"{column}", f"{column}", nbins, min, max), f"{column}")
    bincan = ((max-min)/nbins)

    #hist = rdf.Histo1D((f"{column}", f"{column}", f"{nbins}", f"{min}", f"{max}"), f"{column}")
    c1 = ROOT.TCanvas("c1","c1")
    hist_plot = hist.DrawCopy()
    hist_plot.GetXaxis().SetTitle(branch_dict[column])
    hist_plot.GetYaxis().SetTitle(f"Candidates / ({bincan} MeV)")

    linemax = hist_plot.GetMaximum() + 10
    print(linemax)
    if column == "D1H1D1H2D2H1D2H2KSTH1":
        cutline = ROOT.TLine(4850, 0, 4850, linemax)
        cutline.SetLineColor(ROOT.kRed)
        cutline.SetLineWidth(5)
        cutline.SetLineStyle(ROOT.kDashed)
        cutline.Draw()
        cutline2 = ROOT.TLine(5025, 0, 5025, linemax)
        cutline2.SetLineColor(ROOT.kGreen)
        cutline2.SetLineWidth(5)
        cutline2.SetLineStyle(ROOT.kDashed)
        cutline2.Draw()
        cutline3 = ROOT.TLine(5175, 0, 5175, linemax)
        cutline3.SetLineColor(ROOT.kBlue)
        cutline3.SetLineWidth(5)
        cutline3.SetLineStyle(ROOT.kDashed)
        cutline3.Draw()

    c1.SaveAs(f"{title}.png")

zz_file = ROOT.TFile(f"zz_track_combinations.root")
dtt = zz_file.Get(f"DTT")
rdf_base = RDF(dtt)

# rdf_plot(rdf_base, "D1H1D1H2D2H1D2H2KSTH1", 100, 4200, 5600, "5track_zz")

cut1 = "(abs(D1H1D1H2D2H1D2H2KSTH1 - 5280) < 50)"
cut2 = "(abs(D1H1D1H2D2H1D2H2KSTH1 - 5125) < 50)"
cut3 = "(abs(D1H1D1H2D2H1D2H2KSTH1 - 4950) < 50)"
cutkpp = "(D1H1D1H2KSTH2 < 2050)"

rdf_filter = rdf_base.Filter(f"!{cut1} && !{cut2} && !{cut3} && !{cutkpp}")
                     # .Filter(f"!{cut2}")\
                     # .Filter(f"!{cut3}")\

rdf_3_filter = rdf_base.Filter(f"!{cutkpp}")
                     # .Filter(f"!{cut2}")\
                     # .Filter(f"!{cut3}")\

# rdf_plot(rdf_base, "D1H1D1H2D2H1D2H2KSTH1KSTH2", 100, 4800, 5600, "6track_zz_prefilter")
# rdf_plot(rdf_filter, "D1H1D1H2D2H1D2H2KSTH1KSTH2", 100, 4800, 5600, "6track_zz_postfilter")
rdf_plot(rdf_base, "D2H1D2H2KSTH2", 100, 1900, 3200, "4_track_peak_1")
rdf_plot(rdf_base, "D1H1D1H2KSTH2", 100, 1900, 3200, "4_track_peak_2")
# rdf_plot(rdf_filter, "D1H1D1H2D2H1D2H2KSTH1KSTH2", 100, 4800, 5600, "kpipiandall_cut_8track")
