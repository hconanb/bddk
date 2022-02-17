from fit_essentials import *

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")

z_tree = z_data_file.Get("DecayTreeTuple")
zz_tree = zz_data_file.Get("DecayTreeTuple")
p_tree = p_data_file.Get("DecayTreeTuple")

name = "data_ws"

dws = ROOT.RooWorkspace()

dws.factory("B_M[4800,5600]")
dws.factory("D1_M[1780,1950]")
dws.factory("D2_M[1780,1950]")


z_cut = f"(abs(D1_M - {dpmass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow})"
zz_cut = f"(abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow})"
p_cut = f"(abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {dpmass}) < {dwindow})"

# dws.factory(f"Bernstein:d1_bkg(D1_M, {{d10_s[1,0,10], d11_s[1,0,10], d12_s[1,0,10], d13_s[1,0,10]}})")
# dws.factory(f"Bernstein:d2_bkg(D2_M, {{d20_s[1,0,10], d21_s[1,0,10], d22_s[1,0,10], d23_s[1,0,10]}})")

dws.factory(f"Exponential:d1_bkg(D1_M, d1_s[0, -2, 2])")
dws.factory(f"Exponential:d2_bkg(D2_M, d2_s[0, -2, 2])")

dws.factory(f"Gaussian::d1_a(D1_M, mean_d1[1860,1800,1900], width_a_d1[8.0,0.01,50])")
dws.factory(f"Gaussian::d1_b(D1_M, mean_d1, width_b_d1[7.0,0.01,20])")
dws.factory(f"SUM::d1_fit(d1_a_frac[0.5,0,1]*d1_a, d1_b)")

dws.factory(f"Gaussian::d2_a(D2_M, mean_d2[1860,1800,1900], width_a_d2[8.0,0.01,50])")
dws.factory(f"Gaussian::d2_b(D2_M, mean_d2, width_b_d2[7.0,0.01,20])")
dws.factory(f"SUM::d2_fit(d2_a_frac[0.5,0,1]*d2_a, d2_b)")

dws.factory("PROD::d_sig_fit(d1_fit, d2_fit)")
dws.factory("PROD::d_bkg_fit(d1_bkg, d2_bkg)")

dws.factory(f"SUM::d_fit_all(n_sig[1000,0,1000000]*d_sig_fit, n_bkg[1000,0,1000000]*d_bkg_fit)")

pdf = dws.pdf("d_fit_all")
b_m = dws.var("B_M")
d1_mass = dws.var("D1_M")
d2_mass = dws.var("D2_M")

data_args = ROOT.RooArgSet(d1_mass, d2_mass, b_m)
ds = ROOT.RooDataSet("t1","t1", zz_tree, data_args)
pdf.fitTo(ds)

c1 = ROOT.TCanvas("c1","c1")
d1frame = d1_mass.frame()
ds.plotOn(d1frame)
pdf.plotOn(d1frame)
d1frame.Draw()
c1.SaveAs("d1.pdf")

c2 = ROOT.TCanvas("c2","c2")
d2frame = d2_mass.frame()
ds.plotOn(d2frame)
pdf.plotOn(d2frame)
d2frame.Draw()
c2.SaveAs("d2.pdf")

n_sig = dws.var("n_sig")
n_bkg = dws.var("n_bkg")

spData = ROOT.RooStats.SPlot("sData","An SPlot", ds, pdf, ROOT.RooArgList(n_sig, n_bkg))
sigwargs = spData.GetSWeightVars()
sigwargs.add(b_m)
sigwargs.add(d1_mass)
sigwargs.add(d2_mass)

c4 = ROOT.TCanvas("c4","c4")
d4frame = b_m.frame()
ds.plotOn(d4frame)
d4frame.Draw()
c4.SaveAs("bm.pdf")

ds_new = ROOT.RooDataSet("t1","t1", ds, sigwargs, "", "n_sig_sw")

c3 = ROOT.TCanvas("c3","c3")
d3frame = b_m.frame()
ds_new.plotOn(d3frame, ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
d3frame.Draw()
c3.SaveAs("bm_sw.pdf")

dws.Print()
