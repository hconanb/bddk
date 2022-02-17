import sys
from createXFD import *
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/")
from essentials import *

class spectrum_class:

  def __init__(self, name, d1_string, d2_string, d1_mass, d2_mass, d_strat):
    self.spec = name

    self.d1_string = d1_string
    self.d2_string = d2_string

    self.d1_mass = d1_mass
    self.d2_mass = d2_mass

    if self.spec == "Z_m_pst" or self.spec == "P_z_pst":
        self.d3_string = "(D^{*+} #rightarrow D^{0} #pi+)"
        self.d3_mass = 150
    if self.spec == "Z_mst_p":
        self.d3_string = "(D^{*-} #rightarrow #bar{D^{0}} #pi-)"
        self.d3_mass = 150

    self.d_strat = d_strat

z_c = spectrum_class("Z_m_p", "D^{-}", "D^{+}", dpmass, dpmass, "e_dg_b")
zmst_c = spectrum_class("Z_mst_p", "D^{-}", "D^{+}", dpmass, dpmass, "e_g")
zpst_c = spectrum_class("Z_m_pst", "D^{-}", "D^{+}", dpmass, dpmass, "e_g")
zz_c = spectrum_class("Z_z_z", "#barD^{0}", "D^{0}", d0mass, d0mass, "e_dg_b")
p_c = spectrum_class("P_z_p", "#bar{D^{0}}", "D^{+}", d0mass, dpmass, "e_g")
m_c = spectrum_class("M_m_z", "D^{-}", "D^{0}", dpmass, d0mass, "e_g")
st_c = spectrum_class("P_z_pst","#bar{D^{0}}", "D^{0}", d0mass, d0mass, "e_dg_b")
s_c = spectrum_class("Zs_sm_p","D^{-}_{s}", "D^{+}" ,dsmass, dpmass, "e_g")
n7_c = spectrum_class("norm7","#bar{D^{0}}", "D^{0} #rightarrow k#pi#pi#pi", d0mass, d0mass, "e_dg_b")
n8_c = spectrum_class("norm8","D^{-}", "D^{0} #rightarrow k#pi#pi#pi", dpmass, d0mass, "e_g")


def get_dwindow_values(dwindow_ws, spec):

    if spec == "Z_m_p" or spec == "Z_z_z" or spec == "P_z_pst" or spec == "norm7":

        d1_mstart = dwindow_ws.var(f"mean_{spec}_D1").getValV()
        d2_mstart = dwindow_ws.var(f"mean_{spec}_D1").getValV()
        d1_std = dwindow_ws.var(f"width_a_{spec}_D1").getValV()
        d2_std = dwindow_ws.var(f"width_a_{spec}_D1").getValV()

    if spec == "P_z_p" or spec == "M_m_z" or spec == "Zs_sm_p" or spec == "norm8" or spec == "Z_mst_p" or spec == "Z_m_pst":

        d1_mstart = dwindow_ws.var(f"mean_{spec}_D1").getValV()
        d2_mstart = dwindow_ws.var(f"mean_{spec}_D2").getValV()

        d1_std = dwindow_ws.var(f"width_a_{spec}_D1").getValV()
        d2_std = dwindow_ws.var(f"width_a_{spec}_D2").getValV()

    if spec == "P_z_pst" or spec == "Z_m_pst" or spec == "Z_mst_p":
        d3_mstart = dwindow_ws.var(f"mean_{spec}_D3").getValV()
        d3_std = dwindow_ws.var(f"width_a_{spec}_D3").getValV()
        d3window = d3_std * 2.5

    d1window = d1_std * 2
    d2window = d2_std * 2

    if spec == "P_z_pst" or spec == "Z_m_pst":
        mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window} && abs(D2stmD_M - {d3_mstart}) < {d2window}) "
    elif spec == "Z_mst_p":
        mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window} && abs(D1stmD_M - {d3_mstart}) < {d2window}) "
    else:
        mass_cut = f"(abs(D1_M - {d1_mstart}) < {d1window} && abs(D2_M - {d2_mstart}) < {d2window})"

    return mass_cut

def plot_d_window_fits(sc):

    file = ROOT.TFile(f"d_window_root_files/d_{sc.spec}_mass_fits.root")
    dws = file.Get(f"d_{sc.spec}_mass_fits")
    var_d1 = dws.var("D1_M")
    var_d2 = dws.var("D2_M")



    all_data_sets = dws.data("data_sets")
    all_fit = dws.obj(f"{sc.spec}_d_fit_Pdf")
    all_cats = dws.obj(f"{sc.spec}_d_cats")

    vlist = [var_d1, var_d2]
    slist = [f"{sc.spec}_D1", f"{sc.spec}_D2"]

    if "mst" in sc.spec:
        var_d3 = dws.var("D1stmD_M")
        vlist.append(var_d3)
        slist.append(f"{sc.spec}_D3")

    if "pst" in sc.spec:
        var_d3 = dws.var("D2stmD_M")
        vlist.append(var_d3)
        slist.append(f"{sc.spec}_D3")


    for var, vstring in zip(vlist,slist):
        print(var)
        frame = var.frame()

        all_data_sets.plotOn(frame, ROOT.RooFit.Cut(f"{sc.spec}_d_cats=={sc.spec}_d_cats::{vstring}_spectrum"), ROOT.RooFit.Name(f"{vstring}_data"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, f"{vstring}_spectrum"), ROOT.RooFit.Components(f"{vstring}_fit"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name(f"{vstring}_fit"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, f"{vstring}_spectrum"), ROOT.RooFit.Components(f"{vstring}_bkg"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name(f"{vstring}_bkg"))
        all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, f"{vstring}_spectrum"), ROOT.RooFit.Components(f"{vstring}_signal"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet), ROOT.RooFit.Name(f"{vstring}_signal"))
        if "dg" in sc.d_strat and sc.spec != "st":
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, f"{vstring}_spectrum"), ROOT.RooFit.Components(f"{vstring}_a"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"{vstring}_ga"))
            all_fit.plotOn(frame, ROOT.RooFit.Slice(all_cats, f"{vstring}_spectrum"), ROOT.RooFit.Components(f"{vstring}_b"), ROOT.RooFit.ProjWData(ROOT.RooArgSet(all_cats), all_data_sets), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"{vstring}_gb"))
        d_chi2 = frame.chiSquare(f"{vstring}_fit", f"{vstring}_data")
        if sc.d_strat == "e_dg_b":
            d_mean = dws.obj(f"mean_{sc.spec}_D1")
            d_stda = dws.var(f"width_a_{sc.spec}_D1")
            d_stdb = dws.var(f"width_b_{sc.spec}_D1")
            d_sigfrac = dws.var(f"{vstring}_a_frac")
        if sc.d_strat == "e_dg_s":
            d_mean = dws.obj(f"mean_{vstring}")
            d_stda = dws.var(f"width_a_{sc.spec}_D1")
            d_stdb = dws.var(f"width_b_{vstring}")
            d_sigfrac = dws.var(f"{vstring}_a_frac")
        if "_g" in sc.d_strat or "D3" in vstring:
            d_mean = dws.obj(f"mean_{vstring}")
            d_stda = dws.var(f"width_a_{vstring}")
            print(vstring)

        d_nyield = dws.var(f"n_{vstring}_signal")
        d_nbkg = dws.var(f"n_{vstring}_bkg")

        dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
        dtpave.SetFillStyle(0)
        dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
        dtpave.AddText(f"mean: {round(d_mean.getValV(),3)}")
        if "dg" in sc.d_strat and "D3" not in vstring:
            dtpave.AddText(f"nsig_a: {round(d_nyield.getValV()*d_sigfrac.getValV(),3)}")
            dtpave.AddText(f"nsig_b: {round(d_nyield.getValV()*(1-d_sigfrac.getValV()),3)}")
            dtpave.AddText(f"width_a: {round(d_stda.getValV(),3)}")
            dtpave.AddText(f"width_b: {round(d_stdb.getValV(),3)}")
        if "_g" in sc.d_strat or "D3" in vstring:
            dtpave.AddText(f"nsig_a: {round(d_nyield.getValV(),3)}")
            dtpave.AddText(f"width_a: {round(d_stda.getValV(),3)}")
        dtpave.AddText(f"nbkg: {round(d_nbkg.getValV(),3)}")
        frame.addObject(dtpave)

        p = residualPlot()
        p.pt.cd()
        xaxis = frame.GetXaxis()
        xaxis.SetTickLength(0)
        xaxis.SetNdivisions(0)
        xaxis.SetLabelSize(0)

        leg = ROOT.TLegend(0.75,0.5,0.95,0.87)
        leg.AddEntry(f"{vstring}_data","Data", "P");
        leg.AddEntry(frame.findObject(f"{vstring}_fit"),"Signal + background","L");
        leg.AddEntry(frame.findObject(f"{vstring}_bkg"),"Background only", "L");
        leg.AddEntry(frame.findObject(f"{vstring}_signal"),"Signal only", "L");
        if "dg" in sc.d_strat:
            leg.AddEntry(frame.findObject(f"{vstring}_ga"),"Signal A", "L");
            leg.AddEntry(frame.findObject(f"{vstring}_gb"),"Signal B", "L");
        frame.Draw()
        leg.Draw()
        p.pb.cd()

        hpull = frame.pullHist(f"{vstring}_data", f"{vstring}_fit")
        pull = var.frame()
        pull.addPlotable(hpull, "P")
        pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
        pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize() * p.blabelratio)
        pull.GetYaxis().SetLabelSize(
            ROOT.gStyle.GetLabelSize("Y") * (1 + p.padratio) / (2 * p.padratio)
        )
        pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
        pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize() * p.blabelratio)
        pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset() / p.blabelratio)
        pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength() * p.blabelratio)
        pull.GetYaxis().SetTitle("Pull")
        pull.GetYaxis().SetNdivisions(202, False)
        pull.GetYaxis().SetRangeUser(-3, 3)
        pull.GetYaxis().CenterTitle()
        pull.Draw("AP")

        if vstring == f"{sc.spec}_D1":
            fancy_dstring = sc.d1_string
        if vstring == f"{sc.spec}_D2":
            fancy_dstring = sc.d2_string

        pull.GetXaxis().SetTitle(f"{fancy_dstring} Mass [MeV]")
        save_png(p, "d_window_plots", f"{vstring}", rpflag=1)

def build_d_window_ws(sc, tree_chain):

    dname = f"d_{sc.spec}_mass_fits"
    dws = ROOT.RooWorkspace(dname)
    dfithalf = 50

    dws.factory(f"D1_M[{sc.d1_mass - dfithalf},{sc.d1_mass + dfithalf}]")
    dws.factory(f"D2_M[{sc.d2_mass - dfithalf},{sc.d2_mass + dfithalf}]")
    dws.factory(f"D1stmD_M[140,150]")
    dws.factory(f"D2stmD_M[140,150]")

    dws.factory(f"Exponential:{sc.spec}_D1_bkg(D1_M, c0_{sc.spec}_D1[0, -2, 3])")
    dws.factory(f"Exponential:{sc.spec}_D2_bkg(D2_M, c0_{sc.spec}_D2[0, -2, 3])")

    if sc.spec != "norm7" and sc.spec != "norm8":
        maxn = 200000
    if sc.spec == "norm7" or sc.spec == "norm8":
        maxn = 20000000

    if "P_z_pst" in sc.spec or "Z_m_pst" in sc.spec:
        dws.factory(f"Exponential:{sc.spec}_D3_bkg(D2stmD_M, c0_{sc.spec}_D3[0, -2, 3])")
        dws.factory(f"Gaussian::{sc.spec}_D3_signal(D2stmD_M, mean_{sc.spec}_D3[150,140,160], width_a_{sc.spec}_D3[1,0,17])")
        dws.factory(f"SUM::{sc.spec}_D3_fit(n_{sc.spec}_D3_signal[1000,0,{maxn}]*{sc.spec}_D3_signal, n_{sc.spec}_D3_bkg[1000,0,{maxn}]*{sc.spec}_D3_bkg)")

    if "Z_mst_p" in sc.spec:
        dws.factory(f"Exponential:{sc.spec}_D3_bkg(D1stmD_M, c0_{sc.spec}_D3[0, -2, 3])")
        # dws.factory(f"Gaussian::{sc.spec}_D3_signal(D1stmD_M, mean_{sc.spec}_D3[150,140,160], width_a_{sc.spec}_D3[1,0,17])")
        dws.factory(f"BifurGauss::{sc.spec}_D3_signal(D1stmD_M, mean_{sc.spec}_D3[150,140,160],  width_a_{sc.spec}_D3[1,0,17],  width_b_{sc.spec}_D3[4,0,17])")
        # dws.factory(f"RooBifurGaussExp::{sc.spec}_D3_signal(D1stmD_M, mean_{sc.spec}_D3[150,140,160],width_a_{sc.spec}_D3[1,0,17]),width_a_{sc.spec}_D3[1,0,17]),alpha_1_{spec}_{shape_flag}[2,0.00,10.0],alpha_2_{spec}_{shape_flag}[8,0.00,10.0])")
        dws.factory(f"SUM::{sc.spec}_D3_fit(n_{sc.spec}_D3_signal[1000,0,{maxn}]*{sc.spec}_D3_signal, n_{sc.spec}_D3_bkg[1000,0,{maxn}]*{sc.spec}_D3_bkg)")


    if sc.d_strat == "e_dg_b":
        dws.factory(f"Gaussian::{sc.spec}_D1_a(D1_M, mean_{sc.spec}_D1[{sc.d1_mass},{sc.d1_mass-10},{sc.d1_mass+10}], width_a_{sc.spec}_D1[13,10,17])")
        dws.factory(f"Gaussian::{sc.spec}_D1_b(D1_M, mean_{sc.spec}_D1, width_b_{sc.spec}_D1[6.0,1,9])")
        dws.factory(f"SUM::{sc.spec}_D1_signal({sc.spec}_D1_a_frac[0.5,0,1]*{sc.spec}_D1_a, {sc.spec}_D1_b)")
        dws.factory(f"Gaussian::{sc.spec}_D2_a(D2_M, mean_{sc.spec}_D1, width_a_{sc.spec}_D1)")
        dws.factory(f"Gaussian::{sc.spec}_D2_b(D2_M, mean_{sc.spec}_D1, width_b_{sc.spec}_D1)")
        dws.factory(f"SUM::{sc.spec}_D2_signal({sc.spec}_D2_a_frac[0.5,0,1]*{sc.spec}_D2_a, {sc.spec}_D2_b)")
        dws.factory(f"SUM::{sc.spec}_D1_fit(n_{sc.spec}_D1_signal[1000,0,{maxn}]*{sc.spec}_D1_signal, n_{sc.spec}_D1_bkg[1000,0,{maxn}]*{sc.spec}_D1_bkg)")
        dws.factory(f"SUM::{sc.spec}_D2_fit(n_{sc.spec}_D2_signal[1000,0,{maxn}]*{sc.spec}_D2_signal, n_{sc.spec}_D2_bkg[1000,0,{maxn}]*{sc.spec}_D2_bkg)")

    if sc.d_strat == "e_dg_s":
        dws.factory(f"Gaussian::{sc.spec}_D1_a(D1_M, mean_{sc.spec}_D1[{sc.d1_mass},{sc.d1_mass-10},{sc.d1_mass+10}], width_a_{sc.spec}_D1[15,1,17])")
        dws.factory(f"Gaussian::{sc.spec}_D1_b(D1_M, mean_{sc.spec}_D1, width_b_{sc.spec}_D1[3.0,1,9])")
        dws.factory(f"SUM::{sc.spec}_D1_signal({sc.spec}_D1_a_frac[0.5,0,1]*{sc.spec}_D1_a, {sc.spec}_D1_b)")
        dws.factory(f"Gaussian::{sc.spec}_D2_a(D2_M, mean_{sc.spec}_D2[{sc.d2_mass},{sc.d2_mass-10},{sc.d2_mass+10}], width_a_{sc.spec}_D1)")
        dws.factory(f"Gaussian::{sc.spec}_D2_b(D2_M, mean_{sc.spec}_D2, width_b_{sc.spec}_D2[3.0,1,11])")
        dws.factory(f"SUM::{sc.spec}_D2_signal({sc.spec}_D2_a_frac[0.5,0,1]*{sc.spec}_D2_a, {sc.spec}_D2_b)")
        dws.factory(f"SUM::{sc.spec}_D1_fit(n_{sc.spec}_D1_signal[1000,0,{maxn}]*{sc.spec}_D1_signal, n_{sc.spec}_D1_bkg[1000,0,{maxn}]*{sc.spec}_D1_bkg)")
        dws.factory(f"SUM::{sc.spec}_D2_fit(n_{sc.spec}_D2_signal[1000,0,{maxn}]*{sc.spec}_D2_signal, n_{sc.spec}_D2_bkg[1000,0,{maxn}]*{sc.spec}_D2_bkg)")

    if sc.d_strat == "e_g":
        dws.factory(f"Gaussian::{sc.spec}_D1_signal(D1_M, mean_{sc.spec}_D1[{sc.d1_mass},{sc.d1_mass-10},{sc.d1_mass+10}], width_a_{sc.spec}_D1[13,5,17])")
        dws.factory(f"Gaussian::{sc.spec}_D2_signal(D2_M, mean_{sc.spec}_D2[{sc.d2_mass},{sc.d2_mass-10},{sc.d2_mass+10}], width_a_{sc.spec}_D2[13,5,17])")
        dws.factory(f"SUM::{sc.spec}_D1_fit(n_{sc.spec}_D1_signal[1000,0,{maxn}]*{sc.spec}_D1_signal, n_{sc.spec}_D1_bkg[1000,0,{maxn}]*{sc.spec}_D1_bkg)")
        dws.factory(f"SUM::{sc.spec}_D2_fit(n_{sc.spec}_D2_signal[1000,0,{maxn}]*{sc.spec}_D2_signal, n_{sc.spec}_D2_bkg[1000,0,{maxn}]*{sc.spec}_D2_bkg)")

    all_cats = ROOT.RooCategory(f"{sc.spec}_d_cats", f"{sc.spec}_d_cats")
    all_cats.defineType(f"{sc.spec}_D1_spectrum")
    all_cats.defineType(f"{sc.spec}_D2_spectrum")

    d1_mass_var = dws.var("D1_M")
    d2_mass_var = dws.var("D2_M")

    if "Z_mst_p" in sc.spec:
        d3_mass_var = dws.var("D1stmD_M")
    if "Z_m_pst" in sc.spec:
        d3_mass_var = dws.var("D2stmD_M")

    if "st" not in sc.spec:
        data_args = ROOT.RooArgSet(d1_mass_var, d2_mass_var)
    else:
        data_args = ROOT.RooArgSet(d1_mass_var, d2_mass_var, d3_mass_var)


    data_set = ROOT.RooDataSet(f"{sc.spec}_data", f"{sc.spec}_data", tree_chain, data_args)

    if "pst" not in sc.spec and "mst" not in sc.spec:
        all_data_sets = ROOT.RooDataSet(
            "data_sets",
            "data_sets",
            data_args,
            ROOT.RooFit.Index(all_cats),
            ROOT.RooFit.Import(f"{sc.spec}_D1_spectrum", data_set),
            ROOT.RooFit.Import(f"{sc.spec}_D2_spectrum", data_set),
        )
    if "st" in sc.spec:
        all_data_sets = ROOT.RooDataSet(
            "data_sets",
            "data_sets",
            data_args,
            ROOT.RooFit.Index(all_cats),
            ROOT.RooFit.Import(f"{sc.spec}_D1_spectrum", data_set),
            ROOT.RooFit.Import(f"{sc.spec}_D2_spectrum", data_set),
            ROOT.RooFit.Import(f"{sc.spec}_D3_spectrum", data_set)
        )

    all_fit = ROOT.RooSimultaneous(f"{sc.spec}_d_fit_Pdf", f"{sc.spec}_d_fit_Pdf", all_cats)

    d1_model = dws.pdf(f"{sc.spec}_D1_fit")
    d2_model = dws.pdf(f"{sc.spec}_D2_fit")

    all_fit.addPdf(d1_model, f"{sc.spec}_D1_spectrum")
    all_fit.addPdf(d2_model, f"{sc.spec}_D2_spectrum")

    if "pst" in sc.spec or "mst" in sc.spec:
        d3_model = dws.pdf(f"{sc.spec}_D3_fit")
        all_fit.addPdf(d3_model, f"{sc.spec}_D3_spectrum")

    dws.Print()
    all_fit.fitTo(all_data_sets, ROOT.RooFit.PrintLevel(0))
    dws.Import(all_data_sets)
    dws.Import(all_fit)
    dws.writeToFile(f"d_window_root_files/{dname}.root")

def apply_d_window_cuts(spec, tree, type, trigger):

        rdf_base = RDF(tree)

        if type == "MC" and "norm" not in spec:
            na = spec.split("_")
            tspec = f"{na[1]}_{na[2]}_{na[3]}"
            dwindow_file = ROOT.TFile(f"d_window_root_files/d_{tspec}_mass_fits.root", "READ")
            dwindow_ws = dwindow_file.Get(f"d_{tspec}_mass_fits")
            dwindow_cut = get_dwindow_values(dwindow_ws, tspec)
            dwindow_file.Close()
        if type == "DATA":
            dwindow_file = ROOT.TFile(f"d_window_root_files/d_{spec}_mass_fits.root", "READ")
            dwindow_ws = dwindow_file.Get(f"d_{spec}_mass_fits")
            dwindow_cut = get_dwindow_values(dwindow_ws, spec)
            dwindow_file.Close()
        if "norm" in spec:
            tspec = spec.split("_")[1]
            dwindow_file = ROOT.TFile(f"d_window_root_files/d_{tspec}_mass_fits.root", "READ")
            dwindow_ws = dwindow_file.Get(f"d_{tspec}_mass_fits")
            dwindow_cut = get_dwindow_values(dwindow_ws, tspec)
            dwindow_file.Close()

        rdf_post_dcuts = rdf_base.Filter(dwindow_cut)
        clist = rdf_post_dcuts.GetColumnNames()

        if type == "MC":
            outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/post_d/{spec}_{trigger}_postdcuts.root"
        if type == "DATA":
            outputfile = f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/{year}/post_d/{spec}_postdcuts.root"

        if os.path.exists(outputfile):
            os.remove(outputfile)
            print("First deleting old ", outputfile)
        else:
            print("making ", outputfile)

        rdfsnap = rdf_post_dcuts.Snapshot(f"DecayTreeTuple", outputfile, clist)

for c in [z_c]:

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/DATA_2021/*/pre_d/{c.spec}.root")
    tree_chain = ROOT.TChain(f"DecayTreeTuple_{c.spec}")
    for file in file_list:
        tree_chain.Add(file)

    build_d_window_ws(c, tree_chain)
    # plot_d_window_fits(c)
    #
    # for year in ["2016","2017","2018"]:
    #     file = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/DATA_2021/{year}/pre_d/{c.spec}.root")
    #     tree = file.Get(f"DecayTreeTuple_{c.spec}")
    #     apply_d_window_cuts(c.spec, tree, "DATA", "ToT")

# mc_spec_list = [
#     "01_Z_m_p_11198006",
#     "02_Z_m_p_11198400",
#     "04_Z_m_p_11198401",
#     "norm8_norm8_11198007",
# ]
# for name in mc_spec_list:
#     for year in ["2016","2017","2018"]:
#         for trigger in ["T","nTaT"]:
#
#             file = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/MC_2021/{year}/pre_d/{name}_{trigger}.root")
#             tree = file.Get(f"DecayTreeTuple_{name}")
#
#             apply_d_window_cuts(name, tree, "MC", trigger)



# build_d_window_ws(zz_c)
# build_d_window_ws(p_c)
# build_d_window_ws(m_c)
# build_d_window_ws(st_c)
# build_d_window_ws(s_c)
# build_d_window_ws(n7_c)
# build_d_window_ws(n8_c)
