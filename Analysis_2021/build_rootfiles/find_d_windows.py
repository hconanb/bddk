import sys
from createXFD import *
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *
from rootutils import *

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

def plot_d_window_fits(part):

    if part == "mp":
        fstring = "D^{-} / D^{+}"
    if part == "z":
        fstring = "#bar{D^{0}} / D^{0}"
    if part == "sm":
        fstring = "D_{s}^{-} / D_{s}^{+}"
    if part == "dst":
        fstring = "D^{*+} - D^{0} / D^{*-} - #bar{D^{0}}"
    if part == "d0k3pi":
        fstring = "#bar{D^{0}} #rightarrow K^{+}#pi^{+}#pi^{-}#pi^{-}"

    file = ROOT.TFile(f"d_window_root_files/d_{part}_mass_fits.root")
    dws = file.Get(f"d_{part}_mass_fits")
    var = dws.var("D1_M")

    if part == "dst" or part == "d0k3pi":
        all_data_sets = dws.data("data_d2")
    else:
        all_data_sets = dws.data("data_d1")

    all_fit = dws.obj(f"D_fit")

    frame = var.frame()

    ###Add ROOT.RooFit.Binning(50)
    all_data_sets.plotOn(frame, ROOT.RooFit.Name(f"{part}_data"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D_fit"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"{part}_fit"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D_signal"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"{part}_signal"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D_BKG"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name(f"BKG"))


    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D1_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"{part}_ga"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D1_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"{part}_gb"))

    d_mean = dws.obj(f"mean")
    d_stda = dws.var(f"width_a")

    d_stdb = dws.var(f"width_b")
    d_sigfrac = dws.var(f"D1_a_frac")

    d_nyield = dws.var(f"n_D_signal")
    d_nbkg = dws.var(f"n_D_bkg")

    # d_mean = dws.obj(f"mean")
    # d_nyield = dws.var(f"n_D_signal")
    # d_stda = dws.var(f"width_a")
    # d_nbkg = dws.var(f"n_D_bkg")
    d_chi2 = frame.chiSquare(f"{part}_fit", f"{part}_data")

    dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.AddText(f"mean: {round(d_mean.getValV(),3)}")
    # if part != "dst":
    dtpave.AddText(f"nsig: {round(d_nyield.getValV(),3)}")
    dtpave.AddText(f"nsig_a: {round(d_nyield.getValV()*d_sigfrac.getValV(),3)}")
    dtpave.AddText(f"nsig_b: {round(d_nyield.getValV()*(1-d_sigfrac.getValV()),3)}")

    # dtpave.AddText(f"nsig_a: {round(d_nyield.getValV(),3)}")
    dtpave.AddText(f"width_a: {round(d_stda.getValV(),3)}")

    dtpave.AddText(f"width_b: {round(d_stdb.getValV(),3)}")
    dtpave.AddText(f"nbkg: {round(d_nbkg.getValV(),3)}")
    frame.addObject(dtpave)

    p = residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)

    leg = ROOT.TLegend(0.75,0.5,0.95,0.87)
    leg.AddEntry(f"{part}_data","Data", "P");
    leg.AddEntry(frame.findObject(f"{part}_fit"),"Signal + background","L");
    leg.AddEntry(frame.findObject(f"BKG"),"Background only", "L");
    leg.AddEntry(frame.findObject(f"{part}_signal"),"Signal only", "L");
    # if part != "dst":
    leg.AddEntry(frame.findObject(f"{part}_ga"),"Signal A", "L");
    leg.AddEntry(frame.findObject(f"{part}_gb"),"Signal B", "L");
    frame.Draw()
    leg.Draw()
    p.pb.cd()

    hpull = frame.pullHist(f"{part}_data", f"{part}_fit")
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


    pull.GetXaxis().SetTitle(f"{fstring} Mass [MeV]")
    save_png(p, "d_window_plots", f"{part}", rpflag=1)

def build_d_window_ws(tree_chain_d1, tree_chain_d2, part):

    dname = f"d_{part}_mass_fits"
    dws = ROOT.RooWorkspace(dname)
    dfithalf = 40

    if part == "mp":
        d1_mass = dpmass
        d2_mass = dpmass
    if part == "z" or part == "d0k3pi":
        d1_mass = d0mass
        d2_mass = d0mass
    if part == "sm":
        d1_mass = dsmass
        d2_mass = dsmass
    if part == "dst":
        d1_mass = 145
        d2_mass = 145
        dfithalf = 10
        dws.factory(f"D2stmD_M[{d1_mass - dfithalf},{d1_mass + dfithalf}]")
        dst_mass_var = dws.var("D2stmD_M")

    dws.factory(f"D1_M[{d1_mass - dfithalf},{d1_mass + dfithalf}]")
    dws.factory(f"D2_M[{d2_mass - dfithalf},{d2_mass + dfithalf}]")

    dws.factory(f"Exponential:D_BKG(D1_M, c0_Dz[0, -1, 1])")

    # if part == "dst":
    #     #dws.factory(f"Gaussian::D1_signal(D1_M, mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[15.0,0.0,20.0])")
    #     dws.factory(f"BifurGauss::D1_signal(D1_M, mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[15.0,0.0,20.0], width_b[15.0,0.0,20.0])")
    # if part != "dst":
    dws.factory(f"Gaussian::D1_a(D1_M, mean[{d1_mass},{d1_mass-15},{d1_mass+15}], width_a[15.0,0.1,30.0])")
    dws.factory(f"Gaussian::D1_b(D1_M, mean, width_b[5.0,0.1,30.0])")
    dws.factory(f"SUM::D_signal(D1_a_frac[0.5,0,1]*D1_a, D1_b)")
    # if part == "z" or part == "d0k3pi":
        # dws.factory(f"Gaussian::D_signal(D1_M, mean[{d1_mass},{d1_mass-15},{d1_mass+15}], width_a[15.0,0.1,30.0])")
        # # dws.factory(f"Gaussian::D1_b(D1_M, mean, width_b[5.0,0.1,30.0])")
        # dws.factory(f"SUM::D1_signal(D1_a_frac[0.5,0,1]*D1_a, D1_b)")
    dws.factory(f"SUM::D_fit(n_D_signal[1000,0,100000]*D_signal, n_D_bkg[1000,0,100000]*D_BKG)")

    d1_mass_var = dws.var("D1_M")
    d2_mass_var = dws.var("D2_M")

    all_fit = dws.pdf("D_fit")

    data_args_d1 = ROOT.RooArgSet(d1_mass_var)
    data_args_d2 = ROOT.RooArgSet(d2_mass_var)

    if part == "mp" or part == "z" or part == "sm":
        data_set_d1 = ROOT.RooDataSet(f"data_d1", f"data_d1", tree_chain_d1, data_args_d1)
    if part == "mp" or part == "z" or part == "d0k3pi":
        data_set_d2 = ROOT.RooDataSet(f"data_d2", f"data_d2", tree_chain_d2, data_args_d2)
        dws.Import(data_set_d2, ROOT.RooFit.RenameVariable("D2_M", "D1_M"))
        data_set_d2_fix = dws.data("data_d2")
    if part == "dst":
        data_args_dst = ROOT.RooArgSet(dst_mass_var)
        data_set_d2 = ROOT.RooDataSet(f"data_d2", f"data_d2", tree_chain_d2, data_args_dst)
        dws.Import(data_set_d2, ROOT.RooFit.RenameVariable("D2stmD_M", "D1_M"))
        data_set_d2_fix = dws.data("data_d2")
    if part == "mp" or part == "z":
        data_set_d1.append(data_set_d2_fix)
    if part == "dst" or part == "d0k3pi":
        data_set_d1 = data_set_d2_fix

    all_fit.fitTo(data_set_d1, ROOT.RooFit.PrintLevel(0))
    dws.Import(data_set_d1)
    dws.Import(all_fit)
    dws.writeToFile(f"d_window_root_files/{dname}.root")

def get_files(type, loc, spec):
    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_{type}/*/{loc}/{spec}.root")
    return file_list

zmp = get_files("data", "pre_d", "Z_m_p")
zzz = get_files("data", "pre_d", "Z_z_z")
pzp = get_files("data", "pre_d", "P_z_p")
mmz = get_files("data", "pre_d", "M_m_z")
pzpst = get_files("data", "pre_d", "P_z_pst")
zssmp = get_files("data", "pre_d", "Zs_sm_p")
n7 = get_files("data", "pre_d", "norm7")
n8 = get_files("data", "pre_d", "norm8")

d1_mp_file_list = zmp + mmz + n8
d2_mp_file_list = zmp + pzp + zssmp

d1_z_file_list = zzz + pzp + pzpst + n7
d2_z_file_list = zzz + mmz + pzpst

d1_sm_file_list = zssmp

d2_dst_file_list = pzpst

d2_d0k3pi_file_list = n7 + n8

#"mp","z", "sm", "dst",
#"z","d0k3pi"
for d_flag in ["mp",]:

    tree_chain_d1 = ROOT.TChain(f"DecayTreeTuple")
    tree_chain_d2 = ROOT.TChain(f"DecayTreeTuple")

    if d_flag == "mp":
        for file in d1_mp_file_list:
            tree_chain_d1.Add(file)
        for file in d2_mp_file_list:
            tree_chain_d2.Add(file)

    if d_flag == "z":
        for file in d1_z_file_list:
            tree_chain_d1.Add(file)
        for file in d2_z_file_list:
            tree_chain_d2.Add(file)

    if d_flag == "sm":
        for file in d1_sm_file_list:
            tree_chain_d1.Add(file)

    if d_flag == "dst":
        for file in d2_dst_file_list:
            tree_chain_d2.Add(file)

    if d_flag == "d0k3pi":
        for file in d2_d0k3pi_file_list:
            tree_chain_d2.Add(file)

    build_d_window_ws(tree_chain_d1, tree_chain_d2, d_flag)
    plot_d_window_fits(d_flag)
