import sys
import os

basedir = os.getcwd().split('fits')[0]
sys.path.append(basedir)

from essential_functions import *
from rootutils import residualPlot

RDF = ROOT.ROOT.RDataFrame

def plot_d_window_fits(part, strat, dfithalf):

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

    file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/d_window_files//d_{part}_mass_fits.root")
    dws = file.Get(f"d_{part}_mass_fits")
    var = dws.var("D1_M")

    if part == "dst" or part == "d0k3pi":
        all_data_sets = dws.data("data_d2")
    else:
        all_data_sets = dws.data("data_d1")

    all_fit = dws.obj(f"D_fit")

    frame = var.frame(ROOT.RooFit.Bins(100))

    ###Add ROOT.RooFit.Binning(50)
    all_data_sets.plotOn(frame, ROOT.RooFit.Name(f"{part}_data"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D_fit"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"{part}_fit"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D_signal"), ROOT.RooFit.LineStyle(ROOT.kSolid), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"{part}_signal"))
    all_fit.plotOn(frame, ROOT.RooFit.Components(f"D_BKG"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name(f"BKG"))

    if strat == "DG" or strat == "cb2":
        all_fit.plotOn(frame, ROOT.RooFit.Components(f"D1_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kOrange), ROOT.RooFit.Name(f"{part}_ga"))
        all_fit.plotOn(frame, ROOT.RooFit.Components(f"D1_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kMagenta), ROOT.RooFit.Name(f"{part}_gb"))

    d_mean = dws.obj(f"mean")
    d_stda = dws.var(f"width_a")

    d_stdb = dws.var(f"width_b")
    d_sigfrac = dws.var(f"D1_a_frac")

    d_nyield = dws.var(f"n_D_signal")
    d_nbkg = dws.var(f"n_D_bkg")

    d_chi2 = frame.chiSquare(f"{part}_fit", f"{part}_data")

    dtpave = ROOT.TPaveText(0.125, 0.55, 0.525, 0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.SetBorderSize(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.AddText(f"Mean: {round(d_mean.getValV(),3)}")
    if strat == "G" or strat == "BGG" or strat == "cb1L" or strat == "cb1R":
        dtpave.AddText(f"Signal Candidates: {round(d_nyield.getValV(),0)}")
    if strat == "DG" or strat == "cb2":
        dtpave.AddText(f"nsig_a: {round(d_nyield.getValV()*d_sigfrac.getValV(),0)}")
        dtpave.AddText(f"nsig_b: {round(d_nyield.getValV()*(1-d_sigfrac.getValV()),0)}")
    dtpave.AddText(f"Width: {round(d_stda.getValV(),3)}")
    if strat == "DG" or strat == "BGG" or strat == "cb2":
        dtpave.AddText(f"Width 2: {round(d_stdb.getValV(),3)}")
    if "cb" in strat:
        d_alpha = dws.var("alpha_a")
        dtpave.AddText(f"Alpha: {round(d_alpha.getValV(),2)}")
    if strat == "cb2":
        d_alpha_b = dws.var("alpha_b")
        dtpave.AddText(f"Alpha 2: {round(d_alpha_b.getValV(),2)}")
    dtpave.AddText(f"Background Candidates: {round(d_nbkg.getValV(),0)}")
    dtpave.SetTextSize(0.04)
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
    if strat == "DG" or strat =="cb2" :
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
    save_png(p, "d_window_plots", f"{part}_{strat}_{dfithalf}", rpflag=1)

    print("Getting interval")
    sig_pdf = dws.pdf("D_signal")
    distart = d_mean.getValV()
    diw = round(d_stda.getValV(),0)
    rlist = [3.3*diw, 3.2*diw, 3.1*diw, 3*diw, 2.9*diw, 2.8*diw, 2.7*diw, 2.6*diw, 2.5*diw, 2.4*diw, 2.3*diw, 2.2*diw, 2.1*diw, 2.0*diw, 1.9*diw, 1.8*diw, 1.7*diw, 1.6*diw, 1.5*diw, 1.4*diw, 1.3*diw]
    nset = ROOT.RooArgSet(var)
    normSet = ROOT.RooFit.NormSet(nset)


    f = open(f'/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/dwindow_txt/{part}_{strat}','w')
    for r in rlist:
        var.setRange(f"signal_{r}", distart - r , distart + r)
        igx_sig = sig_pdf.createIntegral(nset, normSet, ROOT.RooFit.Range(f"signal_{r}"))
        print(f"With {r} at gx_Int[x|signal]_Norm[x] = ", igx_sig.getVal())
        f.write(f"With {r} at gx_Int[x|signal]_Norm[x] = {igx_sig.getVal()}\n")
    f.close()

def build_d_window_ws(tree_chain_d1, tree_chain_d2, part, strat, dfithalf):

    dname = f"d_{part}_mass_fits"
    dws = ROOT.RooWorkspace(dname)

    if part == "mp":
        d1_mass = dpmass
        d2_mass = dpmass
    if part == "z" or part == "d0k3pi":
        d1_mass = d0mass
        d2_mass = d0mass
    if part == "dst":
        d1_mass = 145.5
        d2_mass = 145.5
        dfithalf = 4
        dws.factory(f"D2stmD_M[{d1_mass - dfithalf},{d1_mass + dfithalf}]")
        dst_mass_var = dws.var("D2stmD_M")

    dws.factory(f"D1_M[{d1_mass - dfithalf},{d1_mass + dfithalf}]")
    dws.factory(f"D2_M[{d2_mass - dfithalf},{d2_mass + dfithalf}]")

    dws.factory(f"Exponential:D_BKG(D1_M, c0_Dz[0, -2, 2])")
    if strat == "G":
        dws.factory(f"Gaussian::D_signal(D1_M, mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[15.0,0.1,30.0])")
    if strat == "DG":
        dws.factory(f"Gaussian::D1_a(D1_M, mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[10,0.1,30])")
        dws.factory(f"Gaussian::D1_b(D1_M, mean, width_b[6.0,0.1,25])")
        dws.factory(f"SUM::D_signal(D1_a_frac[0.5, 0.1, 1]*D1_a, D1_b)")
    if strat == "BGG":
        dws.factory(f"BifurGauss::D_signal(D1_M,mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[15.0,0.1,30.0],width_b[14.0,0.1,30.0])")
    if strat == "cb1L":
        dws.factory(f"CBShape::D_signal(D1_M,mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[10,0.1,30], alpha_a[1,0.01,3], n_a[50,40,60])")
    if strat == "cb1R":
        dws.factory(f"CBShape::D_signal(D1_M,mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[10,0.1,30], alpha_a[-1,-3.0,-0.01], n_a[150,40,60])")
    if strat == "cb2":
        dws.factory(f"CBShape::D1_a(D1_M,mean[{d1_mass},{d1_mass-10},{d1_mass+10}], width_a[10,0.1,30], alpha_a[1,0.01,3], n_a[50,40,60])")
        dws.factory(f"CBShape::D1_b(D1_M,mean, width_b[10,0.1,30], alpha_b[-1,-3.0,-0.01], n_b[50,40,60])")
        dws.factory(f"SUM::D_signal(D1_a_frac[0.5, 0.1, 1]*D1_a, D1_b)")

    dws.factory(f"SUM::D_fit(n_D_signal[10000,0,1000000]*D_signal, n_D_bkg[01000,0,1000000]*D_BKG)")

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
    if part == "d0k3pi" or part == "dst":
        data_set_d1 = data_set_d2_fix

    all_fit.fitTo(data_set_d1, ROOT.RooFit.PrintLevel(0))

    # Create an integral of gx_Norm[x] over x in range "signal"
    # ROOT.This is the fraction of of pdf gx_Norm[x] which is in the
    # range named "signal"
    # int1 = s1s2.createIntegral(ROOT.RooArgSet(m1, m2), ROOT.RooFit.NormSet(ROOT.RooArgSet(m1, m2)), ROOT.RooFit.Range('sig'))

    if part != "d0k3pi" and part != "dst":
        dws.Import(data_set_d1)
        dws.Import(all_fit)

    dws.writeToFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/d_window_files/{dname}.root")

def get_files(spec):
    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/post_trigger/{spec}.root")
    return file_list

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Find D Window Cuts")
    parser.add_argument("--D_Flag_List", choices = ["mp","z","dst","d0k3pi"],  nargs="+", help = 'D_Flag')
    parser.add_argument("--Strat", help = 'Strat')
    parser.add_argument('--FIT', action='store_true')
    parser.add_argument('--PLOT', action='store_true')

    args = parser.parse_args()
    d_flag_list = args.D_Flag_List
    strat = args.Strat
    fit_flag = args.FIT
    plot_flag = args.PLOT

    zmp = get_files("Z_m_p")
    zzz = get_files("Z_z_z")
    pzp = get_files("P_z_p")
    mmz = get_files("M_m_z")
    pzpst = get_files("P_z_pst")
    n7 = get_files("norm7")
    n8 = get_files("norm8")

    d1_mp_file_list = zmp
    d2_mp_file_list = zmp 
    d1_z_file_list = zzz + pzp + pzpst + n7
    d2_z_file_list = zzz + mmz + pzpst
    d2_d0k3pi_file_list = n7 + n8
    d2_dst_file_list = pzpst

    tree_chain_d1 = ROOT.TChain(f"DecayTreeTuple")
    tree_chain_d2 = ROOT.TChain(f"DecayTreeTuple")

    for d_flag in d_flag_list:

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

        if d_flag == "dst":
            for file in d2_dst_file_list:
                tree_chain_d2.Add(file)

        if d_flag == "d0k3pi":
            for file in d2_d0k3pi_file_list:
                tree_chain_d2.Add(file)

        dfithalf = 50

        if fit_flag:
            build_d_window_ws(tree_chain_d1, tree_chain_d2, d_flag, strat, dfithalf)
        if plot_flag:
            plot_d_window_fits(d_flag, strat, dfithalf)
