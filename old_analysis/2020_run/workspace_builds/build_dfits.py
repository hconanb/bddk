import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
from rootutils import *

dstrat = "e_bg"
normdstrat = "dg"

def plot_d_fit(tag, dvar, dstring, fancy_dstring, frame):

    uid = f"{tag}_{dstring}"
    p = residualPlot()
    p.pt.cd()
    xaxis = frame.GetXaxis()
    xaxis.SetTickLength(0)
    xaxis.SetNdivisions(0)
    xaxis.SetLabelSize(0)
    frame.Draw()
    p.pb.cd()

    # legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
    # legend.SetHeader(title,"C")
    # legend.AddEntry("pdf","Total Fit","lp")
    # legend.AddEntry("fit0","Fully Reconstructed Peak","l")
    # legend.AddEntry("fit1","1 missing particle peak","l")
    # legend.AddEntry("fit2","2 missing particle peak","l")
    # legend.AddEntry("bkg","Background","l")
    # legend.Draw()

    hpull = frame.pullHist(f"{uid}_data", f"{uid}_fit")
    pull = dvar.frame(ROOT.RooFit.Range(f"{uid}_range"));
    pull.addPlotable(hpull, "P");
    pull.SetLineWidth(ROOT.gStyle.GetFrameLineWidth())
    pull.GetXaxis().SetLabelSize(ROOT.gStyle.GetLabelSize()*p.blabelratio)
    pull.GetYaxis().SetLabelSize(ROOT.gStyle.GetLabelSize("Y")*(1+p.padratio)/(2*p.padratio))
    pull.GetXaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
    pull.GetYaxis().SetTitleSize(ROOT.gStyle.GetTitleSize()*p.blabelratio)
    pull.GetYaxis().SetTitleOffset(ROOT.gStyle.GetTitleOffset()/p.blabelratio)
    pull.GetXaxis().SetTickLength(ROOT.gStyle.GetTickLength()*p.blabelratio)
    pull.GetYaxis().SetTitle("Pull")
    pull.GetYaxis().SetNdivisions(202,False)
    pull.GetYaxis().SetRangeUser(-3,3)
    pull.GetYaxis().CenterTitle()
    pull.Draw("AP")

    pull.GetXaxis().SetTitle(f"{fancy_dstring} Mass [MeV]")
    now = datetime.datetime.now()
    if not os.path.exists(f'plots/{now.month}_{now.day}/'):
        os.makedirs(f'plots/{now.month}_{now.day}/')
    p.save(f"plots/{now.month}_{now.day}/{uid}_{dstrat}.png")
def build_d_fit(tag, dvar, dstring, mydmass, data_set):
    uid = f"{tag}_{dstring}"
    dvar.setRange(f"{uid}_range", mydmass-50, mydmass+50)
    if "norm" not in tag:
        if dstrat == "e_g":
            dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -2, 2])")
            dws.factory(f"Gaussian::{uid}_signal({dstring}_M, mean_{uid}[{mydmass},{mydmass}-10,{mydmass}+10], width_{uid}[10.0,0.01,20])")
            dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal,n_{uid}_bkg[10000,0,1000000]*{uid}_bkg)")
        if dstrat == "e_bg":
            dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -2, 2])")
            dws.factory(f"BifurGauss::{uid}_signal({dstring}_M, mean_{uid}[{mydmass},{mydmass}-10,{mydmass}+10], width_{uid}_a[9.0,0.01,20], width_{uid}_b[10.0,0.01,20])")
            dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[1000,0,100000]*{uid}_signal,n_{uid}_bkg[10000,0,1000000]*{uid}_bkg)")
    else:
        if normdstrat == "dg":
            dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -3, 3])")
            dws.factory(f"Gaussian::{uid}_signal_a({dstring}_M, mean_{uid}[{mydmass},{mydmass}-10,{mydmass}+10], width_{uid}_a[20.0,0.01,40])")
            dws.factory(f"Gaussian::{uid}_signal_b({dstring}_M, mean_{uid}, width_{uid}_b[5.0,0.01,40])")
            dws.factory(f"SUM::{uid}_signal(signal_a_frac[0.5,0,1]*{uid}_signal_a, {uid}_signal_b)")
            if "norm7" in uid:
                print("ASDFASDF")
                dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[2000,0,1000000]*{uid}_signal,n_{uid}_bkg[5000000,0,100000000]*{uid}_bkg)")
            else:
                print("fajsdghf")
                dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[10000,0,1000000]*{uid}_signal,n_{uid}_bkg[100000,0,10000000]*{uid}_bkg)")
        if normdstrat == "cb":
            dws.factory(f"Exponential:{uid}_bkg({dstring}_M, c0_{uid}[0, -3, 3])")
            dws.factory(f"CBShape::{uid}_signal({dstring}_M,mean_{uid}[{mydmass},{mydmass}-10,{mydmass}+10],width_{uid}[10.0,0.01,40],alpha_{uid}[1,0.1,100], n_{uid}[3,1,100])")
            dws.factory(f"SUM::{uid}_fit(n_{uid}_signal[100000,0,1000000]*{uid}_signal,n_{uid}_bkg[100000,0,10000000]*{uid}_bkg)")


    print(f"built {uid}")

    d_model = dws.pdf(f"{uid}_fit")
    d_model.fitTo(data_set, ROOT.RooFit.Range(f"{uid}_range"), ROOT.RooFit.PrintLevel(0))
    d_frame = dvar.frame(ROOT.RooFit.Range(f"{uid}_range"))

    data_set.plotOn(d_frame, ROOT.RooFit.CutRange(f"{uid}_range"), ROOT.RooFit.Name(f"{uid}_data"))

    nData_d = data_set.sumEntries("", f"{uid}_range")

    d_model.plotOn(d_frame, ROOT.RooFit.NormRange(f"{uid}_range"), ROOT.RooFit.Normalization(nData_d, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Range(f"{uid}_range"), ROOT.RooFit.Components(f"{uid}_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name(f"{uid}_fit"))
    d_model.plotOn(d_frame, ROOT.RooFit.NormRange(f"{uid}_range"), ROOT.RooFit.Normalization(nData_d, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Range(f"{uid}_range"), ROOT.RooFit.Components(f"{uid}_signal"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"{uid}_signal"))
    d_model.plotOn(d_frame, ROOT.RooFit.NormRange(f"{uid}_range"), ROOT.RooFit.Normalization(nData_d, ROOT.RooAbsReal.NumEvent), ROOT.RooFit.Range(f"{uid}_range"), ROOT.RooFit.Components(f"{uid}_bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"{uid}_bkg"))

    # d_model.plotOn(d_frame, ROOT.RooFit.Components(f"{uid}_fit"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.Name(f"{uid}_fit"))
    # d_model.plotOn(d_frame, ROOT.RooFit.Components(f"{uid}_signal"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name(f"{uid}_signal"))
    # d_model.plotOn(d_frame, ROOT.RooFit.Components(f"{uid}_bkg"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"{uid}_bkg"))


    d_chi2 = d_frame.chiSquare(f"{uid}_fit", f"{uid}_data")
    d_mean = dws.obj(f"mean_{uid}")
    d_nyield = dws.var(f"n_{uid}_signal")
    d_nbkg = dws.var(f"n_{uid}_bkg")

    dtpave = ROOT.TPaveText(0.65,0.65,0.85,0.85, "NB NDC")
    dtpave.SetFillStyle(0)
    dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
    dtpave.AddText(f"mean: {d_mean.getValV()}")
    dtpave.AddText(f"nsig: {d_nyield.getValV()}")
    dtpave.AddText(f"nbkg: {d_nbkg.getValV()}")
    d_frame.addObject(dtpave)
    return d_frame

dname = "d_mass_fits"
dws = ROOT.RooWorkspace(dname)

dws.factory(f"B_DTF_M[{bmin},{bmax}]")
dws.factory("D1_M[0,6000]")
dws.factory("D2_M[0,6000]")
dws.factory("D1_DIRA_ORIVX[-5,5]")
dws.factory("D2_DIRA_ORIVX[-5,5]")
dws.factory("D2st_M[0,10000]")

b_dtf_m = dws.var("B_DTF_M")
d1_mass = dws.var("D1_M")
d2_mass = dws.var("D2_M")
d1_dira = dws.var("D1_DIRA_ORIVX")
d2_dira = dws.var("D2_DIRA_ORIVX")
d2st_mass = dws.var("D2st_M")

data_args = ROOT.RooArgSet(d1_mass, d2_mass, d1_dira, d2_dira, b_dtf_m)
st_data_args = ROOT.RooArgSet(d1_mass, d2_mass, d1_dira, d2_dira, b_dtf_m, d2st_mass)

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")
m_data_file = ROOT.TFile(data_basepath+"m_spectrum.root")
st_data_file = ROOT.TFile(data_basepath+"st_spectrum_newfd.root")
s_data_file = ROOT.TFile(data_basepath+"s_spectrum.root")
norm7_data_file  = ROOT.TFile(data_basepath+"norm7_data.root")
norm8_data_file  = ROOT.TFile(data_basepath+"norm8_data.root")

z_tree = z_data_file.Get("DecayTreeTuple")
zz_tree = zz_data_file.Get("DecayTreeTuple")
p_tree = p_data_file.Get("DecayTreeTuple")
m_tree = m_data_file.Get("DecayTreeTuple")
st_tree = st_data_file.Get("DecayTreeTuple")
s_tree = s_data_file.Get("DecayTreeTuple")
norm7_tree = norm7_data_file.Get("DecayTreeTuple")
norm8_tree = norm8_data_file.Get("DecayTreeTuple")

z_data_set = ROOT.RooDataSet("z_data","z_data", z_tree, data_args, dira_cut)
zz_data_set = ROOT.RooDataSet("zz_data","zz_data", zz_tree, data_args, dira_cut)
p_data_set = ROOT.RooDataSet("p_data","p_data", p_tree, data_args, dira_cut)
m_data_set = ROOT.RooDataSet("m_data","m_data", m_tree, data_args, dira_cut)
st_data_set = ROOT.RooDataSet("st_data","st_data", st_tree, st_data_args, dira_cut)
s_data_set = ROOT.RooDataSet("s_data","s_data", s_tree, data_args, dira_cut)
norm7_data_set = ROOT.RooDataSet("norm7_data","norm7_data", norm7_tree, data_args, dira_cut)
norm8_data_set = ROOT.RooDataSet("norm8_data","norm8_data", norm8_tree, data_args, dira_cut)

print("got data sets")
tag_list = ["z","zz","p","m","st","s","norm7","norm8"]
data_set_list = [z_data_set, zz_data_set, p_data_set, m_data_set, st_data_set, s_data_set, norm7_data_set, norm8_data_set]

for tag, data_set in zip(tag_list, data_set_list):
    if tag == "z":
        myd1mass = dpmass
        myd2mass = dpmass
        d1id = "D^{-}"
        d2id = "D^{+}"
    if tag == "zz":
        myd1mass = d0mass
        myd2mass = d0mass
        d1id = "#bar{D^{0}}"
        d2id = "D^{0}"
    if tag == "p":
        myd1mass = d0mass
        myd2mass = dpmass
        d1id = "#bar{D^{0}}"
        d2id = "D^{+}"
    if tag == "m":
        myd1mass = dpmass
        myd2mass = d0mass
        d1id = "D^{-}"
        d2id = "D^{0}"
    if tag == "st":
        myd1mass = d0mass
        myd2mass = d0mass
        d1id = "D^{0}"
        d2id = "(D^{*+} #rightarrow D^{0} #pi+)"
    if tag == "s":
        myd1mass = dsmass
        myd2mass = dpmass
        d1id = "D^{-}_s"
        d2id = "D^{+}"
    if tag == "norm7":
        myd1mass = d0mass
        myd2mass = d0mass
        d1id = "#bar{D^{0}}"
        d2id = "D^{0} #rightarrow k#pi#pi#pi"
    if tag == "norm8":
        myd1mass = dpmass
        myd2mass = d0mass
        d1id = "D^{-}"
        d2id = "D^{0} #rightarrow K#pi#pi#pi"

    d1f = build_d_fit(tag, d1_mass, "D1", myd1mass, data_set)
    plot_d_fit(tag, d1_mass, "D1", d1id, d1f)
    d2f = build_d_fit(tag, d2_mass, "D2", myd2mass, data_set)
    plot_d_fit(tag, d2_mass, "D2", d2id, d2f)

dws.writeToFile(f"dmass/{dname}.root")
dws.Print()
