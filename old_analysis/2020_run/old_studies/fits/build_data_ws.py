import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
p_data_file = ROOT.TFile(data_basepath+"p_spectrum.root")
m_data_file = ROOT.TFile(data_basepath+"m_spectrum.root")
st_data_file = ROOT.TFile(data_basepath+"st_spectrum_newfd.root")
s_data_file = ROOT.TFile(data_basepath+"s_spectrum.root")

z_tree = z_data_file.Get("DecayTreeTuple")
zz_tree = zz_data_file.Get("DecayTreeTuple")
p_tree = p_data_file.Get("DecayTreeTuple")
m_tree = m_data_file.Get("DecayTreeTuple")
st_tree = st_data_file.Get("DecayTreeTuple")
s_tree = s_data_file.Get("DecayTreeTuple")

dws = build_data_ws(data_ws_name)

b_dtf_m = dws.var("B_DTF_M")
b_dtf_m.setRange("myrange", bmin, bmax)
d1_mass = dws.var("D1_M")
d2_mass = dws.var("D2_M")
d1_dira = dws.var("D1_DIRA_ORIVX")
d2_dira = dws.var("D2_DIRA_ORIVX")
d1_fdx2 = dws.var("D1_FDCHI2_ORIVX")
d2_fdx2 = dws.var("D2_FDCHI2_ORIVX")
d2st_mass = dws.var("D2st_M")

data_args = ROOT.RooArgSet(d1_mass, d2_mass, d1_dira, d2_dira, d1_fdx2, d2_fdx2, b_dtf_m)
st_data_args = ROOT.RooArgSet(d1_mass, d2_mass, d1_dira, d2_dira, d1_fdx2, d2_fdx2, b_dtf_m, d2st_mass)

all_cats = ROOT.RooCategory("all_cats","all_cats")
all_cats.defineType("z_spectrum")
all_cats.defineType("zz_spectrum")
all_cats.defineType("p_spectrum")
all_cats.defineType("m_spectrum")
all_cats.defineType("st_spectrum")
all_cats.defineType("s_spectrum")

z_data_set = ROOT.RooDataSet("z_data","z_data", z_tree, data_args, z_m_cut)
zz_data_set = ROOT.RooDataSet("zz_data","zz_data", zz_tree, data_args, zz_m_cut)
p_data_set = ROOT.RooDataSet("p_data","p_data", p_tree, data_args, p_m_cut)
m_data_set = ROOT.RooDataSet("m_data","m_data", m_tree, data_args, m_m_cut)
st_data_set = ROOT.RooDataSet("st_data","st_data", st_tree, st_data_args, st_m_cut)
s_data_set = ROOT.RooDataSet("s_data","s_data", s_tree, data_args, s_m_cut)

# else:
#     z_data_set = ROOT.RooDataSet("z_data","z_data", z_tree, data_args, test_z_cut)
#     zz_data_set = ROOT.RooDataSet("zz_data","zz_data", zz_tree, data_args, test_zz_cut)
#     p_data_set = ROOT.RooDataSet("p_data","p_data", p_tree, data_args, test_p_cut)
#     m_data_set = ROOT.RooDataSet("m_data","m_data", m_tree, data_args, test_m_cut)
#     st_data_set = ROOT.RooDataSet("st_data","st_data", st_tree, data_args, test_st_cut)

all_data_sets = ROOT.RooDataSet("all_data_sets", "all_data_sets", data_args, ROOT.RooFit.Index(all_cats),
            ROOT.RooFit.Import(
                        "z_spectrum",
                        z_data_set),
            ROOT.RooFit.Import(
                        "zz_spectrum",
                        zz_data_set),
            ROOT.RooFit.Import(
                        "p_spectrum",
                        p_data_set),
            ROOT.RooFit.Import(
                        "m_spectrum",
                        m_data_set),
            ROOT.RooFit.Import(
                        "st_spectrum",
                        st_data_set),
            ROOT.RooFit.Import(
                        "s_spectrum",
                        s_data_set),
               )

getattr(dws,'import')(all_cats)
getattr(dws,'import')(all_data_sets)

dws.writeToFile(f"{data_ws_basepath}/{data_ws_name}.root")
dws.Print()
