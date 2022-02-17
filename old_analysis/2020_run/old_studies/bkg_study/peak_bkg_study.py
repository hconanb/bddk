from bkg_essentials import *
from itertools import combinations
RDF = ROOT.ROOT.RDataFrame

tla = ROOT.TLatex()
tla.SetTextFont(32)
tla.SetTextColor(1)
tla.SetTextSize(0.03)
tla.SetTextAlign(12)

def build_track_strings(tl, n, tag):
    name_list = []
    form_list = []
    if n == 2:
        for i in tl :
            t1 = i[0]
            t2 = i[1]
            name = f"{t1}_{t2}_M_{tag}"
            form = f"sqrt( pow({t1}_PE + {t2}_PE, 2) - ( pow({t1}_PX + {t2}_PX,2) + pow({t1}_PY + {t2}_PY, 2) + pow({t1}_PZ + {t2}_PZ,2)))"
            name_list.append(name)
            form_list.append(form)
    if n == 3:
        for i in tl :
            t1 = i[0]
            t2 = i[1]
            t3 = i[2]
            name = f"{t1}_{t2}_{t3}_M_{tag}"
            form = f"sqrt( pow({t1}_PE + {t2}_PE + {t3}_PE, 2) - ( pow({t1}_PX + {t2}_PX + {t3}_PX,2) + pow({t1}_PY + {t2}_PY + {t3}_PY, 2) + pow({t1}_PZ + {t2}_PZ + {t3}_PZ,2)))"
            name_list.append(name)
            form_list.append(form)
    return name_list, form_list;
def build_final_rdf(rdf, nl, fl):
    rdf_list = [rdf]
    for i in range(len(nl)):
        rdf_temp = rdf_list[i].Define(nl[i], fl[i])
        rdf_list.append(rdf_temp)
    return rdf_list[-1]
def plot_2(rdf, nl):
    for name in nl:
        if "D1H1" in name and "D1H2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 1700, 2100), name)
        elif "D1H1" in name and "D2H1" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 600, 2600), name)
        elif "D1H1" in name and "D2H2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 600, 3300), name)
        elif "D1H1" in name and "KSTH1" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 900, 2000), name)
        elif "D1H1" in name and "KSTH2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 300, 2700), name)

        elif "D1H2" in name and "D2H1" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 500, 3500), name)
        elif "D1H2" in name and "D2H2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 200, 3500), name)
        elif "D1H2" in name and "KSTH1" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 700, 2800), name)
        elif "D1H2" in name and "KSTH2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 200, 2400), name)

        elif "D2H1" in name and "D2H2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 1700, 2100), name)
        elif "D2H1" in name and "KSTH1" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 900, 3300), name)
        elif "D2H1" in name and "KSTH2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 200, 2300), name)

        elif "D2H2" in name and "KSTH1" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 600, 2000), name)
        elif "D2H2" in name and "KSTH2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 200, 1800), name)

        elif "KSTH1" in name and "KSTH2" in name:
            m2_hist = rdf.Histo1D((name, name, 100, 600, 1200), name)

        canvas = ROOT.TCanvas("c1","c1")
        newhist = m2_hist.DrawCopy()
        newhist.GetXaxis().SetTitle(name)
        canvas.Update()
        saveplot(canvas, name)
def plot_3(rdf, nl):
    for name in nl:
        m2_hist = rdf.Histo1D((name, name, 100, 1000, 3500), name)
        canvas = ROOT.TCanvas("c1","c1")
        newhist = m2_hist.DrawCopy()
        newhist.GetXaxis().SetTitle(name)
        canvas.Update()
        saveplot(canvas, name)

zz_data_file = ROOT.TFile(data_basepath+"zz_spectrum.root")
zz_tree = zz_data_file.Get("DecayTreeTuple")

z_data_file = ROOT.TFile(data_basepath+"z_spectrum.root")
z_tree = zz_data_file.Get("DecayTreeTuple")

Track_List_8 = {"D1H1","D1H2","D1H3","D2H1","D2H2","D2H3","KSTH1","KSTH2"}
Track_List_7 = {"D1H1","D1H2","D2H1","D2H2","D2H3","KSTH1","KSTH2"}
Track_List_6 = {"D1H1","D1H2","D2H1","D2H2","KSTH1","KSTH2"}

track6_2_name_list, track6_2_form_list = build_track_strings(combinations(Track_List_6, 2), 2, "zz")
track6_3_name_list, track6_3_form_list = build_track_strings(combinations(Track_List_6, 3), 3, "zz")

zz_base = RDF(zz_tree)
zz_base_pc = zz_base.Filter(f"(abs(D1_M - {d0mass}) < {dwindow}) && (abs(D2_M - {d0mass}) < {dwindow}) && D1_DIRA_ORIVX>0 && D2_DIRA_ORIVX>0")

rdf_zz_2 = build_final_rdf(zz_base, track6_2_name_list, track6_2_form_list)
rdf_zz_2_pc = build_final_rdf(zz_base_pc, track6_2_name_list, track6_2_form_list)

rdf_zz_3 = build_final_rdf(zz_tree, track6_3_name_list, track6_3_form_list)

# rdf_z_2 = build_final_rdf(z_tree, track8_2_name_list, track8_2_form_list)
#rdf_z_3 = build_final_rdf(z_tree, track8_3_name_list, track8_2_form_list)

#plot_2(rdf_zz_2, track6_2_name_list)
#plot_2(rdf_zz_2_pc, track6_2_name_list)
plot_3(rdf_zz_3, track6_3_name_list)
