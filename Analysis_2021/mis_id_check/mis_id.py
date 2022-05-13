import sys
sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/")
from essentials import *


class spectrum_class:

  def __init__(self, name, d1_string, d2_string, d1_mass, d2_mass, b_string):
    self.spec = name

    self.d1_string = d1_string
    self.d2_string = d2_string

    self.d1_mass = d1_mass
    self.d2_mass = d2_mass

    if self.spec == "st":
        self.d3_string = "(D^{*+} #rightarrow D^{0} #pi+)"
        self.d3_mass = 150

    self.b_string = b_string

z_c = spectrum_class("Z_m_p", "D^{-}", "D^{+}", dpmass, dpmass, "B #rightarrow D-D+K*0")
zz_c = spectrum_class("Z_z_z", "#barD^{0}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D^{0} K*0")
p_c = spectrum_class("P_z_p", "#bar{D^{0}}", "D^{+}", d0mass, dpmass, "B #rightarrow #bar{D^{0}} D+ K*0")
m_c = spectrum_class("M_m_z", "D^{-}", "D^{0}", dpmass, d0mass,  "B #rightarrow D- D^{0} K*0")
st_c = spectrum_class("P_z_pst","#bar{D^{0}}", "D^{0}", d0mass, d0mass, "B #rightarrow #bar{D^{0}} D*+ K*0")


def misid_check(sc):
    spec = sc.spec

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/post_d/{spec}.root")

    if spec == "M_m_z":
        d1_flag = "mp"
        d2_flag = "z"
    tc = ROOT.TChain(f"DecayTreeTuple")
    for file in file_list:
        tc.Add(file)

    rdf_base = RDF(tc)

    D1H3_PE_new = "sqrt((pow(493.677,2)) + pow(D1H3_PX,2) + pow(D1H3_PY,2) + pow(D1H3_PZ,2))"
    PE2_d1h1d1h2d1h3_new = "(pow(D1H1_PE + D1H2_PE + D1H3_PE_new,2))"
    PX2_d1h1d1h2d1h3 = "(pow(D1H1_PX + D1H2_PX + D1H3_PX,2))"
    PY2_d1h1d1h2d1h3 = "(pow(D1H1_PY + D1H2_PY + D1H3_PY,2))"
    PZ2_d1h1d1h2d1h3 = "(pow(D1H1_PZ + D1H2_PZ + D1H3_PZ,2))"
    M2_d1h1d1h2d1h3 = f"{PE2_d1h1d1h2d1h3_new} - ({PX2_d1h1d1h2d1h3} + {PY2_d1h1d1h2d1h3} + {PZ2_d1h1d1h2d1h3})"
    D1_M_new = f"sqrt({M2_d1h1d1h2d1h3})"

    rdf_data_base = rdf_base.Define("D1H3_PE_new", D1H3_PE_new) \
                            .Define("PE_d1h1d1h2d1h3_new", PE2_d1h1d1h2d1h3_new) \
                            .Define("D1_M_new", D1_M_new)

    clist = rdf_data_base.GetColumnNames()
    outputfile = f"root_files/{spec}_mid_check.root"
    print(f"Starting snapshot for {outputfile}")

    rdfsnap = rdf_data_base.Snapshot(f"{spec}_mid", outputfile, clist)

for sc in [m_c]:
    misid_check(sc)
    # fit_est(sc, strat = "2peak")
    # fit_est(sc, strat = "3peak")
