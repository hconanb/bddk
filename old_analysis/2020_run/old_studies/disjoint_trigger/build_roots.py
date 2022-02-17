import os, sys
import array
import numpy as np
import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
#
print ("Max tree size",ROOT.TTree.GetMaxTreeSize())
ROOT.TTree.SetMaxTreeSize(2000000000000000000)
print ("Updated tree size",ROOT.TTree.GetMaxTreeSize())

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
#ROOT.ROOT.EnableImplicitMT()
opts.fMode = "UPDATE"

trigger_TOS = "B_L0Global_TOS == 1"
trigger_TIS = "(B_L0HadronDecision_TIS == 1 ||  B_L0MuonDecision_TIS == 1 ||  B_L0ElectronDecision_TIS == 1 ||  B_L0PhotonDecision_TIS == 1)"
trigger_HLT1 = "B_Hlt1Global_TOS == 1"
trigger_HLT2 = "B_Hlt2Global_TOS == 1"

fl_dict = {
#"ToT" : f"({trigger_TOS} || {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
"nTaT" : f"(!{trigger_TOS} && {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
"T" : f"({trigger_TOS}) && {trigger_HLT1} && {trigger_HLT2}",

#"norm_ToT" : f"({trigger_TOS} || {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2} && {loose_filters_norm}",
# "norm_nTaT" : f"(!{trigger_TOS} && {trigger_TIS}) && {trigger_HLT1} && {trigger_HLT2}",
# "norm_T" : f"({trigger_TOS}) && {trigger_HLT1} && {trigger_HLT2}",
}

def mc_filter(id, spec, eid, flines = ["nTaT", "T"]):

    thisname = f"{id}_{spec}_{eid}"

    inpath = f"/data6/hbernste/mc/{eid}"
    outputfile = f"/data6/hbernste/filtered/mc/{thisname}.root"

    if os.path.exists(outputfile):
        os.remove(outputfile)
        print("First deleting old ", outputfile)

    print ('Now building new ', outputfile)

    file_list = []

    tc = ROOT.TChain(f"{thisname}/DecayTreeTuple")
    tcmc = ROOT.TChain("MCDecayTreeTuple/MCDecayTree")

    for path, dirs, files in os.walk(inpath):
        for filename in files:
            if ('ntuple.root') in filename:
                fp = path+"/"+filename
                t = ROOT.TFile.Open(fp)
                if (t.GetListOfKeys().Contains(thisname) and t.GetListOfKeys().Contains('MCDecayTreeTuple')):
                    tc.Add(fp)
                    tcmc.Add(fp)
                else:
                    print (f"Not adding {fp} to tc or tcmc")
                t.Close()

    rdf = RDF(tc)
    rdf_MC = RDF(tcmc)


    TRUE_PE_REC = "(pow(D1_TRUEP_E + D2_TRUEP_E + KST_TRUEP_E,2))"
    TRUE_PX_REC = "(pow(D1_TRUEP_X + D2_TRUEP_X + KST_TRUEP_X,2))"
    TRUE_PY_REC = "(pow(D1_TRUEP_Y + D2_TRUEP_Y + KST_TRUEP_Y,2))"
    TRUE_PZ_REC = "(pow(D1_TRUEP_Z + D2_TRUEP_Z + KST_TRUEP_Z,2))"
    TRUE_BM_REC = f"sqrt({TRUE_PE_REC} - ({TRUE_PX_REC} + {TRUE_PY_REC} + {TRUE_PZ_REC}))"

    if "norm" not in thisname:
        rdf_dtf = rdf.Define("B_DTF_M",'B_dtf_M[0]') \
                     .Define("TRUE_BM_REC", TRUE_BM_REC) \
                     .Define("Res",f"B_DTF_M - TRUE_BM_REC")


    if "norm" in thisname:
        rdf_dtf = rdf.Define("B_DTF_M",'B_dtf_M[0]')

    print(f"Building MC DecayTreeTuple\n")
    mc_clist = rdf_MC.GetColumnNames()
    rdfsnap_MC = rdf_MC.Snapshot(f"MCDecayTreeTuple", outputfile, mc_clist, opts)
    print(f"Added MC DecayTreeTuple\n")

    for flag in flines:
        print(f"Adding DecayTreeTuple_{flag}\n")
        rdf_clean = rdf_dtf.Filter(fl_dict[flag])
        clist = rdf_clean.GetColumnNames()
        rdfsnap = rdf_clean.Snapshot(f"DecayTreeTuple_{flag}", outputfile, clist, opts)
        print(f"Added DecayTreeTuple_{flag}\n")

    print(f"Created file {outputfile}\n")
def data_filter(inflag, flines = ["nTaT", "T"]):

    inputfile = f"{data_basepath}/{inflag}"



    inputfile = ROOT.TFile(f"{data_basepath}/{inflag}")
    tree = inputfile.Get("DecayTreeTuple")
    outputfile = f"{TT_data_basepath}"+"TT_{i}"

    print("got all files")

    rdf = RDF(tree)

    print ('Now building ', outputfile)
    for flag in flines:
        print(f"the flag is {flag}\n")
        print(f"Adding DecayTreeTuple_{flag}\n")
        rdf_tt = rdf.Filter(fl_dict[flag])
        clist = rdf_tt.GetColumnNames()
        rdfsnap = rdf_tt.Snapshot(f"DecayTreeTuple_{flag}", outputfile, clist, opts)
        print(f"Added DecayTreeTuple_{flag}\n")

    print(f"Created file {outputfile}\n")

# mc_filter("norm7","12197009", ["norm7"])
# mc_filter("norm8","11198030",["norm8"])

# mc_filter("1", "z", "11198000")
# mc_filter("2a", "z", "11198400")
# mc_filter("2b", "p", "11198005")
# #
# mc_filter("4a", "z", "11198401")
# mc_filter("4b", "p", "11198410")
# # # mc_filter("4c", "11198410", ["st"])
# mc_filter("4d", "zz", "11198020")

#12:13 1/11/2021 Stop

# mc_filter("5", "12197020", ["p"])
# mc_filter("6", "12197410", ["p"])
#
# mc_filter("7a", "12197400", ["p"])
# mc_filter("7b", "12197022", ["zz"])
# # mc_filter("7c", "12197022", ["st"])
#
# mc_filter("8a", "12197401", ["p"])
# mc_filter("8b_pi", "12197420", ["zz"])
# mc_filter("8b_g", "12197220", ["zz"])
#
# mc_filter("9", "11196000", ["zz"])
# mc_filter("10_pi", "11196410", ["zz"])
# mc_filter("10_g", "11196200", ["zz"])
#
# mc_filter("12_pi", "11196420", ["zz"])
# mc_filter("12_g", "11196210", ["zz"])
# mc_filter("12_pig", "11196620", ["zz"])
#
# mc_filter("13", "13198040", ["s"])
# mc_filter("14", "13198200", ["s"])
# mc_filter("15", "13198400", ["s"])
# mc_filter("16", "13198600", ["s"])

tlist = ["z_spectrum.root",
        "zz_spectrum.root",
        "p_spectrum.root",
        "m_spectrum.root",
        "st_spectrum.root",
        "s_spectrum.root",
        "norm7_data.root",
        "norm8_data.root"
        ]

for i in tlist:
    data_filter(i)
