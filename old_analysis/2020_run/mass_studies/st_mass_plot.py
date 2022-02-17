import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
import uncertainties
from createXFD import *

RDF = ROOT.ROOT.RDataFrame

CreateXFD(data_basepath+"st_spectrum.root", "DecayTreeTuple", data_basepath + "st_spectrum_newfd.root", "DecayTreeTuple")
CreateXFD(TT_data_basepath+"TT_st_spectrum.root", "DecayTreeTuple_nTaT", TT_data_basepath + "st_spectrum_newfd_nTaT.root", "DecayTreeTuple")
CreateXFD(TT_data_basepath+"TT_st_spectrum.root", "DecayTreeTuple_T", TT_data_basepath + "st_spectrum_newfd_T.root", "DecayTreeTuple")

# dx2 = "pow(B_ENDVERTEX_X - D2_ENDVERTEX_X,2)"
# dy2 = "pow(B_ENDVERTEX_Y - D2_ENDVERTEX_Y,2)"
# dz2 = "pow(B_ENDVERTEX_Z - D2_ENDVERTEX_Z,2)"
#
# ux2 = "pow(B_ENDVERTEX_XERR,2) + pow(D2_ENDVERTEX_XERR,2)"
# uy2 = "pow(B_ENDVERTEX_YERR,2) + pow(D2_ENDVERTEX_YERR,2)"
# uz2 = "pow(B_ENDVERTEX_ZERR,2) + pow(D2_ENDVERTEX_ZERR,2)"

# st_data_file = ROOT.TFile(data_basepath+"st_spectrum.root")
# st_tree = st_data_file.Get("DecayTreeTuple")
#
#
# outHistFile = ROOT.TFile("test.root","RECREATE")
# outHistFile.cd()
# st_test_tree = st_tree.CopyTree(st_m_cut)
# st_test_tree.Write()
# outHistFile.Close()
# fd_new = f"sqrt({dx2} + {dy2} + {dz2})"
# fd_new_u = f"sqrt({dx2}*{ux2} + {dy2}*{uy2} + {dz2}*{uz2})/{fd_new}"
#
#
# st_m_cut = f"(abs(D1_M - {d0mass}) < {dwindow} && abs(D2_M - {d0mass}) < {dwindow})"

# rdf_base = RDF(st_tree)
# rdf = rdf_base.Define("DstmD", "D2st_M - D2_M")\
#               .Define("fd_new", fd_new)\
#               .Define("fd_newx2", fd_new_u)\
#               .Filter(st_m_cut)
#               .Foreach([](int i){ std::cout << i << std::endl;}, {"myIntColumn"});
#
#
#
# dbins = 100
# dmin = 1820
# dmax = 1920
# dstmin = 1900
# dstmax = 2200
#
# # D1_M = rdf.Histo1D((f"D1_M", f"D1_M", dbins, dmin, dmax), "D1_M")
# # D2_M = rdf.Histo1D((f"D2_M", f"D2_M", dbins, dmin, dmax), "D2_M")
# # D2_st_M = rdf.Histo1D((f"D2st_M", f"D2st_M", dbins, dstmin, dstmax), "D2st_M")
# # D2mm = rdf.Histo2D(("", "", dbins, dmin, dmax, dbins, 120, 220), "D2_M", "DstmD")
# # D2mmm = rdf.Histo1D(("", "", dbins, 130, 220), "DstmD")
# dstfdx2 = rdf.Histo1D(("", "", dbins, -3, 20), "fd_newx2")
#
# # c1 = ROOT.TCanvas("c1","c1")
# # D2m = D2mmm.DrawCopy()
# # D2m.GetXaxis().SetTitle("D^{*+} - D^{0} Mass [MeV]")
# # c1.SaveAs(f"st_M.png")
#
# c2 = ROOT.TCanvas("c2","c2")
# dstfdx2_c = dstfdx2.DrawCopy()
# dstfdx2_c.GetXaxis().SetTitle("D^{*+} fdx2")
# c2.SaveAs(f"x2st_M.png")
