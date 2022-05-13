import sys
import os

basedir = os.getcwd().split('fits')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

RDF = ROOT.ROOT.RDataFrame
opts = ROOT.ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"

from hep_ml import reweight
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import roc_auc_score, roc_curve
from hep_ml.metrics_utils import ks_2samp_weighted
import matplotlib.pyplot as plt
import warnings

fit_base_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/d_window_files/d_mp_mass_fits.root")

fit_ws = fit_base_file.Get(f"d_mp_mass_fits")
fit_var = fit_ws.var("D1_M")
fit_pdf = fit_ws.pdf("D_signal")
data_args = ROOT.RooArgSet(fit_var)


file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/pre_d/01_Z_m_p*.root"
DecayTree_List = ["DecayTreeTuple"]
tchain = grab_files_and_chain(file_path, DecayTree_List)
data_set = ROOT.RooDataSet(f"events", f"events", tchain, data_args)

c1 = ROOT.TCanvas("c1","c1")

frame = fit_var.frame()

data_set.plotOn(frame, ROOT.RooFit.Name("events"))
fit_pdf.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name(f"fit"))

frame.Draw()

save_pdf(c1, f"fit_d_check", f"z", rpflag = 0)
