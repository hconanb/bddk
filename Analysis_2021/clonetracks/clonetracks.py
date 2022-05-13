#!/usr/bin/python
import sys
import os
import ROOT
import pandas as pd
import uproot
import glob
import datetime
import itertools
RDF = ROOT.ROOT.RDataFrame


def twotracktheta(p):
    """arguments:
    plist -- list of particles, e.g., ["K1", "K2", ...]
    """
    x = f"({p[0]}_PX*{p[1]}_PX)"
    y = f"({p[0]}_PY*{p[1]}_PY)"
    z = f"({p[0]}_PZ*{p[1]}_PZ)"
    n = f"({x}+{y}+{z})"
    x0 = f"pow({p[0]}_PX,2)"
    y0 = f"pow({p[0]}_PY,2)"
    z0 = f"pow({p[0]}_PZ,2)"
    s0 = f"sqrt({x0}+{y0}+{z0})"
    x1 = f"pow({p[1]}_PX,2)"
    y1 = f"pow({p[1]}_PY,2)"
    z1 = f"pow({p[1]}_PZ,2)"
    s1 = f"sqrt({x1}+{y1}+{z1})"
    f = f"acos({n}/({s0}*{s1}))"
    print(f)
    return f
def convertTuple(tup, spec):
    pname = ""
    for p in tup:
        pname += p
    return pname
def plot_clone_tracks(spec, tchain):
    rdf = RDF(tchain)
    if "norm" not in spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "KSTH1", "KSTH2"]
    if "norm" in spec:
        tl = ["D1H1", "D1H2", "D2H1", "D2H2", "D2H3", "D2H4","K"]
    if "Z_m_p" in spec:
        tl = tl + ["D1H3", "D2H3"]
    if "P_z_p" in spec:
        tl = tl + ["D2H3"]
        # Final_n_List = []
    all_n_combinations = itertools.combinations(tl, 2)
    tupflag = 0
    arr = ""
    for tup in all_n_combinations:
        tname = convertTuple(tup, spec)
        arr = f"cc_{tname}, {arr}"
        rdf = rdf.Define(f"cc_{tname}", twotracktheta(tup))
    arr = "{" + arr + "}"
    print(arr)
    rdf = rdf.Define("minpairvector",f"std::array<double, 28> v{arr}; return v;")
    # clist = rdf.GetColumnNames()
    # rdfsnap = rdf.Snapshot(f"DecayTreeTuple", f"{spec}.root", clist)

for spec in ["Z_m_p"]:
    inputfile = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_data/*/post_d/{spec}.root")
    tchain = ROOT.TChain()
    for file_name in inputfile:
        tchain.Add(f"{file_name}/DecayTreeTuple")
    plot_clone_tracks(spec, tchain)
