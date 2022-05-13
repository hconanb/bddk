import os

gc_onflag = 0
fit_strat = "dp_f"

#spec is spectrum, tuples are (mc_id, fit shape)
specs = ["Z_m_p","Z_z_z","P_z_p","M_m_z","P_z_pst","Zs_sm_p"]
#,"zz","p","m","st"]

run_name = f"zo_{fit_strat}_{gc_onflag}"
# plot_data(run_name, specs)

sl = " ".join(str(x) for x in specs)

os.system(f"python3 data_plot_script.py\
            --run_name {run_name}\
            --fit_strat {fit_strat}\
            --specs {sl}\
            --data_only_flag 1\
            ")
def plot_post_d_data(spec):

    file_list = glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_data/*/post_d/{spec}_postdcuts.root")
    tchain = ROOT.TChain("DecayTreeTuple")
    for file_name in file_list:
        tchain.Add(file_name)

        if spec == "Z_m_p":
            title = "D^{-} D^{+} K^{*0}"
        if spec == "Z_z_z":
            title = "#bar{D^{0}} D^{0} K^{*0}"
        if spec == "P_z_p":
            title = "#bar{D^{0}} D^{+} K^{*0}"
        if spec == "M_m_z" :
            title = "D^{-} D^{0} K^{*0}"
        if spec == "P_z_pst" :
            title = "#bar{D^{0}} (D^{*+} #rightarrow D^{0} #pi+) K^{*0}"
        if spec == "Zs_sm_p" :
            title = "D_{s}^{-} D^{+} K^{*0}"

    rdf = RDF(tchain)
    xmax = 5600
    xmin = 4800
    bins = 100
    hist_bdtfm = rdf.Histo1D((f"bdtfm_{spec}", f"bdtfm_{spec}", bins, xmin, xmax), 'B_DTF_M')
    xr = xmax - xmin
    c1 = ROOT.TCanvas("c1","c1")
    hist_bdtfm.GetXaxis().SetTitle(f"m({title}) with DTF Constraints [MeV]")
    hist_bdtfm.GetYaxis().SetTitle(f"Events / ({xr/bins})")
    hist_bdtfm.Draw("E")

    save_png(c1, "b_post_d", f"{spec}_b_post_d", rpflag = 0)
