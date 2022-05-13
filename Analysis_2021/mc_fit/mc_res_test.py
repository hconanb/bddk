import ROOT as ROOT
from Analysis_2021.essentials import get_free_shapes, save_png, DrawStack
import glob as glob
from array import array
analysis_path = "/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021"
RDF = ROOT.ROOT.RDataFrame

def grab_file_list(id_list):
    file_list = []
    for file_id in id_list:
        file_list = file_list + glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_MC/*/final_sample/{file_id}.root")
    tchain = ROOT.TChain()
    for file_name in file_list:
        tchain.Add(f"{file_name}/DecayTreeTuple_T")
        tchain.Add(f"{file_name}/DecayTreeTuple_nTaT")
    return tchain

def build_tmc_ws(spec_list):

    tmr_list = []
    spec = spec_list[0][0]

    for tuple in spec_list:
        id_list = [tuple[0]]
        mp = tuple[1]

        tc = grab_file_list(id_list)

        spec = id_list[0]

        mc_rdf = RDF(tc)
        mc_rdf = mc_rdf.Define("truth_m_rec", "B_TrueMass - B_DTF_M")

        nbins = 50
        if mp == 1:
            xmin = 5050
            xmax = 5200
        if mp == 2:
            xmin = 4900
            xmax = 5050

        h_bdtfm = mc_rdf.Histo1D((f"{spec}_B_DTF_M", f"{spec}_B_DTF_M", nbins, xmin, xmax), "B_DTF_M")
        h_btruem = mc_rdf.Histo1D((f"{spec}_B_TrueMass", f"{spec}_B_TrueMass", nbins, xmin, xmax), "B_TrueMass")
        h_tmr = mc_rdf.Histo1D((f"{spec}_Btmr", f"{spec}_Btmr", 100, -30, 30), "truth_m_rec")
        tmr_list.append(h_tmr.GetPtr())

        # c1 = ROOT.TCanvas("c1","c1")
        # h_bdtfm_c = h_bdtfm.DrawCopy("")
        # h_btruem_c = h_btruem.DrawCopy("Same")
        #
        # h_bdtfm_c.SetLineColor(ROOT.kRed)
        # h_btruem_c.SetLineColor(ROOT.kBlue)
        #
        # legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
        # #legend.SetHeader(title,"C")
        # legend.AddEntry(f"{spec}_B_DTF_M",f"{spec}_B_DTF_M","l")
        # legend.AddEntry(f"{spec}_B_DTF_M",f"{spec}_B_DTF_M","l")
        #
        # save_png(c1, f"mc_res_test", f"{spec}_TruthM_BDTFM", rpflag = 0)


    for i in tmr_list:
        i.Scale(tmr_list[0].Integral()/i.Integral())

    c2 = ROOT.TCanvas("c2","c2")
    c2.cd()
    hs = DrawStack(ROOT.gPad, tmr_list, legend=(0.20, 0.5, 0.40, 0.9), drawopts="nostack plc pmc")
    save_png(c2, f"mc_res_test", f"{spec}_Tmr", rpflag = 0)

def create_new_tree(spec, data_id, mc_id, mc_file, mp):

        data_fr_path = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/data_fit/fit_data_files/fit_test_noFix_noSplit.root"
        data_fr_file = ROOT.TFile(data_fr_path)
        fws = data_fr_file.Get(f"fit_test_noFix_noSplit")
        fit_pdf = fws.pdf(f"{spec}_{data_id}_fit")
        fit_var = fws.var("B_DTF_M")

        # print(fit_pdf, fit_var)
        real_data_set = fit_pdf.generate(ROOT.RooArgSet(fit_var), int(20000))
        print(real_data_set)
        mc_tc = ROOT.TChain(f"DecayTreeTuple")
        mc_fl = grab_file_list([mc_file])
        for file_name in mc_fl:
            mc_tc.Add(file_name)
        #
        # print(mc_fl)

        if mp == 1:
            xmin = 5050
            xmax = 5200
        if mp == 2:
            xmin = 4900
            xmax = 5050

        mcws = ROOT.RooWorkspace(f"MC")
        mcws.factory(f"B_TrueMass[-1,100000]")
        tmvar = mcws.var("B_TrueMass")

        truth_mc_data_set = ROOT.RooDataSet(f"events", f"events", mc_tc, ROOT.RooArgSet(tmvar))

        outputfile = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_fit/res_test_files/{spec}_{mc_id}.root"
        f = ROOT.TFile(outputfile, "RECREATE")
        tree = ROOT.TTree("tree", "An Example Tree")
        #renamed to B_DTF_M for ease later
        btm = array('f', [0.])
        tree.Branch("B_DTF_M", btm, 'B_DTF_M/F')

        mean_var = fws.var(f"mean_{spec}_{data_id}")
        real_mean = mean_var.getValV()

        for i in range(0, int(truth_mc_data_set.sumEntries())):
            # print (i, truth_mc_data_set.sumEntries())

            temp_mc_argset = truth_mc_data_set.get(i)
            temp_mc_val = temp_mc_argset.getRealValue("B_TrueMass")

            temp_data_argset = real_data_set.get(i)
            temp_data_val = temp_data_argset.getRealValue("B_DTF_M")

            btm[0] = temp_mc_val + temp_data_val - real_mean
            tree.Fill()

        tree.Write("", ROOT.TObject.kOverwrite)
        f.Close()

def compare_res(spec, mc_id, mc_file, mp):

        rec_mc_tc = ROOT.TChain(f"DecayTreeTuple")
        rec_mc_fl = grab_file_list([mc_file])
        for file_name in rec_mc_fl:
            rec_mc_tc.Add(file_name)

        rec_rdf = RDF(rec_mc_tc)
        print(rec_rdf.Count().GetValue())

        new_res_path = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2021/mc_fit/res_test_files/{spec}_{mc_id}.root"
        new_res_file = ROOT.TFile(new_res_path)
        new_res_tree = new_res_file.Get("tree")

        res_rdf = RDF(new_res_tree)
        print(res_rdf.Count().GetValue())

        nbins = 50
        if mp == 1:
            xmin = 5050
            xmax = 5200
        if mp == 2:
            xmin = 4900
            xmax = 5050

        h_rec = rec_rdf.Histo1D((f"{spec}_{mc_id}_rec", f"{spec}_{mc_id}_rec", nbins, xmin, xmax), "B_DTF_M")
        h_res = res_rdf.Histo1D((f"{spec}_{mc_id}_nres", f"{spec}_{mc_id}_nres", nbins, xmin, xmax), "B_DTF_M")

        c2 = ROOT.TCanvas("c2","c2")
        c2.cd()
        hs = DrawStack(ROOT.gPad, [h_rec.GetPtr(), h_res.GetPtr()], legend=(0.20, 0.5, 0.40, 0.9), drawopts="nostack plc pmc")
        save_png(c2, f"mc_res_test", f"{spec}_{mc_id}_comp", rpflag = 0)
