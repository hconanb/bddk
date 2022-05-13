import sys
import os

basedir = os.getcwd().split('fits')[0]
sys.path.append(basedir)

from rootutils import residualPlot
from essential_functions import *

def fit_mc(spec, tchain, varoi, shape_list):

    mc_ws = ROOT.RooWorkspace(spec)

    mean_start =  id_to_b_values_dict[spec][0]
    b_window =  id_to_b_values_dict[spec][1]

    b_min = mean_start - b_window
    b_max = mean_start + b_window

    mc_ws.factory(f"{varoi}[{b_min}, {b_max}]")

    for shape_flag in shape_list:
        get_roofit_pdf(mc_ws, spec, shape_flag, varoi)

    mc_ws.Print()

    fit_var = mc_ws.var(varoi)
    fit_args = ROOT.RooArgSet(fit_var)
    data_set = ROOT.RooDataSet(f"{spec}_events", f"{spec}_events", tchain, fit_args)

    fit_ws = ROOT.RooWorkspace(f"fit_ws")

    for shape_flag in shape_list:

        fit_pdf = mc_ws.pdf(f"{spec}_fit_{shape_flag}")
        fit_result = fit_pdf.fitTo(data_set, ROOT.RooFit.PrintLevel(0), ROOT.RooFit.Save())

        fit_ws.Import(fit_pdf)
        fit_ws.Import(fit_result)

    fit_ws.Import(data_set)
    fit_ws.Print()

    fit_output_file = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/mc_files/{id_to_spec_dict[spec]}.root"
    fit_ws.writeToFile(fit_output_file)

def plot_mc(spec, varoi, shape_list, binning):

    for shape_flag in shape_list:

        fit_input_file = f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/Analysis_2022/fits/mc_files/{id_to_spec_dict[spec]}.root"
        plot_base_file = ROOT.TFile(fit_input_file)
        fit_ws = plot_base_file.Get(f"fit_ws")

        fit_var = fit_ws.var(varoi)

        frame = fit_var.frame(ROOT.RooFit.Title(spec), ROOT.RooFit.Bins(binning))

        data_set = fit_ws.data(f"{spec}_events")

        data_set.plotOn(frame, ROOT.RooFit.Name("events"))

        fit_pdf = fit_ws.pdf(f"{spec}_fit_{shape_flag}")
        fit_pdf.plotOn(frame, ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kViolet), ROOT.RooFit.Name("total_fit") )

        if shape_flag == "DG" or "_fr" in shape_flag:
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit_{shape_flag}_a"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Name("fa"))
            fit_pdf.plotOn(frame, ROOT.RooFit.Components(f"{spec}_fit_{shape_flag}_b"), ROOT.RooFit.LineStyle(ROOT.kDashed), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Name("fb"))

        d_chi2 = frame.chiSquare("total_fit", "events")

        dtpave = ROOT.TPaveText(0.20, 0.65, 0.40, 0.85, "NB NDC")
        dtpave.SetFillStyle(0)
        dtpave.AddText(f"#chi^{{2}}: {round(d_chi2, 3)}")
        dtpave.Draw()
        frame.addObject(dtpave)

        p = residualPlot()
        p.pt.cd()
        xaxis = frame.GetXaxis()
        xaxis.SetTickLength(0)
        xaxis.SetNdivisions(0)
        xaxis.SetLabelSize(0)
        frame.Draw()

        legend = ROOT.TLegend(0.70, 0.4, 0.90, 0.93)
        shape_on_l = shape_flag.split("_fr")[0]
        legend.AddEntry("total_fit",f"Fit Choice: {shape_on_l}","l")
        legend.AddEntry("total_fit",f"MC Scheme: {id_to_scheme_dict[id_to_spec_dict[spec]]}","l")
        if shape_flag == "DG":
            legend.AddEntry("fa","Gauss 1","l")
            legend.AddEntry("fb","Gauss 2","l")
        if "_fr" in shape_flag:
            shape_1 = shape_flag.split("Add")[0]
            shape_2 = shape_flag.split("Add")[1].split("_fr")[0]
            legend.AddEntry("fa",f"{shape_1}","l")
            legend.AddEntry("fb",f"{shape_2}","l")
        legend.AddEntry("data","Simulated Events","lep")
        legend.SetTextSize(0.050)
        legend.SetFillStyle(1001)
        legend.Draw()

        p.pb.cd()
        hpull = frame.pullHist("events","total_fit")
        pull = fit_ws.var(f"{varoi}").frame();
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

        title = get_decay_string(id_to_meson_string_dict[spec], "Rec")
        pull.GetXaxis().SetTitle(f"{title} [MeV]")

        save_pdf(p, f"mc_fits", f"{spec}_{shape_flag}", rpflag = 1)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Fit a PDF to final MC samples")
    parser.add_argument('--Spec_List', choices=id_to_spec_dict.keys(),  nargs="+", help = 'Spec_List')
    parser.add_argument('--Varoi', choices=["B_DTF_M","B_M"], default = "B_DTF_M", help = 'Variable to Fit')
    parser.add_argument("--Shape_List", nargs="+")
    parser.add_argument("--Binning", type = float, default = 50)
    parser.add_argument('--FIT', action='store_true')
    parser.add_argument('--PLOT', action='store_true')
    args = parser.parse_args()

    spec_list = args.Spec_List
    varoi = args.Varoi
    shape_list = args.Shape_List
    binning = args.Binning
    fit_flag = args.FIT
    plot_flag = args.PLOT

    if spec_list is None:
        spec_list = ids_to_bestfit_dict.keys()

    for spec in spec_list:

        file_path = f"/mnt/c/Users/Harris/Desktop/rootfiles/2022_filtered_mc/*/final_sample/{id_to_spec_dict[spec]}.root"
        DecayTree_List = ["DecayTreeTuple_SIG"]
        tchain = grab_files_and_chain(file_path, DecayTree_List)

        if shape_list is None:
            shape = ids_to_bestfit_dict[spec]
            if fit_flag:
                    fit_mc(spec, tchain, varoi, [shape])
            if plot_flag:
                    plot_mc(spec, varoi, [shape], binning)
        else:
            if fit_flag:
                    fit_mc(spec, tchain, varoi, shape_list)
            if plot_flag:
                    plot_mc(spec, varoi, shape_list, binning)
