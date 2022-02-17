import sys

sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *

def grab_file_list(type, event_list):
    file_list = []
    for event in event_list:
        file_list = file_list + glob.glob(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type}_2021/*/post_d/*{event}*.root")
    return file_list

def build_mc_ws(ws_name, event_list, shape_list, mean_guess, mean_window, fit_window, wflag = 0):

    print(f"mcws is {ws_name}")
    mcws = ROOT.RooWorkspace(ws_name)
    b_min = mean_guess - fit_window
    b_max = mean_guess + fit_window

    mcws.factory(f"B_DTF_M[{b_min}, {b_max}]")
    for shape_flag in shape_list:
        s_list = [ws_name, shape_flag, mean_guess, mean_window]
        get_free_shapes(mcws, ws_name, "MC", s_list)
    mcws.Print()
    b_dtf_m = mcws.var("B_DTF_M")

    if wflag == 0:
        tc = ROOT.TChain(f"DecayTreeTuple")
        fl = grab_file_list("MC", event_list)
        for file_name in fl:
            tc.Add(file_name)
            print(file_name)
        data_args = ROOT.RooArgSet(b_dtf_m)
        data_set = ROOT.RooDataSet(f"{ws_name}_events", f"{ws_name}_events", tc, data_args)
        mcws.Import(data_set)
        mcws.writeToFile(f"{analysis_path}/mc_fit/base_mc_files/{ws_name}.root")
        print(f"Saved: {analysis_path}/mc_fit/base_mc_files/{ws_name}.root")
        mcws.Print()
    #
    # # if wflag > 0:
    # #
    # #     data_set_list = []
    # #
    # #     tc_1 = ROOT.TChain(f"DecayTreeTuple")
    # #     tc_2 = ROOT.TChain(f"DecayTreeTuple")
    # #
    # #     fl_1 = grab_file_list("MC", [event_list[0]])
    # #     for file_name in fl_1:
    # #         tc_1.Add(file_name)
    # #
    # #     fl_2 = grab_file_list("MC", [event_list[1]])
    # #     for file_name in fl_2:
    # #         tc_2.Add(file_name)
    # #
    # #     if len(event_list) == 3:
    # #         fl_3 = grab_file_list("MC", [event_list[2]])
    # #         for file_name in fl_3:
    # #             tc_3.Add(file_name)
    # #
    # #
    # #     vallist = np.linspace(0, 1, wflag, endpoint = False)
    # #
    # #     name_valist = list(range(0, wflag))
    # #
    # #     if len(event_list) == 2:
    # #         for wval_1, nid in zip(vallist, name_valist):
    # #
    # #             wval_2 = 1 - wval_1
    # #
    # #             data_args_1 = ROOT.RooArgSet(b_dtf_m)
    # #             data_args_2 = ROOT.RooArgSet(b_dtf_m)
    # #
    # #             data_set_1 = ROOT.RooDataSet(f"{name_mc_ws}_events_{nid}", f"{name_mc_ws}_events_{nid}", tc_1, data_args_1)
    # #             data_set_2 = ROOT.RooDataSet(f"{name_mc_ws}_events_2_{nid}_temp", f"{name_mc_ws}_events_2_{nid}_temp", tc_2, data_args_2)
    # #
    # #             wvar_1 = ROOT.RooRealVar("wvar_1","wvar_1", wval_1)
    # #             wvar_2 = ROOT.RooRealVar("wvar_2","wvar_2", wval_2)
    # #
    # #             data_set_1.addColumn(wvar_1)
    # #             data_set_2.addColumn(wvar_2)
    # #
    # #             data_args_1.add(wvar_1)
    # #             data_args_2.add(wvar_2)
    # #
    # #             new_data_set_1 = ROOT.RooDataSet(data_set_1.GetName(), data_set_1.GetTitle(), data_set_1, data_args_1, "B_DTF_M > 0", "wvar_1");
    # #             new_data_set_2 = ROOT.RooDataSet(data_set_2.GetName(), data_set_2.GetTitle(), data_set_2, data_args_2, "B_DTF_M > 0", "wvar_2");
    # #
    # #             new_data_set_1.append(new_data_set_2)
    # #             data_set_list.append(new_data_set_1)
    #
    #
    #
    # else:
    #     for ds in data_set_list:
    #         mcws.Import(ds)
    #     mcws.writeToFile(f"../mc_ws_root_files/{name_mc_ws}.root")
    #     mcws.Print()

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--ws_name')
parser.add_argument('--event_list', nargs='+')
parser.add_argument('--shape_list', nargs='+')
parser.add_argument('--mean_guess', type = float)
parser.add_argument('--mean_window', type = float)
parser.add_argument('--fit_window', type = float)

args = parser.parse_args()
ws_name = args.ws_name
event_list = args.event_list
shape_list = args.shape_list
mean_guess  = args.mean_guess
mean_window  = args.mean_window
fit_window  = args.fit_window

build_mc_ws(ws_name, event_list, shape_list, mean_guess, mean_window, fit_window)
