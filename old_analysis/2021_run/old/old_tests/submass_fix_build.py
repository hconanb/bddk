import sys

sys.path.append("/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run")
from essentials import *

fixlist = ["B_M", "B_VtxChi2_", "B_VtxnDoF", "B_VtxM", "B_IP", "B_IPChi2_"]


def get_dwindow_values_fordf(df, tag):
    dwindow_file = ROOT.TFile(
        f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/ws_root_files/d_mass_fits.root",
        "READ",
    )
    dwindow_ws = dwindow_file.Get("d_mass_fits")
    d1_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
    d2_mstart = dwindow_ws.var(f"mean_{tag}_D2").getValV()
    df_mask = df.loc[
        (abs(df["D1_M"] - d1_mstart) < dwindow)
        & (abs(df["D2_M"] - d2_mstart) < dwindow)
    ]
    df_next = df_mask.copy()
    dwindow_file.Close()
    return df_next


def get_fix_id(row, D1KSTH1, M02, M03, M12, M13):
    # FIX_0_ID = ""
    # FIX_1_ID = ""
    # FIX_2_ID = ""
    # FIX_3_ID = ""
    FIX_ID = "0"
    if row[D1KSTH1] == row[M02]:
        # FIX_0_ID = "D0bar"
        # FIX_1_ID = "D0"
        # FIX_2_ID = "K"
        # FIX_3_ID = "pi"
        FIX_ID = "0123"
    elif row[D1KSTH1] == row[M03]:
        # FIX_0_ID = "D0bar"
        # FIX_1_ID = "D0"
        # FIX_2_ID = "pi"
        # FIX_3_ID = "K"
        FIX_ID = "0132"
    elif row[D1KSTH1] == row[M12]:
        # FIX_0_ID = "D0"
        # FIX_1_ID = "D0bar"
        # FIX_2_ID = "K"
        # FIX_3_ID = "pi"
        FIX_ID = "1023"
    elif row[D1KSTH1] == row[M13]:
        # FIX_0_ID = "D0"
        # FIX_1_ID = "D0bar"
        # FIX_2_ID = "pi"
        # FIX_3_ID = "K"
        FIX_ID = "1032"
    for var in fixlist:
        row[f"{var}01_Fixed"] = row[f"{var}01"]
        row[f"{var}02_Fixed"] = row[f"{var}{FIX_ID[0]}{FIX_ID[2]}"]
        row[f"{var}03_Fixed"] = row[f"{var}{FIX_ID[0]}{FIX_ID[3]}"]
        row[f"{var}12_Fixed"] = row[f"{var}{FIX_ID[1]}{FIX_ID[2]}"]
        row[f"{var}13_Fixed"] = row[f"{var}{FIX_ID[1]}{FIX_ID[3]}"]
        row[f"{var}23_Fixed"] = row[f"{var}23"]
        row[f"{var}012_Fixed"] = row[f"{var}01{FIX_ID[2]}"]
        row[f"{var}013_Fixed"] = row[f"{var}01{FIX_ID[3]}"]
        row[f"{var}023_Fixed"] = row[f"{var}{FIX_ID[0]}23"]
        row[f"{var}123_Fixed"] = row[f"{var}{FIX_ID[1]}23"]
    return row


def invmass(row, plist):
    """arguments:
    row -- pandas DataFrame row
    plist -- list of particles, e.g., ["K1", "K2", ...]
    """

    def getSum(x):
        return pow(sum(row[f"{p}_{x}"] for p in plist), 2)

    out = getSum("PE")
    for x in ("PX", "PY", "PZ"):
        out -= getSum(x)
    return np.sqrt(out)


def submass_harris(infilename, intreename, flag):
    spec = "zz"
    print(f"reading '{intreename}' from '{infilename}'...")
    tl = ["D1", "D2", "KSTH1", "KSTH2"]
    tl_sb = ["0", "1", "2", "3"]
    columns_to_read_in = [
        "D{1,2}_P{X,Y,Z,E}",
        "KSTH*{1,2}_P{X,Y,Z,E}",
        "D1_M",
        "D2_M",
        "B_*",
    ]
    nlist = [2]
    nsblist = [2, 3]
    df_base = rp.read_root(infilename, intreename, columns_to_read_in)
    df = get_dwindow_values_fordf(df_base, spec)
    Final_n_List = []
    for n in nlist:
        all_n_combinations = itertools.combinations(tl, n)
        for tup in all_n_combinations:
            name = ""
            for p in tup:
                name += p
            tuplist = list(tup)
            df[name] = df.apply(invmass, args=(tuplist,), axis=1, result_type="expand")
    for n in nsblist:
        all_n_sb_combinations = itertools.combinations(tl_sb, n)
        for tup in all_n_sb_combinations:
            name = ""
            for p in tup:
                name += p
            for var in fixlist:
                new_name = f"{var}{name}_Fixed"
                Final_n_List.append(new_name)
    print(f"Fixing Submasses for {infilename}")
    df = df.apply(
        get_fix_id,
        args=("D1KSTH1", "B_M02", "B_M03", "B_M12", "B_M13"),
        axis=1,
        result_type="expand",
    )
    print(f"Saving {flag}_allyear_fixed.root")
    df[Final_n_List].to_root(
        f"rootfiles/{flag}_allyear_fixed.root", key=f"{intreename}", store_index=False
    )
    print(f"Saving {flag}_allyear_fixed.root")


zz_mc_list = ["04", "07", "08", "09", "10", "12"]
for id in zz_mc_list:
    zz_file_list = glob.glob(f"{root_basepath}/{id}_zz*.root")
    submass_harris(zz_file_list, "DecayTreeTuple", f"{id}_MC")
