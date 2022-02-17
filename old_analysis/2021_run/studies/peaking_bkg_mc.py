import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run')
from essentials import *
from uncertainties import correlated_values
from uncertainties.umath import sqrt

zz_mc_list = ["04", "07", "08", "09", "10", "12"]


def convertTuple(tup, spec):
    pname = ""
    for p in tup:
        pname += p
    return pname
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
def get_dwindow_values_fordf(df, tag):
    dwindow_file = ROOT.TFile(f"/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis/2021_run/ws_root_files/d_mass_fits.root","READ")
    dwindow_ws = dwindow_file.Get("d_mass_fits")
    d1_mstart = dwindow_ws.var(f"mean_{tag}_D1").getValV()
    d2_mstart = dwindow_ws.var(f"mean_{tag}_D2").getValV()
    df_mask = df.loc[(abs(df['D1_M'] - d1_mstart) < dwindow) & (abs(df['D2_M'] - d2_mstart) < dwindow)]
    df_next = df_mask.copy()
    dwindow_file.Close()
    return df_next
def submass_harris(infilename, intreename, flag):
    spec = "zz"
    if os.path.exists(outputroot):
        os.remove(outputroot)
    print(f"reading '{intreename}' from '{infilename}'...")
    tl = ["D1H1","D1H2","D2H1","D2H2","KSTH1","KSTH2"]
    columns_to_read_in = ["*H{1,2}_P{X,Y,Z,E}","D1_M","D2_M","B_DTF_M"]
    nlist = [2, 3, 4, 5, 6]
    df_base = rp.read_root(
    infilename,
    intreename,
    columns_to_read_in,)
    df = get_dwindow_values_fordf(df_base, spec)
    for n in nlist:
        Final_n_List = []
        all_n_combinations = itertools.combinations(tl, n)
        for tup in all_n_combinations:
            name = convertTuple(tup, spec)
            print(name)
            tuplist = list(tup)
            df[name] = df.apply(invmass,  args=(tuplist,), axis=1, result_type='expand')
            Final_n_List.append(name)
        print(f"made {n} track for {spec}")
        df[Final_n_List].to_root(f"rootfiles/{flag}_{spec}_track_combinations.root", key=f"DTT", store_index=False, mode='a')
# submass_harris(data_basepath + f"zz_spectrum.root", "DecayTreeTuple", "Data")
# for id in zz_mc_list:
#     zz_file_list = glob.glob(f"{root_basepath}/{id}_zz*.root")
#     submass_harris(zz_file_list, "DecayTreeTuple", f"{id}_MC")
submass_harris(data_basepath+"zz_spectrum.root", "DecayTreeTuple", f"zz_Data")
