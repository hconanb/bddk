import sys
sys.path.append('/mnt/c/Users/Harris/Google Drive/LHCb/bddk/analysis')
from essentials import *
from uncertainties import correlated_values
from uncertainties.umath import sqrt
import root_pandas as rp
import itertools

z_branch_dict = {
"D1H1" : "Kp_Dm",
"D1H2" : "pim1_Dm",
"D2H1" : "Km_Dp",
"D2H2" : "pip1_Dp",
"KSTH1" : "Kp_Kst0",
"KSTH2" : "pim_Kst0",
}
zz_branch_dict = {
"D1H1" : "Kp_D0bar",
"D1H2" : "pim_D0bar",
"D2H1" : "Km_D",
"D2H2" : "pip1_Dp",
"D2H3" : "pip2_Dp",
"KSTH1" : "Kp_Kst0",
"KSTH2" : "pim_Kst0",
}
p_branch_dict = {
"D1H1" : "Kp_D0bar",
"D1H2" : "pim_D0bar",
"D2H1" : "Km_Dp",
"D2H2" : "pip1_Dp",
"D2H3" : "pip2_Dp",
"KSTH1" : "Kp_Kst0",
"KSTH2" : "pim_Kst0",
}
p_branch_dict = {
"D1H1" : "Kp_Dm",
"D1H2" : "pim1_Dm",
"D1H3" : "pim2_Dm",
"D2H2" : "pip1_Dp",
"D2H3" : "pip2_Dp",
"KSTH1" : "Kp_Kst0",
"KSTH2" : "pim_Kst0",
}
def convertTuple(tup, spec):

    # if sspectrum == "zz":
    #     d1 = "#bar}D^}0}}"
    #     d2 = "D^{0}"
    #     d1h1 = "K^{+}"
    #     d1h2 = "#pi^{-}"
    #     d2h1 = "K^{-}"
    #     d2h2 = "#pi^{+}"
    pname = ""
    # for p in tup:
    #     pname += z_branch_dict[p]+"_"
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
speclist = ["p","m"] # ["m", "st"]
#["z",["zz", "p",
intreename = "DecayTreeTuple"

for spec in speclist:
    outputroot = f"{spec}_track_combinations.root"
    if os.path.exists(outputroot):
        os.remove(outputroot)
    infilename = data_basepath +f"{spec}_spectrum.root"
    if spec == "z":
        tl = ["D1H1","D1H2","D1H3","D2H1","D2H2","D2H3","KSTH1","KSTH2"]
        columns_to_read_in = ["*H{1,2,3}_P{X,Y,Z,E}","D1_M","D2_M","B_DTF_M"]
        nlist = [2, 3, 4, 5, 6, 7, 8]
    if spec == "p" or spec == "st":
        tl = ["D1H1","D1H2","D2H1","D2H2","D2H3","KSTH1","KSTH2"]
        columns_to_read_in = ["*H{1,2,3}_P{X,Y,Z,E}","D1_M","D2_M","B_DTF_M"]
        nlist = [2, 3, 4, 5, 6, 7]
    if spec == "m":
        tl = ["D1H1","D1H2","D1H3","D2H1","D2H2","KSTH1","KSTH2"]
        columns_to_read_in = ["*H{1,2,3}_P{X,Y,Z,E}","D1_M","D2_M","B_DTF_M"]
        nlist = [2, 3, 4, 5, 6, 7]
    if spec == "zz":
        tl = ["D1H1","D1H2","D2H1","D2H2","KSTH1","KSTH2"]
        columns_to_read_in = ["*H{1,2}_P{X,Y,Z,E}","D1_M","D2_M","B_DTF_M"]
        nlist = [2, 3, 4, 5, 6]
    print(f"reading '{intreename}' from '{infilename}'...")
    df_base = rp.read_root(
    infilename,
    intreename,
    columns_to_read_in,
)
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
        df[Final_n_List].to_root(f"{spec}_track_combinations.root", key=f"DTT", store_index=False, mode='a')
