def CreateXFD(infilename, intreename, outfilename, outtreename):
    """write TTree with X_myFD and X_myFDERR branches

    arguments:
    infilename -- ROOT file to read
    intreename -- name of TTree to read from infilename
    outfilename -- ROOT file to write
    outtreename -- name of TTree to write to outfilename
    """

    from uncertainties import correlated_values
    from uncertainties.umath import sqrt
    import root_pandas as rp


    # -- define column transformations
    def corrval(row, xbase, part):
        "correlate values using covariance matrix"

        return correlated_values(
            [
                getattr(row, f"{part}_{xbase}_X"),
                getattr(row, f"{part}_{xbase}_Y"),
                getattr(row, f"{part}_{xbase}_Z"),
            ],
            [
                getattr(row, f"{part}_{xbase}_COV_")[0:3],
                getattr(row, f"{part}_{xbase}_COV_")[3:6],
                getattr(row, f"{part}_{xbase}_COV_")[6:9]
            ],
        )
    def FD(row):
        "calculate distance from myPV to myENDV"
        return sqrt(
            (row.B_myENDV_X - row.D2_myENDV_X) ** 2
            + (row.B_myENDV_Y - row.D2_myENDV_Y) ** 2
            + (row.B_myENDV_Z - row.D2_myENDV_Z) ** 2
        )
    def ns(row, colname):
        "return tuple of ufloat.n, .s"
        v = getattr(row, colname)
        fdx = (v.n/v.s) ** 2
        return (v.n, v.s, fdx)

    columns_to_read_in= ["B_DTF_M", "D1_M", "D2_M", "D2st_M", "D1_DIRA_ORIVX", "D2_DIRA_ORIVX", "D1_FDCHI2_ORIVX"] + [f"{part}_ENDVERTEX_{coord}" for coord in ("X", "Y", "Z", "COV_") for part in ("B", "D2")]
    columns_to_read_out = ["B_DTF_M", "D1_M", "D2_M", "D2st_M", "D1_DIRA_ORIVX", "D2_DIRA_ORIVX", "D1_FDCHI2_ORIVX"] + ["D2_FDCHI2_ORIVX"]
    # -- read input
    print(f"reading '{intreename}' from '{infilename}'...")
    df = rp.read_root(
        infilename,
        intreename,
        columns_to_read_in,
    )
    # df = df[1:]

    # -- transform columns
    print("applying transformations...")
    df[["B_myENDV_X", "B_myENDV_Y", "B_myENDV_Z"]] = df.apply(
        corrval, args=["ENDVERTEX", "B"], axis=1, result_type="expand"
    )
    df[["D2_myENDV_X", "D2_myENDV_Y", "D2_myENDV_Z"]] = df.apply(
        corrval, args=["ENDVERTEX", "D2"], axis=1, result_type="expand"
    )
    df["X_myFD_u"] = df.apply(FD, axis=1)
    df[["X_myFD", "X_myFDERR", "D2_FDCHI2_ORIVX"]] = df.apply(
        ns, args=["X_myFD_u"], axis=1, result_type="expand"
    )

    # -- write output
    print(f"writing to '{outfilename}'...")
    df[columns_to_read_out].to_root(outfilename, key=outtreename, store_index=False)
    print("done")
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=CreateXFD.__doc__,
    )
    parser.add_argument("infilename")
    parser.add_argument("intreename")
    parser.add_argument("outfilename")
    args = parser.parse_args()

    CreateXFD(args.infilename, args.intreename, args.outfilename, args.intreename)
