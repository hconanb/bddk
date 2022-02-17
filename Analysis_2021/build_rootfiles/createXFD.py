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
    print(infilename)
    # -- define column transformations
    def corrval(row, xbase, part):

        # if mc_flag == 0:
        #     return correlated_values(
        #         [
        #             float(getattr(row, f"{part}_{xbase}_X")),
        #             float(getattr(row, f"{part}_{xbase}_Y")),
        #             float(getattr(row, f"{part}_{xbase}_Z")),
        #         ],
        #         [
        #             getattr(row, f"{part}_{xbase}_COV_")[0:3],
        #             getattr(row, f"{part}_{xbase}_COV_")[3:6],
        #             getattr(row, f"{part}_{xbase}_COV_")[6:9],
        #         ],
        #     )
        # if mc_flag == 1:
            # print (getattr(row, f"{part}_{xbase}_COV_")[0], getattr(row, f"{part}_{xbase}_COV_")[1], getattr(row, f"{part}_{xbase}_COV_")[2])
        return correlated_values(
            [
                float(getattr(row, f"{part}_{xbase}_X")),
                float(getattr(row, f"{part}_{xbase}_Y")),
                float(getattr(row, f"{part}_{xbase}_Z")),
            ],
            [
                getattr(row, f"{part}_{xbase}_COV_")[0],
                getattr(row, f"{part}_{xbase}_COV_")[1],
                getattr(row, f"{part}_{xbase}_COV_")[2],
            ],
        )

    def FD_d1(row):
        "calculate distance from myPV to myENDV"
        return sqrt(
                (row.B_myENDV_X - row.D1_myENDV_X) ** 2
            + (row.B_myENDV_Y - row.D1_myENDV_Y) ** 2
            + (row.B_myENDV_Z - row.D1_myENDV_Z) ** 2
        )

    def FD_d2(row):
        "calculate distance from myPV to myENDV"
        return sqrt(
                (row.B_myENDV_X - row.D2_myENDV_X) ** 2
            + (row.B_myENDV_Y - row.D2_myENDV_Y) ** 2
            + (row.B_myENDV_Z - row.D2_myENDV_Z) ** 2
        )

    def pd1(row):
        "calculate distance from myPV to myENDV"
        return sqrt((row.D1_PX) ** 2 + (row.D1_PY) ** 2 + (row.D1_PZ) ** 2)

    def pd2(row):
        "calculate distance from myPV to myENDV"
        return sqrt((row.D2_PX) ** 2 + (row.D2_PY) ** 2 + (row.D2_PZ) ** 2)

    def fd1_dot_p(row):
        dot = (
            (row.D1_PX * (row.D1_myENDV_X - row.B_myENDV_X))
            + (row.D1_PY * (row.D1_myENDV_Y - row.B_myENDV_Y))
            + (row.D1_PZ * (row.D1_myENDV_Z - row.B_myENDV_Z))
        )
        return dot

    def fd2_dot_p(row):
        dot = (
            (row.D2_PX * (row.D2_myENDV_X - row.B_myENDV_X))
            + (row.D2_PY * (row.D2_myENDV_Y - row.B_myENDV_Y))
            + (row.D2_PZ * (row.D2_myENDV_Z - row.B_myENDV_Z))
        )
        return dot

    def ns(row, colname):
        "return tuple of ufloat.n, .s"
        v = getattr(row, colname)
        fdx = (v.n / v.s) ** 2
        return (v.n, v.s, fdx)

    def ns2(row, x, p, dot):
        xb = getattr(row, x)
        pb = getattr(row, p)
        dotb = getattr(row, dot)
        return dotb.n / (xb.n * pb)

    df = rp.read_root(
        infilename,
        intreename,
        ["*"]
    )
    # df = df[1:]
    # -- transform columns
    # print("applying transformations...")
    l = df.columns

    df[["B_myENDV_X", "B_myENDV_Y", "B_myENDV_Z"]] = df.apply(
        corrval, args=["ENDVERTEX", "B"], axis=1, result_type="expand"
    )

    if "mst" in intreename:
        df[["D1_myENDV_X", "D1_myENDV_Y", "D1_myENDV_Z"]] = df.apply(
            corrval, args=["ENDVERTEX", "D1"], axis=1, result_type="expand"
        )
        df["X_myFD_u"] = df.apply(FD_d1, axis=1)
        df[["p2"]] = df.apply(pd1, axis=1)
        df[["B_D1_DOT"]] = df.apply(
            fd1_dot_p,
            axis=1,
        )
        df[["D1_DIRA_ORIVX_FIXED"]] = df.apply(
            ns2, axis=1, args=["X_myFD_u", "p2", "B_D1_DOT"], result_type="expand"
        )
        df["D1stmD_M"] = df["D1st_M"] - df["D1_M"]


    if "pst" in outfilename:
        df[["D2_myENDV_X", "D2_myENDV_Y", "D2_myENDV_Z"]] = df.apply(
            corrval, args=["ENDVERTEX", "D2"], axis=1, result_type="expand"
        )
        df["X_myFD_u"] = df.apply(FD_d2, axis=1)
        df["p2"] = df.apply(pd2, axis=1)
        df["B_D2_DOT"] = df.apply(
            fd2_dot_p,
            axis=1,
        )
        df["D2_DIRA_ORIVX_FIXED"] = df.apply(
            ns2, axis=1, args=["X_myFD_u", "p2", "B_D2_DOT"], result_type="expand"
        )
        df["D2stmD_M"] = df["D2st_M"] - df["D2_M"]

    df['B_dtf_M'] = df.B_dtf_M.apply(lambda x: x[0])
    df['B_DTF_M'] = df['B_dtf_M']

    # ["B_D2_DOT"]/(df["p2"]*df["X_myFD_u"])
    # df[["X_myFD", "X_myFDERR", "D2_FDCHI2_ORIVX"]] = df.apply(
    #     ns, args=["X_myFD_u"], axis=1, result_type="expand"
    # )

    # -- write output
    # print(f"writing to '{outfilename}'...")
    # if mc_flag == 1:
    #     columns_to_read_out = columns_to_read_out + mc_to_read_in + st_id_l
    # if mc_flag == 0:
    #     columns_to_read_out = columns_to_read_out + data_to_read_in
    df.to_root(outfilename, key=outtreename, store_index=False)
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
