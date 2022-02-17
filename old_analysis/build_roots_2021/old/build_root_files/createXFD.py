def CreateXFD(infilename, intreename, outfilename, outtreename, mc_flag = 0):
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
        # print("ye")
        # "correlate values using covariance matrix"
        # print([
        #     float(getattr(row, f"{part}_{xbase}_X")),
        #     float(getattr(row, f"{part}_{xbase}_Y")),
        #     float(getattr(row, f"{part}_{xbase}_Z")),
        # ]),
        # print([
        #     getattr(row, f"{part}_{xbase}_COV_")[0:3],
        #     # getattr(row, f"{part}_{xbase}_COV_")[3:6],
        #     # getattr(row, f"{part}_{xbase}_COV_")[6:9],
        # ])

        if mc_flag == 0:
            return correlated_values(
                [
                    float(getattr(row, f"{part}_{xbase}_X")),
                    float(getattr(row, f"{part}_{xbase}_Y")),
                    float(getattr(row, f"{part}_{xbase}_Z")),
                ],
                [
                    getattr(row, f"{part}_{xbase}_COV_")[0:3],
                    getattr(row, f"{part}_{xbase}_COV_")[3:6],
                    getattr(row, f"{part}_{xbase}_COV_")[6:9],
                ],
            )
        if mc_flag == 1:
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

    def FD(row):
        "calculate distance from myPV to myENDV"
        return sqrt(
                (row.B_myENDV_X - row.D2_myENDV_X) ** 2
            + (row.B_myENDV_Y - row.D2_myENDV_Y) ** 2
            + (row.B_myENDV_Z - row.D2_myENDV_Z) ** 2
        )

    def p2(row):
        "calculate distance from myPV to myENDV"
        return sqrt((row.D2_PX) ** 2 + (row.D2_PY) ** 2 + (row.D2_PZ) ** 2)

    def fd_dot_p(row):
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

    columns_to_read_in = (
        [
            "B_dtf_M",
            "D1_M",
            "D2_M",
            "D2st_M",
            "D1_DIRA_ORIVX",
            "D2_DIRA_ORIVX",
            "D1_FDCHI2_ORIVX",
            "eventNumber",
            "runNumber",
            "nCandidate",

        ]
        + [
            f"{part}_ENDVERTEX_{coord}"
            for coord in ("X", "Y", "Z", "COV_")
            for part in ("B", "D2")
        ]
        + [
            "D1_PE",
            "D1_PX",
            "D1_PY",
            "D1_PZ",
            "D2_PE",
            "D2_PX",
            "D2_PY",
            "D2_PZ",
            "KSTH1_PE",
            "KSTH1_PX",
            "KSTH1_PY",
            "KSTH1_PZ",
            "KSTH2_PE",
            "KSTH2_PX",
            "KSTH2_PY",
            "KSTH2_PZ",
        ]
    )
    st_id_l = [
        "D1H2_TRUEID",
        "D1H1_TRUEID",
        "KSTH2_TRUEID",
        "D2H2_TRUEID",
        "KSTH1_TRUEID",
        "D2_TRUEID",
        "D2st_TRUEID",
        "D1_TRUEID",
        "D2H1_TRUEID",
        "B_TRUEID",
        "KST_TRUEID",
        "D2H3_TRUEID",
        "KST_MC_GD_MOTHER_KEY",
        "KSTH2_MC_MOTHER_KEY",
        "D1_MC_GD_GD_MOTHER_KEY",
        "KST_MC_MOTHER_KEY",
        "D2st_MC_MOTHER_KEY",
        "D2H3_MC_GD_GD_MOTHER_KEY",
        "D2_MC_MOTHER_KEY",
        "D2_MC_GD_MOTHER_KEY",
        "D2H1_MC_GD_GD_MOTHER_KEY",
        "KSTH1_MC_GD_GD_MOTHER_KEY",
        "D2st_MC_GD_MOTHER_KEY",
        "D1H2_MC_GD_GD_MOTHER_KEY",
        "D1_MC_GD_MOTHER_KEY",
        "KSTH1_MC_MOTHER_KEY",
        "KSTH2_MC_GD_GD_MOTHER_KEY",
        "D2H3_MC_MOTHER_KEY",
        "B_MC_MOTHER_KEY",
        "D1H1_MC_GD_MOTHER_KEY",
        "D1H2_MC_MOTHER_KEY",
        "D1H1_MC_GD_GD_MOTHER_KEY",
        "D1_MC_MOTHER_KEY",
        "D2_MC_GD_GD_MOTHER_KEY",
        "KST_MC_GD_GD_MOTHER_KEY",
        "KSTH1_MC_GD_MOTHER_KEY",
        "D1H1_MC_MOTHER_KEY",
        "D2H1_MC_MOTHER_KEY",
        "B_MC_GD_MOTHER_KEY",
        "D2H2_MC_GD_GD_MOTHER_KEY",
        "D2H1_MC_GD_MOTHER_KEY",
        "B_MC_GD_GD_MOTHER_KEY",
        "D2st_MC_GD_GD_MOTHER_KEY",
        "D2H2_MC_GD_MOTHER_KEY",
        "D2H2_MC_MOTHER_KEY",
        "KSTH2_MC_GD_MOTHER_KEY",
        "D2H3_MC_GD_MOTHER_KEY",
        "D1H2_MC_GD_MOTHER_KEY",
        "KSTH2_MC_GD_GD_MOTHER_ID",
        "D2_MC_MOTHER_ID",
        "D1_MC_MOTHER_ID",
        "D1_MC_GD_GD_MOTHER_ID",
        "B_MC_GD_GD_MOTHER_ID",
        "KSTH1_MC_GD_GD_MOTHER_ID",
        "KSTH2_MC_MOTHER_ID",
        "D1_MC_GD_MOTHER_ID",
        "D1H2_MC_GD_MOTHER_ID",
        "D2H3_MC_MOTHER_ID",
        "D2_MC_GD_GD_MOTHER_ID",
        "D2st_MC_GD_MOTHER_ID",
        "D2st_MC_GD_GD_MOTHER_ID",
        "KSTH2_MC_GD_MOTHER_ID",
        "D1H2_MC_MOTHER_ID",
        "D2H1_MC_MOTHER_ID",
        "D2H1_MC_GD_GD_MOTHER_ID",
        "KST_MC_GD_MOTHER_ID",
        "D2H3_MC_GD_GD_MOTHER_ID",
        "D2H2_MC_GD_GD_MOTHER_ID",
        "D2st_MC_MOTHER_ID",
        "KST_MC_GD_GD_MOTHER_ID",
        "B_MC_GD_MOTHER_ID",
        "D2H1_MC_GD_MOTHER_ID",
        "D2H3_MC_GD_MOTHER_ID",
        "KST_MC_MOTHER_ID",
        "D1H2_MC_GD_GD_MOTHER_ID",
        "D1H1_MC_GD_GD_MOTHER_ID",
        "KSTH1_MC_GD_MOTHER_ID",
        "KSTH1_MC_MOTHER_ID",
        "B_MC_MOTHER_ID",
        "D2H2_MC_MOTHER_ID",
        "D2_MC_GD_MOTHER_ID",
        "D1H1_MC_GD_MOTHER_ID",
        "D1H1_MC_MOTHER_ID",
        "D2H2_MC_GD_MOTHER_ID",
    ]


    mc_to_read_in = [
            'D1_TRUEP_E',
            'D2_TRUEP_E',
            'KST_TRUEP_E',
            'D1_TRUEP_X',
            'D2_TRUEP_X',
            'KST_TRUEP_X',
            'D1_TRUEP_Y',
            'D2_TRUEP_Y',
            'KST_TRUEP_Y',
            'D1_TRUEP_Z',
            'D2_TRUEP_Z',
            'KST_TRUEP_Z',
            'KST_M',
            'D1H1_ProbNNk',
            'D2H1_ProbNNk',
            'KSTH1_ProbNNk',
            'B_L0Global_TOS',
            'B_Hlt1Global_TOS',
            'B_Hlt2Global_TOS',
            'B_L0HadronDecision_TOS',
            'B_L0MuonDecision_TOS',
            'B_L0ElectronDecision_TOS',
            'B_L0PhotonDecision_TOS',
            'B_L0HadronDecision_TIS',
            'B_L0MuonDecision_TIS',
            'B_L0ElectronDecision_TIS',
            'B_L0PhotonDecision_TIS',
            "RD_org_eventNumber",
            "RD_org_runNumber",
            "B_Hlt1TrackMVADecision_TOS",
            "B_Hlt1TwoTrackMVADecision_TOS",
            "B_Hlt2Topo2BodyDecision_TOS",
            "B_Hlt2Topo3BodyDecision_TOS",
            "B_Hlt2Topo4BodyDecision_TOS",
            "B_Hlt1TrackMVADecision_TIS",
            "B_Hlt1TwoTrackMVADecision_TIS",
            "B_Hlt2Topo2BodyDecision_TIS",
            "B_Hlt2Topo3BodyDecision_TIS",
            "B_Hlt2Topo4BodyDecision_TIS"
    ]
    data_to_read_in = [
            'KST_M',
            'D1H1_ProbNNk',
            'D2H1_ProbNNk',
            'KSTH1_ProbNNk',
            'B_L0Global_TOS',
            'B_Hlt1Global_TOS',
            'B_Hlt2Global_TOS',
            'B_L0HadronDecision_TOS',
            'B_L0MuonDecision_TOS',
            'B_L0ElectronDecision_TOS',
            'B_L0PhotonDecision_TOS',
            'B_L0HadronDecision_TIS',
            'B_L0MuonDecision_TIS',
            'B_L0ElectronDecision_TIS',
            'B_L0PhotonDecision_TIS',
    ]
    columns_to_read_out = [
        "D1_M",
        "D2_M",
        "D2st_M",
        "DstmD_M",
        "D1_DIRA_ORIVX",
        "D2_DIRA_ORIVX",
        "D1_FDCHI2_ORIVX",
        "D2_DIRA_ORIVX_FIXED",
        "D1_PE",
        "D1_PX",
        "D1_PY",
        "D1_PZ",
        "D2_PE",
        "D2_PX",
        "D2_PY",
        "D2_PZ",
        "KSTH1_PE",
        "KSTH1_PX",
        "KSTH1_PY",
        "KSTH1_PZ",
        "KSTH2_PE",
        "KSTH2_PX",
        "KSTH2_PY",
        "KSTH2_PZ",
        'B_dtf_M',
        # 'B_DTF_M'
    ]
    # -- read input
    # print(f"reading '{intreename}' from '{infilename}'...")
    if mc_flag == 1:
        columns_to_read_in = columns_to_read_in + mc_to_read_in + st_id_l
    if mc_flag == 0:
        columns_to_read_in = columns_to_read_in + data_to_read_in
    df = rp.read_root(
        infilename,
        intreename,
        columns_to_read_in,
        # flatten=["B_dtf_M"]
    )
    # df = df[1:]
    # -- transform columns
    # print("applying transformations...")
    l = df.columns
    # for i in l:
    #     if "MOTHER_ID" in i or "TRUEID" in i or "MOTHER_KEY" in i:
    #         print(i)
    df[["B_myENDV_X", "B_myENDV_Y", "B_myENDV_Z"]] = df.apply(
        corrval, args=["ENDVERTEX", "B"], axis=1, result_type="expand"
    )
    df[["D2_myENDV_X", "D2_myENDV_Y", "D2_myENDV_Z"]] = df.apply(
        corrval, args=["ENDVERTEX", "D2"], axis=1, result_type="expand"
    )

    df["X_myFD_u"] = df.apply(FD, axis=1)

    df[["p2"]] = df.apply(p2, axis=1)

    df[["B_D2_DOT"]] = df.apply(
        fd_dot_p,
        axis=1,
    )
    df['B_dtf_M'] = df.B_dtf_M.apply(lambda x: x[0])
    # df.drop(columns=['B_DTF_M'])
    df['B_DTF_M'] = df['B_dtf_M']
    df[["D2_DIRA_ORIVX_FIXED"]] = df.apply(
        ns2, axis=1, args=["X_myFD_u", "p2", "B_D2_DOT"], result_type="expand"
    )

    # ["B_D2_DOT"]/(df["p2"]*df["X_myFD_u"])
    # df[["X_myFD", "X_myFDERR", "D2_FDCHI2_ORIVX"]] = df.apply(
    #     ns, args=["X_myFD_u"], axis=1, result_type="expand"
    # )

    df["DstmD_M"] = df["D2st_M"] - df["D2_M"]
    # -- write output
    # print(f"writing to '{outfilename}'...")
    if mc_flag == 1:
        columns_to_read_out = columns_to_read_out + mc_to_read_in + st_id_l
    if mc_flag == 0:
        columns_to_read_out = columns_to_read_out + data_to_read_in
    df[columns_to_read_out].to_root(outfilename, key=outtreename, store_index=False)
    # print("done")


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
