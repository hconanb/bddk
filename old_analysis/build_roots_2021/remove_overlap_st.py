def zz_olapfit(name, year, trigger):

    file_base = ROOT.TFile(f"/mnt/c/Users/Harris/Desktop/rootfiles/2021_filtered_root/{type_flag}/{spec_name}_{year}_{trigger_flag}_spectrum_filtered.root")

    tree_zz = file_base.Get(f"DecayTreeTuple")

    st_name = spec_name.replace("_zz_","_st_")

    tree_st = file_base.Get(f"DecayTreeTuple")

    rdf_zz = RDF(tree_zz)
    rdf_st = RDF(tree_st)

    np_zz = rdf_zz.AsNumpy()
    np_st = rdf_st.AsNumpy()

    zz_en_base = np.char.mod('%d', np_zz["eventNumber"])
    zz_nCan_base = np.char.mod('%d', np_zz["runNumber"])

    st_en_base = np.char.mod('%d', np_st["eventNumber"])
    st_nCan_base = np.char.mod('%d', np_st["runNumber"])

    zz_all = np.core.defchararray.add(np.core.defchararray.add(zz_en_base, "_"),zz_nCan_base)
    st_all = np.core.defchararray.add(np.core.defchararray.add(st_en_base, "_"),st_nCan_base)

    olap_list = np.intersect1d(zz_all,st_all)
    fline = "B_M > 0"
    for olap_entry in olap_list:
        eno = olap_entry.split("_")[0]
        ncan = olap_entry.split("_")[1]
        fline = f"{fline} && (eventNumber != {eno} || (eventNumber == {eno} && nCandidate != {ncan}))"
         # if old_evvent.eventNumer != eno or (old_event.eventNumber = eno and old_event.nCandidate != ncan):
    rdf_zz_next = rdf_zz.Filter(fline)

    print(len(zz_en_base))
    print(len(olap_list))
    print(rdf_zz_next.Count().GetValue())

    clist = rdf_zz_next.GetColumnNames()

    if os.path.exists(file_name.replace("ntuple.root", "ntuple_zz_fix.root")):
        os.remove(file_name.replace("ntuple.root", "ntuple_zz_fix.root"))
        print("First deleting old ", file_name.replace("ntuple.root", "ntuple_zz_fix.root"))
    else:
        print("making ", file_name.replace("ntuple.root", "ntuple_zz_fix.root"))
    rdfsnap = rdf_zz_next.Snapshot(f"{spec_name}/DecayTreeTuple", file_name.replace("ntuple.root", "ntuple_zz_fix.root"), clist, opts)
#
