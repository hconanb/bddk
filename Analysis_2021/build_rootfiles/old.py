def mult_can_checks(rdf):

    npy_temp = rdf.AsNumpy(columns=["eventNumber","runNumber"])
    df_temp = pandas.DataFrame(npy_temp)
    unique_events = df_temp[~df_temp.duplicated()].value_counts()

    n_unique_events = unique_events.size
    n_unique_events_err = math.sqrt(n_unique_events)
    n_unique_events_ufloat = ufloat(n_unique_events, n_unique_events_err)

    nodups = df_temp[~df_temp.duplicated(keep=False)].value_counts()
    dups = df_temp[df_temp.duplicated(keep=False)].value_counts()
    both = nodups.append(dups)
    can_per_event = both.mean()
    can_per_event_err = both.std()
    cen_per_event_ufloat = ufloat(can_per_event, can_per_event_err)

    return n_unique_events_ufloat, can_per_event_ufloat
    
