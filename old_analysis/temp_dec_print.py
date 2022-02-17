dm_decay        = "${D1}(D- -> ${D1H1}K+ ${D1H2}pi- ${D1H3}pi-)"
dm_decay_cc     = "${D1}(D+ -> ${D1H1}K- ${D1H2}pi+ ${D1H3}pi+)"

dp_decay        = "${D2}(D+ -> ${D2H1}K- ${D2H2}pi+ ${D2H3}pi+)"
dp_decay_cc     = "${D2}(D- -> ${D2H1}K+ ${D2H2}pi- ${D2H3}pi-)"

d0bar_decay     = "${D1}(D0 -> ${D1H1}K+ ${D1H2}pi-)"
d0bar_decay_cc  = "${D1}(D0 -> ${D1H1}K- ${D1H2}pi+)"

d0_decay        = "${D2}(D0 -> ${D2H1}K- ${D2H2}pi+)"
d0_decay_cc     = "${D2}(D0 -> ${D2H1}K+ ${D2H2}pi-)"

dstm_decay      = "${{D1st}}(D*(2010)- -> {} ${{D1H3}}pi-)".format(d0bar_decay)
dstm_decay_cc   = "${{D1st}}(D*(2010)+ -> {} ${{D1H3}}pi+)".format(d0bar_decay_cc)

dstp_decay      = "${{D2st}}(D*(2010)+ -> {} ${{D2H3}}pi+)".format(d0_decay)
dstp_decay_cc   = "${{D2st}}(D*(2010)- -> {} ${{D2H3}}pi-)".format(d0_decay_cc)

dsm_decay       = "${D1}(D- -> ${D1H1}K- ${D1H2}K+ ${D1H3}pi-)"
dsm_decay_cc    = "${D1}(D+ -> ${D1H1}K+ ${D1H2}K- ${D1H3}pi+)"

d04_decay       = "${D2}(D0 -> ${D2H1}K- ${D2H2}pi+ ${D2H3}pi+ ${D2H4}pi-)"
d04_decay_cc    = "${D2}(D0 -> ${D2H1}K+ ${D2H2}pi- ${D2H3}pi- ${D2H4}pi+)"

kst_decay       = "${KST}(K*(892)0 -> ${KSTH1}K+ ${KSTH2}pi-)"
kst_decay_cc    = "${KST}(K*(892)~0 -> ${KSTH1}K- ${KSTH2}pi+)"

filtercode = "(NINTREE( (ABSID=='pi+') & (PIDK >= 0) ) == 0) & (NINTREE( (ABSID=='K+') & (PIDK < 4) ) == 0) & (NINTREE( (ABSID=='p+') & (PIDp < 0) ) == 0) & (NINTREE( (ABSID=='K*(892)0') & ((M>1050) | (M<750))) == 0) & (NINGENERATION( (VFASPF(VCHI2PDOF) > 5), 1) == 0) & (MINTREE( ((ABSID=='D+')|(ABSID=='D0')), VFASPF(VZ)) - VFASPF(VZ) > -2.0 )"

rec_Descriptors = {
    "Z_m_p"       : "(${{B}}(B0 ->{}{}{}))".format(dm_decay,dp_decay,kst_decay,dm_decay_cc,dp_decay_cc,kst_decay_cc),

    "Z_mst_p"     : "(${{B}}(B0 ->{}{}{}))".format(dstm_decay,dp_decay,kst_decay,dstm_decay_cc,dp_decay_cc,kst_decay_cc),
    "Z_m_pst"     : "(${{B}}(B0 ->{}{}{}))".format(dm_decay,dstp_decay,kst_decay,dm_decay_cc,dstp_decay_cc,kst_decay_cc),
    "Z_mst_pst"   : "(${{B}}(B0 ->{}{}{}))".format(dstm_decay,dstp_decay,kst_decay,dstm_decay_cc,dstp_decay_cc,kst_decay_cc),

    "Z_z_z"       : "(${{B}}(B0 ->{}{}{}))".format(d0bar_decay,d0_decay,kst_decay,d0bar_decay_cc,d0_decay_cc,kst_decay_cc),

    "P_z_p"       : "(${{B}}(B+ ->{}{}{}))".format(d0bar_decay,dp_decay,kst_decay,d0bar_decay_cc,dp_decay_cc,kst_decay_cc),
    "P_z_pst"     : "(${{B}}(B+ ->{}{}{}))".format(d0bar_decay,dstp_decay,kst_decay,d0bar_decay_cc,dstp_decay_cc,kst_decay_cc),

    "M_m_z"       : "(${{B}}(B- ->{}{}{}))".format(dm_decay, d0_decay, kst_decay, dm_decay_cc, d0_decay_cc, kst_decay_cc),
    "M_mst_z"     : "(${{B}}(B- ->{}{}{}))".format(dm_decay, d0_decay, kst_decay, dm_decay_cc, d0_decay_cc, kst_decay_cc),

    "Zs_ms_p"      : "(${{B}}(B0 ->{}{}{}))".format(dsm_decay,dp_decay,kst_decay,dsm_decay_cc,dp_decay_cc,kst_decay_cc),

    "norm7" : "(${{B}}(B+ ->{} {} ${{K}} K+) || ${{B}}(B- ->{} {} ${{K}} K-))".format(d0bar_decay,d04_decay,d0bar_decay_cc,d04_decay_cc),
    "norm8" : "(${{B}}(B0 ->{} {} ${{K}} K+) || ${{B}}(B0 ->{} {} ${{K}} K-))".format(dm_decay,d04_decay,dm_decay_cc,d04_decay_cc),
}

for i in rec_Descriptors:
    print(rec_Descriptors[i])
