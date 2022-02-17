from Gaudi.Configuration import *
from Configurables import *
from GaudiConf import IOHelper
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence, DataOnDemand
from PhysConf.Filters import LoKi_Filters
from PhysConf.Selections import FilterSelection, PrintSelection
from DecayTreeTuple.Configuration import *
from MicroDSTConf.TriggerConfUtils import configureL0AndHltDecoding
from Configurables import TrackSmearState as SMEAR

import glob
import sys, os
sys.path.append(os.getcwd())
from steering import *

def addtopinfo(branch):
    loki_top = branch.addTupleTool("LoKi::Hybrid::TupleTool/loki_top")
    loki_top.Variables = {
        "PHI" : "PHI",
        "ETA" : "ETA",
        "MIPCHI2DV":"MIPCHI2DV(PRIMARY)",
        "BPVVDCHI2":"BPVVDCHI2",
        "BPVDIRA":"BPVDIRA",
        "BPVLTIME":"BPVLTIME()",
        "BPVZ":"BPV(VZ)",
        "VCHI2PDOF":"VFASPF(VCHI2PDOF)",
    }

def addmidinfo(branch):
    loki_mid = branch.addTupleTool("LoKi::Hybrid::TupleTool/loki_mid")
    loki_mid.Variables = {
        "PHI" : "PHI",
        "ETA" : "ETA",
        "MIPCHI2DV":"MIPCHI2DV(PRIMARY)",
        "VCHI2PDOF":"VFASPF(VCHI2PDOF)"
    }

def addtrkinfo(branch):
    loki_trk = branch.addTupleTool("LoKi::Hybrid::TupleTool/loki_trk")
    loki_trk.Variables = {
        "PHI" : "PHI",
        "ETA" : "ETA",
        "MIPCHI2DV":"MIPCHI2DV(PRIMARY)",
        "TRCHI2DOF":"TRCHI2DOF",
        "TRGHOSTPROB":"TRGHOSTPROB",
    }

def addtools (t, spectrum):
    ToolList = [
        'TupleToolGeometry',
        'TupleToolPrimaries',
        'TupleToolEventInfo',
        'TupleToolAngles',
        'TupleToolL0Calo',
        'TupleToolRecoStats',
        'TupleToolKinematic',
        'TupleToolTrackInfo',
        'TupleToolPropertime',
        'TupleToolPid',
        ]
    triglist = [
        "L0HadronDecision",
        "L0MuonDecision",
        "L0ElectronDecision",
        "L0PhotonDecision",
        'Hlt1TrackMVADecision',
        'Hlt1TwoTrackMVADecision',
        "Hlt2Topo2BodyDecision",
        "Hlt2Topo3BodyDecision",
        "Hlt2Topo4BodyDecision",
        "Hlt1TrackMVALooseDecision",
        "Hlt1TwoTrackMVALooseDecision",
        "Hlt1CalibHighPTLowMultTrksDecision",
        "Hlt1IncPhiDecision",
        "Hlt2ForwardDecision",
        "Hlt2XcMuXForTauB2XcFakeMuDecision",
        ]
    if mytype != 'Data':
        truth = t.addTupleTool("TupleToolMCTruth")
        truth.addTupleTool("MCTupleToolHierarchy")
        t.addTupleTool("TupleToolMCBackgroundInfo")
        t.addTupleTool("MCTupleToolRedecay")


    ####DTF Constrain
    t.ToolList += ToolList
    t.B.addTupleTool('TupleToolDecayTreeFitter/dtf')
    t.B.addTupleTool('TupleToolDecayTreeFitter/dtf_c')

    t.B.dtf.constrainToOriginVertex = False
    t.B.dtf.Verbose = True

    t.B.dtf_c.constrainToOriginVertex = True
    t.B.dtf_c.Verbose = True

    if spectrum == 'Z_m_p':
        t.B.dtf.daughtersToConstrain += ['D+', 'D-']
        t.B.dtf_c.daughtersToConstrain += ['D+', 'D-']

    if spectrum == "Z_z_z" or spectrum == "norm7":
        t.B.dtf.daughtersToConstrain += ['D0', 'D~0']
        t.B.dtf_c.daughtersToConstrain += ['D0', 'D~0']

    if spectrum == "P_z_p" or spectrum == "M_m_z" or spectrum == "norm8":
        t.B.dtf.daughtersToConstrain += ['D+', 'D-', 'D0', 'D~0']
        t.B.dtf_c.daughtersToConstrain += ['D+', 'D-', 'D0', 'D~0']

    if spectrum == "P_z_pst" or spectrum == "M_mst_z" or spectrum == "Z_mst_pst":
        t.B.dtf.daughtersToConstrain += ['D*(2010)+', 'D*(2010)-', 'D0', 'D~0']
        t.B.dtf_c.daughtersToConstrain += ['D*(2010)+', 'D*(2010)-', 'D0', 'D~0']

    if spectrum == 'Zs_sm_p':
        t.B.dtf.Substitutions = {
            "Beauty -> Charm ^(D+ -> K+ K- pi+) Meson" : "D_s+",
            "Beauty -> Charm ^(D- -> K- K+ pi-) Meson" : "D_s-",
        }
        t.B.dtf.daughtersToConstrain += ['D_s+', 'D_s-', 'D+', 'D-']

        t.B.dtf_c.Substitutions = {
            "Beauty -> Charm ^(D+ -> K+ K- pi+) Meson" : "D_s+",
            "Beauty -> Charm ^(D- -> K- K+ pi-) Meson" : "D_s-",
        }
        t.B.dtf_c.daughtersToConstrain += ['D_s+', 'D_s-', 'D+', 'D-']


    if spectrum == 'Z_mst_p' or spectrum == 'Z_m_pst':
        t.B.dtf.daughtersToConstrain += ['D*(2010)+', 'D*(2010)-', 'D+', 'D-', 'D0', 'D~0']
        t.B.dtf_c.daughtersToConstrain += ['D*(2010)+', 'D*(2010)-', 'D+', 'D-', 'D0', 'D~0']

    trig = t.B.addTupleTool("TupleToolTISTOS")
    trig.TriggerList = triglist
    trig.Verbose = True

    trig = t.D1.addTupleTool("TupleToolTISTOS")
    trig.TriggerList = triglist
    trig.Verbose = True

    trig = t.D2.addTupleTool("TupleToolTISTOS")
    trig.TriggerList = triglist
    trig.Verbose = True

    #adds info for tuples with a Kst
    if 'norm' not in spectrum:
        if '11198030' not in spectrum and '11198007' not in spectrum:
            trig = t.KST.addTupleTool("TupleToolTISTOS")
            trig.TriggerList = triglist
            trig.Verbose = True
            addmidinfo( t.KST )
            addtrkinfo( t.KSTH1 )
            addtrkinfo( t.KSTH2 )

    sinfo_str = "INFO( {}, -1.0)"
    ## angles are just 1.5, 1.7, 1 every time, so don't bother saving them

    addtopinfo( t.B )
    for n in [1,2,3]:
        for j,v in enumerate([ "MULT", "PTASYM" ]):
            t.B.loki_top.Variables["CONE{}_{}".format(v,n)] = sinfo_str.format(9000 + 10*n + j)
    #print t.B.loki_top.Variables
    addmidinfo( t.D1 )
    addmidinfo( t.D2 )
    addtrkinfo( t.D1H1 )
    addtrkinfo( t.D1H2 )
    addtrkinfo( t.D2H1)
    addtrkinfo( t.D2H2)
    try:
        addmidinfo( t.D1st)
        addmidinfo( t.D2st)
        addtrkinfo( t.D1H3)
        addtrkinfo( t.D2H3)
        addtrkinfo( t.D1H4)
        addtrkinfo( t.D2H4)
    except:
        pass

    t.addTool(TupleToolSubMass, name="TupleToolSubMass")
    t.TupleToolSubMass.SetMax = 4
    t.TupleToolSubMass.EndTreePIDs = [411, -411, 421, -421]
    t.ToolList += ["TupleToolSubMass"]

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
    "Z_m_p"       : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dm_decay,dp_decay,kst_decay,dm_decay_cc,dp_decay_cc,kst_decay_cc),

    "Z_mst_p"     : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dstm_decay,dp_decay,kst_decay,dstm_decay_cc,dp_decay_cc,kst_decay_cc),
    "Z_m_pst"     : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dm_decay,dstp_decay,kst_decay,dm_decay_cc,dstp_decay_cc,kst_decay_cc),
    "Z_mst_pst"   : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dstm_decay,dstp_decay,kst_decay,dstm_decay_cc,dstp_decay_cc,kst_decay_cc),

    "Z_z_z"       : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(d0bar_decay,d0_decay,kst_decay,d0bar_decay_cc,d0_decay_cc,kst_decay_cc),

    "P_z_p"       : "(${{B}}(B+ ->{}{}{}) || ${{B}}(B- ->{}{}{}))".format(d0bar_decay,dp_decay,kst_decay,d0bar_decay_cc,dp_decay_cc,kst_decay_cc),
    "P_z_pst"     : "(${{B}}(B+ ->{}{}{}) || ${{B}}(B- ->{}{}{}))".format(d0bar_decay,dstp_decay,kst_decay,d0bar_decay_cc,dstp_decay_cc,kst_decay_cc),

    "M_m_z"       : "(${{B}}(B- ->{}{}{}) || ${{B}}(B+ ->{}{}{}))".format(dm_decay, d0_decay, kst_decay, dm_decay_cc, d0_decay_cc, kst_decay_cc),
    "M_mst_z"     : "(${{B}}(B- ->{}{}{}) || ${{B}}(B+ ->{}{}{}))".format(dstm_decay, d0_decay, kst_decay, dstm_decay_cc, d0_decay_cc, kst_decay_cc),

    "Zs_sm_p"      : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dsm_decay,dp_decay,kst_decay,dsm_decay_cc,dp_decay_cc,kst_decay_cc),

    "norm7" : "(${{B}}(B+ ->{} {} ${{K}} K+) || ${{B}}(B- ->{} {} ${{K}} K-))".format(d0bar_decay,d04_decay,d0bar_decay_cc,d04_decay_cc),
    "norm8" : "(${{B}}(B0 ->{} {} ${{K}} K+) || ${{B}}(B0 ->{} {} ${{K}} K-))".format(dm_decay,d04_decay,dm_decay_cc,d04_decay_cc),
}

strip_dict = {
    'Z_m_p'       : 'B02DDKstBeauty2CharmLine',
    'Zs_sm_p'     : 'B02DDKstBeauty2CharmLine',

    'Z_mst_p'    : 'B02DstDKstBeauty2CharmLine',
    'Z_m_pst'    : 'B02DstDKstBeauty2CharmLine',

    'Z_mst_pst' : 'B02DstDstKstBeauty2CharmLine',

    'Z_z_z'    : 'B02D0D0KstD02HHD02HHBeauty2CharmLine',

    'P_z_p'     : 'B2DD0KstBeauty2CharmLine',
    'M_m_z'     : 'B2DD0KstBeauty2CharmLine',

    'P_z_pst'    : 'B2DstD0KstBeauty2CharmLine',
    'M_mst_z'    : 'B2DstD0KstBeauty2CharmLine',

    'norm7' : 'B2D0D0KD02HHD02K3PiBeauty2CharmLine',
    'norm8' : 'B02D0DKD02K3PiBeauty2CharmLine',
    }

class mytuple():
    def __init__(info, name, id, spectrum):
        info.name = name
        info.id = id
        info.spectrum = spectrum

t = []

t.append(mytuple("01", "11198006", ["Z_m_p"]))

t.append(mytuple("02", "11198400", ["Z_m_p"]))
t.append(mytuple("02", "11198005", ["P_z_p"]))

t.append(mytuple("04", "11198401", ["Z_m_p"]))
t.append(mytuple("04", "11198410", ["P_z_p"]))

# t.append(mytuple("04", "11198022", ["Z_z_z", "P_z_p", "M_m_z", "P_z_pst", "M_mst_z", "Z_mst_pst"]))
t.append(mytuple("04", "11198023", ["Z_z_z", "P_z_p", "M_m_z", "P_z_pst"]))

t.append(mytuple("05", "12197023", ["P_z_p"]))
t.append(mytuple("06", "12197410", ["P_z_p"]))

t.append(mytuple("07", "12197400", ["P_z_p"]))
# t.append(mytuple("07", "12197024", ["Z_z_z", "P_z_pst"]))
t.append(mytuple("07", "12197045", ["Z_z_z", "P_z_pst"]))

t.append(mytuple("08", "12197401", ["P_z_p"]))
# t.append(mytuple("08", "12197422", ["Z_z_z", "P_z_pst"]))
t.append(mytuple("08", "12197423", ["Z_z_z", "P_z_pst"]))

t.append(mytuple("09", "11196019", ["Z_z_z"]))
t.append(mytuple("10", "11196413", ["Z_z_z"]))
t.append(mytuple("12", "11196414", ["Z_z_z"]))

t.append(mytuple("13", "13198040", ["Zs_sm_p"]))
t.append(mytuple("14", "13198200", ["Zs_sm_p"]))
t.append(mytuple("15", "13198400", ["Zs_sm_p"]))
t.append(mytuple("16", "13198600", ["Zs_sm_p"]))

t.append(mytuple("norm7", "12197008", ["norm7"]))
t.append(mytuple("norm8", "11198007", ["norm8"]))


#12197045

dl = ["Z_m_p",
      "Zs_sm_p",
      "Z_mst_pst",
      "Z_z_z",
      "P_z_p",
      "M_m_z",
      "P_z_pst",
      "norm7",
      "norm8"]

t.append(mytuple("data", "Data", dl))

def getthistuple(type):
    for i in t:
        if i.id == mytype:
            return i

tt = getthistuple(mytype)

DaVinci().InputType = myst
DaVinci().RootInTES = "/Event/{}".format(mystream)

if tt.id != "Data":
    MCDTT = MCDecayTreeTuple("MCDecayTreeTuple")
    MCDTT_2 = MCDecayTreeTuple("MCDecayTreeTuple2")

    if tt.id == "11198007" or tt.id == "12197008":
        MCDTT.setDescriptorTemplate("${B}[Beauty => ${C1}Charm ${C2}Charm ${K}K+]CC")
    if tt.id != "11198007" and tt.id != "12197008":
        MCDTT.setDescriptorTemplate("${B}[Beauty => ${C1}Charm ${C2}Charm ${Kst}K*(892)0]CC")
        tempd = "(${{B}}(Beauty ->{}{}{}) || ${{B}}(Beauty ->{}{}{}))".format(dm_decay,dp_decay,kst_decay,dm_decay_cc,dp_decay_cc,kst_decay_cc).replace('->','=>')
        MCDTT_2.setDescriptorTemplate(tempd)

    MCDTT.ToolList+=["MCTupleToolKinematic","MCTupleToolHierarchy"]
    MCDTT_2.ToolList+=["MCTupleToolKinematic","MCTupleToolHierarchy"]

    DaVinci().UserAlgorithms.append(MCDTT)
    DaVinci().UserAlgorithms.append(MCDTT_2)


    smear = SMEAR('smear')
    DaVinci().UserAlgorithms.append(smear)

for spectrum in tt.spectrum:
    for n in [1,2,3]:
        for j,v in enumerate([ "ANGLE", "MULT", "PTASYM" ]):
            rinfo_str = "v{var}{num} = SINFO( {code}, RELINFO('/Event/{stream}/Phys/{{}}/P2ConeVar{num}','CONE{var}',-1.0 ), True)"
            preamble = []
            sinfo_codes = []
            preamble.append( rinfo_str.format(var = v, num = n, stream = mystream, code = 9000 + 10*n + j) )
            sinfo_codes.append( "(v{var}{num} > -100.0)".format(var=v, num=n) )
    sinfo_code = " & ".join(sinfo_codes)
    loc = AutomaticData(Location = '/Phys/{}/Particles'.format(strip_dict[spectrum]))
    relsel = FilterSelection("relsel_{}".format(spectrum), loc, Preambulo= [rstr.format(strip_dict[spectrum]) for rstr in preamble], Code=sinfo_code)
    DaVinci().UserAlgorithms.append(relsel)
    if spectrum == "M_m_z":
        switchkst_bu = SubstitutePID( 'switchkst_bu',
                                      Code = "DECTREE('(B+ -> D+ D0 K*(892)0) || (B- ->D- D0 K*(892)~0)')",
                                      Substitutions = {
                                          'B+ -> D+ D0 (Meson -> X+ ^pi-)' : 'K-',
                                          'B+ -> D+ D0 (Meson -> ^K+ X-)' : 'pi+',
                                          'B+ -> D+ D0 ^K*(892)0' : 'K*(892)~0',
                                          'B- -> D- D0 (Meson -> X- ^pi+)' : 'K+',
                                          'B- -> D- D0 (Meson -> ^K- X+)' : 'pi-',
                                          'B- -> D- D0 ^K*(892)~0' : 'K*(892)0',
                                      })
        switchkst_bu.MaxChi2PerDoF = -1
        switchkst_sel_bu = Selection(
            'switchkst_sel_bu',
            Algorithm = switchkst_bu,
            RequiredSelections=[relsel]
        )
        switchSeq_bu = SelectionSequence('switchSeq_bu', TopSelection=switchkst_sel_bu)
        DaVinci().UserAlgorithms.append(switchSeq_bu)
        dttsel = FilterSelection("bu_dd0_switch_sel", switchkst_sel_bu, Code=filtercode)
    if spectrum == "M_mst_z":
        switchkst_bu2 = SubstitutePID( 'switchkst_bu2',
                                      Code = "DECTREE('(B+ -> D*(2010)+ D0 K*(892)0) || (B- -> D*(2010)- D0 K*(892)~0)')",
                                      Substitutions = {
                                          'B+ -> D*(2010)+ D0 (Meson -> X+ ^pi-)' : 'K-',
                                          'B+ -> D*(2010)+ D0 (Meson -> ^K+ X-)' : 'pi+',
                                          'B+ -> D*(2010)+ D0 ^K*(892)0' : 'K*(892)~0',
                                          'B- -> D*(2010)- D0 (Meson -> X- ^pi+)' : 'K+',
                                          'B- -> D*(2010)- D0 (Meson -> ^K- X+)' : 'pi-',
                                          'B- -> D*(2010)- D0 ^K*(892)~0' : 'K*(892)0',
                                      })
        switchkst_bu2.MaxChi2PerDoF = -1
        switchkst_sel_bu2 = Selection(
            'switchkst_sel_bu2',
            Algorithm = switchkst_bu2,
            RequiredSelections=[relsel]
        )
        switchSeq_bu2 = SelectionSequence('switchSeq_bu2', TopSelection=switchkst_sel_bu2)
        DaVinci().UserAlgorithms.append(switchSeq_bu2)
        dttsel = FilterSelection("bu_dd0_switch_sel2", switchkst_sel_bu2, Code=filtercode)
    if spectrum != "M_mst_z" and spectrum != "M_m_z":
        dttsel = FilterSelection("dttsels_{}".format(spectrum), relsel, Code=filtercode)
    print("the dttsel for {} is {}".format(spectrum, dttsel))
    DTT = DecayTreeTuple("{0}_{1}_{2}".format(tt.name, spectrum, tt.id))
    DTT.Inputs = [dttsel.outputLocation()]
    DTT.TupleName = "DecayTreeTuple"
    DTT.setDescriptorTemplate(rec_Descriptors[spectrum])
    addtools(DTT, spectrum)
    DaVinci().UserAlgorithms.append(dttsel)
    DaVinci().UserAlgorithms.append(DTT)

#!!!!!!!!!!!!!!!!!!!!11
# GET RID OF THIS FOR GANGA!!!!
# DaVinci().TupleFile = "ntuple.root"
# DaVinci().DataType  = "2016"
# DaVinci().EvtMax = 4000
# DaVinci().Simulation = True
# DaVinci().Lumi = False
# DaVinci().CondDBtag = "sim-20170721-2-vc-mu100"
# DaVinci().DDDBtag = "dddb-20170721-3"
# # # # #!!!!!!!!!!!!!!!!!!!!!!!!!!!11:
# # #
# from GaudiConf import IOHelper
# # #11198006
# IOHelper().inputFiles([
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000082_7.AllStreams.mdst',
# ], clear=True)
# # # # #11198000
# # # #'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00087018/0000/00087018_00001255_7.AllStreams.mdst'
# # ], clear=True)
# IOHelper().inputFiles([
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125739/0000/00125739_00000009_7.AllStreams.mdst'
# ], clear=True)


# #11198006
#'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000082_7.AllStreams.mdst',
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000081_7.AllStreams.mdst'
# '/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000082_7.AllStreams.mdst',
# '/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000083_7.AllStreams.mdst',
# '/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000084_7.AllStreams.mdst',
# '/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000086_7.AllStreams.mdst',
# '/eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000088_7.AllStreams.mdst'
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/MC/2016/ALLSTREAMS.MDST/00125679/0000/00125679_00000082_7.AllStreams.mdst',

# IOHelper().inputFiles([
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRON.MDST/00103102/0000/00103102_00007316_1.bhadron.mdst',
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRON.MDST/00103102/0000/00103102_00009250_1.bhadron.mdst',
# 'root://eoslhcb.cern.ch//eos/lhcb/grid/prod/lhcb/LHCb/Collision16/BHADRON.MDST/00103102/0000/00103102_00007380_1.bhadron.mdst',
# ], clear=True)
