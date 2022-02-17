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
        "BPVVDCHI2":"BPVVDCHI2()",
        "BPVDIRA":"BPVDIRA()",
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
def addtools (t):
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
        #'TupleToolSubmass'
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
        ]
    if mytype != 'Data':
        truth = t.addTupleTool("TupleToolMCTruth")
        truth.addTupleTool("MCTupleToolHierarchy")
        t.addTupleTool("TupleToolMCBackgroundInfo")
        t.addTupleTool("MCTupleToolRedecay")
    t.ToolList += ToolList
    t.B.addTupleTool('TupleToolDecayTreeFitter/dtf')
    t.B.dtf.constrainToOriginVertex = True
    t.B.dtf.Verbose = True
    name = str(t.name())
    print (name, "\n")
    if '_st_' in name:
        t.B.dtf.daughtersToConstrain += ["D*(2010)+", "D*(2010)-", 'D0', 'D~0']
        print ("_st", "\n")
    elif '_s_' in name:
        t.B.dtf.Substitutions = {
            "Beauty -> Charm ^(D+ -> K+ K- pi+) Meson" : "D_s+",
            "Beauty -> Charm ^(D- -> K- K+ pi-) Meson" : "D_s-",
        }
        t.B.dtf.daughtersToConstrain += ['D_s+', 'D_s-', 'D+', 'D-']
        print ("_s", "\n")
    elif '_p_' in name or '_m_' in name:
        t.B.dtf.daughtersToConstrain += ['D+', 'D-', 'D0', 'D~0']
        print ("p", "\n")
    elif '_zz_' in name:
        t.B.dtf.daughtersToConstrain += ['D0', 'D~0']
    elif '_z_' in name:
        t.B.dtf.daughtersToConstrain += ['D+', 'D-']

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
    if 'norm' not in name:
        if '11198030' not in name:
            trig = t.KST.addTupleTool("TupleToolTISTOS")
            trig.TriggerList = triglist
            trig.Verbose = True
            addmidinfo( t.KST )
            addtrkinfo( t.KSTH1 )
            addtrkinfo( t.KSTH2 )
    addtopinfo( t.B )
    sinfo_str = "INFO( {}, -1.0)"
    ## angles are just 1.5, 1.7, 1 every time, so don't bother saving them
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

    t.addTool(TupleToolSubMass,name="TupleToolSubMass")
    t.TupleToolSubMass.SetMax = 4
    t.ToolList += ["TupleToolSubMass"]
##########################################
#Fixes a bug that adds superfulous branches to our ntuple
def addBranches(self, branches):
    'Simplified adding of branches a little bit'
    'takes a dictionary of {branch: decay descriptor}, returns a dictionary of {branch: configurable instances}'
    from Configurables import TupleToolDecay
    if 'Branches' not in dir(self):
        raise TypeError('you\'re trying to add branches to something which doesn\'t support branching, ' + str(type(self)))
    if not isinstance(branches, type({})):
        raise TypeError('expected a dictionary of branches, got a ' + str(type(branches)) + ' instead')

    if self.Branches is None:
        self.Branches = {}

    instances = {}
    for branch in branches:
        # check for whitespace
        #for char in ""\t\r\s":
            #if char in branch:
                #raise NameError('You have tried to add a branch named \'' + branch + '\',which contains whitespace. This is not permitted.')
        self.Branches[branch] = branches[branch]
        self.addTool(TupleToolDecay, branch)
        instances[branch] = getattr(self, branch)

    return instances
def setDescriptorTemplate(self, template):
    'Bored of typing decay descriptors and adding carat symbols?'
    'Use some python string template magic to set your decay descriptor'
    'and define your branches all in one go without excess typing!'
    'from https://gitlab.cern.ch/snippets/147'
    if 'Decay' not in dir(self):
        raise TypeError('You\'re trying to set the decay descriptor of something that doesn\'t have one, ' + str(type(self)))
    if 'Branches' not in dir(self):
        raise TypeError('You\'re trying to define branches on something that doesn\'t support them, ' + str(type(self)))

    from string import Template
    # The argument 'template' is a Python string template
    # e.g. '${D}[D0 -> ${kaon}K- ${pion}pi+]CC'
    # Here ['D', 'kaon', 'pion'] are the branch names you want
    dd = Template(template)

    # This parses the temlate to get the list of branch names,
    # i.e. ['D', 'kaon', 'pion']
    particles = [y[1] if len(y[1]) else y[2] for y in dd.pattern.findall(dd.template) if len(y[1]) or len(y[2])]

    # To form the decay descriptor, we need to mark all the particles
    # except for the top-level particle, which is included by default
    mapping = {p: '^' for p in particles }
    mapping[particles[0]] = ''

    # Make the descriptor
    # '[D0 -> ^K- ^pi+]CC'
    self.Decay = dd.substitute(mapping)

    # Now make the branches
    branches = {}
    for p in particles:
        # Need a version of the descriptor where particle 'p' is marked but nothing else is.
        mapping = {q: '^' if p == q else '' for q in particles }
        branches[p] = dd.substitute(mapping)

    # Finally, add the branches to the DecayTreeTuple
    return self.addBranches(branches)
# DecayTreeTuple.addBranches = addBranches
# MCDecayTreeTuple.addBranches = addBranches
# DecayTreeTuple.setDescriptorTemplate = setDescriptorTemplate
# MCDecayTreeTuple.setDescriptorTemplate = setDescriptorTemplate

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
    "z"       : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dm_decay,dp_decay,kst_decay,dm_decay_cc,dp_decay_cc,kst_decay_cc),
    "zstm"     : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dstm_decay,dp_decay,kst_decay,dstm_decay_cc,dp_decay_cc,kst_decay_cc),
    "zstp"     : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dm_decay,dstp_decay,kst_decay,dm_decay_cc,dstp_decay_cc,kst_decay_cc),
    "zstst"   : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dstm_decay,dstp_decay,kst_decay,dstm_decay_cc,dstp_decay_cc,kst_decay_cc),
    "zz"       : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(d0bar_decay,d0_decay,kst_decay,d0bar_decay_cc,d0_decay_cc,kst_decay_cc),
    "p"       : "(${{B}}(B+ ->{}{}{}) || ${{B}}(B- ->{}{}{}))".format(d0bar_decay,dp_decay,kst_decay,d0bar_decay_cc,dp_decay_cc,kst_decay_cc),
    "m"       : "(${{B}}(B- ->{}{}{}) || ${{B}}(B+ ->{}{}{}))".format(dm_decay, d0_decay, kst_decay, dm_decay_cc, d0_decay_cc, kst_decay_cc),
    "st"     : "(${{B}}(B+ ->{}{}{}) || ${{B}}(B- ->{}{}{}))".format(d0bar_decay,dstp_decay,kst_decay,d0bar_decay_cc,dstp_decay_cc,kst_decay_cc),
    "s"      : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dsm_decay,dp_decay,kst_decay,dsm_decay_cc,dp_decay_cc,kst_decay_cc),
    #"Bs_DsmDstp"    : "(${{B}}(B0 ->{}{}{}) || ${{B}}(B0 ->{}{}{}))".format(dsm_decay,dstp_decay,kst_decay,dsm_decay_cc,dstp_decay_cc,kst_decay_cc),
    "norm7" : "(${{B}}(B+ ->{} {} ${{K}} K+) || ${{B}}(B- ->{} {} ${{K}} K-))".format(d0bar_decay,d04_decay,d0bar_decay_cc,d04_decay_cc),
    "norm8" : "(${{B}}(B0 ->{} {} ${{K}} K+) || ${{B}}(B0 ->{} {} ${{K}} K-))".format(dm_decay,d04_decay,dm_decay_cc,d04_decay_cc),
}
strip_dict = {
    'z'     : 'B02DDKstBeauty2CharmLine',
    's'     : 'B02DDKstBeauty2CharmLine',
    'zst'    : 'B02DstDKstBeauty2CharmLine',
    # 'zstst'  : 'B02DstDstKstBeauty2CharmLine',
    'norm8' : 'B02D0DKD02K3PiBeauty2CharmLine',
    'zz'    : 'B02D0D0KstD02HHD02HHBeauty2CharmLine',
    'p'     : 'B2DD0KstBeauty2CharmLine',
    'm'     : 'B2DD0KstBeauty2CharmLine',
    'st'    : 'B2DstD0KstBeauty2CharmLine',
    'norm7' : 'B2D0D0KD02HHD02K3PiBeauty2CharmLine',
    }

class mytuple():
    def __init__(info, name, id, spectrum):
        info.name = name
        info.id = id
        info.spectrum = spectrum

t = []

t.append(mytuple("01",  "11198000", ["z"]))
t.append(mytuple("02", "11198400", ["z"]))
t.append(mytuple("02", "11198005", ["p"]))

t.append(mytuple("04", "11198401", ["z"]))
t.append(mytuple("04", "11198410", ["p"]))
t.append(mytuple("04", "11198020", ["zz", "st"]))
# t.append(mytuple("04", "11198020", ["st"]))

t.append(mytuple("05", "12197020", ["p"]))
t.append(mytuple("06", "12197410", ["p"]))

t.append(mytuple("07", "12197400",["p"]))
t.append(mytuple("07", "12197022", ["zz", "st"]))
# t.append(mytuple("07", "12197022", ["st"]))

t.append(mytuple("08", "12197401", ["p"]))
t.append(mytuple("08", "12197420", ["zz","st"]))
# t.append(mytuple("08", "12197420", ["st"]))
# t.append(mytuple("8b_g", "12197220", ["zz"]))
# t.append(mytuple("8c_g", "12197220", ["st"]))

t.append(mytuple("09", "11196000", ["zz"]))
t.append(mytuple("10", "11196410", ["zz"]))
# t.append(mytuple("10_g", "11196200", ["zz"]))
#
t.append(mytuple("12", "11196420", ["zz"]))
# t.append(mytuple("12_g", "11196210", ["zz"]))
# t.append(mytuple("12_pig", "11196620", ["zz"]))
#
t.append(mytuple("13", "13198040",["s"]))
t.append(mytuple("14", "13198200", ["s"]))
t.append(mytuple("15", "13198400", ["s"]))
t.append(mytuple("16", "13198600", ["s"]))
#
t.append(mytuple("norm7", "12197009", ["norm7"]))
t.append(mytuple("norm8", "11198030", ["norm8"]))
#
t.append(mytuple("data", "Data", ["z","zz","p","m","st","s"]))

def getthistuple(type):
    for i in t:
        if i.id == mytype:
            return i
tt = getthistuple(mytype)

DaVinci().InputType = myst

if myst == 'mDST':
    DaVinci().RootInTES = '/Event/{}'.format(mystream)
if tt.id != "Data":
    MCDTT = MCDecayTreeTuple("MCDecayTreeTuple")
    if "norm" in tt.name:
        MCDTT.setDescriptorTemplate("${B}[Beauty => ${C}Charm ${C}Charm ${K}K+]CC")
        print("lol")
    else:
        MCDTT.setDescriptorTemplate("${B}[Beauty => ${C}Charm ${C}Charm ${Kst}(K*(892)0 => ${K}K+ ${pi}pi-)]CC")
    MCDTT.ToolList+=["MCTupleToolKinematic","MCTupleToolHierarchy"]
    DaVinci().UserAlgorithms.append(MCDTT)
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
    if myst == 'mDST':
        loc = AutomaticData(Location = 'Phys/{}/Particles'.format(strip_dict[spectrum]))
    if myst == 'DST':
        loc = AutomaticData(Location = '/Event/{}/Phys/{}/Particles'.format(mystream, strip_dict[spectrum]))
    relsel = FilterSelection("relsel_{}".format(spectrum), loc, Preambulo= [rstr.format(strip_dict[spectrum]) for rstr in preamble], Code=sinfo_code)
    DaVinci().UserAlgorithms.append(relsel)
    if spectrum == "m":
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
    else:
        dttsel = FilterSelection("dttsels_{}".format(spectrum), relsel, Code=filtercode)
    DTT = DecayTreeTuple("{0}_{1}_{2}".format(tt.name, spectrum, tt.id))
    DTT.Inputs = [dttsel.outputLocation()]
    DTT.TupleName = "DecayTreeTuple"
    DTT.setDescriptorTemplate(rec_Descriptors[spectrum])
    addtools(DTT)
    DaVinci().UserAlgorithms.append(dttsel)
    DaVinci().UserAlgorithms.append(DTT)

# #!!!!!!!!!!!!!!!!!!!!11
#GET RID OF THIS FOR GANGA!!!!
# DaVinci().TupleFile = "ntuple.root"
# DaVinci().DataType  = "2011"
# DaVinci().EvtMax = 1000
# DaVinci().Simulation = True
# DaVinci().Lumi = False
# DaVinci().CondDBtag = "sim-20170721-2-vc-mu100"
# DaVinci().DDDBtag = "dddb-20170721-3"
# DaVinci().DDDBtag = "dddb-20160318-1"
# DaVinci().CondDBtag = "sim-20160614-1-vc-mu100"
# # #!!!!!!!!!!!!!!!!!!!!!!!!!!!11:
