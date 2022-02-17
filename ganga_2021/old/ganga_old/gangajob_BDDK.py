#-----------------------------------------------
# Author   : Harris Bernstein
# Date     : WIP
# Comments : B0 -> D(*)0,+,s D(*)0,-,s K*
# Using    : Run1 and Run2 Data and available MC
#-----------------------------------------------
"""Job submission script for running over b quark particle to -> (Charm Charm Kst) MC and Data BF"""
import os, sys, argparse
from argparse import RawTextHelpFormatter
sys.path.append(os.getcwd())
from paths_CharmCharmKst import *

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
parser.add_argument('--TYPE',
                    help = "Data - Run over real data \n"
                           "######## - Specify MC file. i.e, 11198000 \n",
                    choices=["Data",
                             '11198000',
                             '11198020',
                             '11198005',
                             '11198400',
                             '11198401',
                             '11198410',
                             "12197400",
                             "12197022",
                             "12197401",
                             "12197420",
                             "11196000",
                             "11196410",
                             "11196420",
                             "13198040",
                             "13198200",
                             "13198400",
                             "13198600",
                             "11198030",
                             "12197009",
                             "12197220",
                             "11196200",
                             "11196210",
                             "11196620"])
parser.add_argument('--YEAR',default='2016',help="do you want 2011, 2012, 2015, 2016, 2017, or 2018)?",
                    choices=['2011','2012','2015','2016','2017','2018'])
parser.add_argument('--TESTING', default='yes',help="Are you testing or running for real?",
                    choices=['yes','no'])
args = parser.parse_args()
mytype = args.TYPE
myyear = args.YEAR
mytesting = args.TESTING
myApplication = GaudiExec()
myApplication.directory = "/afs/cern.ch/user/h/hbernste/DaVinci_v50r6"
myApplication.platform ='x86_64-centos7-gcc8-opt'

if mytype == 'Data':
    paths = path_dict[mytype][myyear]

if mytype != 'Data':
    print(mytype)
    paths = path_dict["mc"][mytype][myyear][0]
    mypol1 = paths[0].split("Mag")[1].split("-")[0]
    mypol2 = paths[1].split("Mag")[1].split("-")[0]
    if mypol1 != 'Up':
        print("\n")
        print("BAAAAAD POL U")
        print("\n")
    if mypol2 != 'Down':
        print("\n")
        print("BAAAAAD POL D")
        print("\n")
    if mypol1 == mypol2:
        print("\n")
        print("BAAAAAD POL")
        print("\n")
if mytesting == 'yes':
    paths = [paths[0]]
    print (paths)

for path in paths:
    Files = []
    if mytype == 'Data':
        mypol = path.split("Mag")[1].split("/")[0]
    else:
        mypol = path.split("Mag")[1].split("-")[0]
    lastcontain = path.split("/")[-1]
    if 'STRIP' in lastcontain:
        mystream = lastcontain.split('.')[0]+'.'+lastcontain.split('.')[1]
        myst = lastcontain.split('.')[2]
    else:
        mystream = lastcontain.split('.')[0]
        myst = lastcontain.split('.')[1]

    if mystream == "ALLSTREAMS":
        mystream = "AllStreams"
    if mystream == "BHADRON":
        mystream = "Bhadron"
    if mystream == "B2DDKST.STRIP":
        mystream = "b2ddkst.Strip"

    if myst == 'DST':
        myst = 'DST'
    if myst == 'MDST':
        myst = 'mDST'

    name = mytype+"_"+myyear+"_"+mypol
    print(name)

    j = Job(
        backend=Dirac(),
        name=name,
        parallel_submit=True,
        splitter=SplitByFiles(filesPerJob=4, maxFiles = -1),
        do_auto_resubmit=False,
        outputfiles=[DiracFile(namePattern='*.root')], #has to match what's in optsfile
        )

    j.splitter.ignoremissing = True
    #j.backend.settings['BannedSites'] = ['LCG.RAL.uk']

    if mytesting == 'yes':
        j.splitter=SplitByFiles(filesPerJob=1, maxFiles=200)
        extraOpts = 'from Configurables import DaVinci\nDaVinci().EvtMax = -1;'
    elif mytesting == 'no':
        #j.splitter=SplitByFiles(filesPerJob=15, maxFiles = -1)
        extraOpts = 'from Configurables import DaVinci\nDaVinci().EvtMax = -1;'
        if mytype != 'Data':
            j.splitter=SplitByFiles(filesPerJob=2, maxFiles = -1)

    extraOpts += 'DaVinci().DataType = "'+myyear+'";'
    extraOpts += 'DaVinci().PrintFreq = 50000;'
    extraOpts += 'DaVinci().TupleFile = "ntuple.root";'

    if mytype != 'Data':
        if mypol == 'Up':
            extraOpts += path_dict["mc"][mytype][myyear][1][0]
            extraOpts += path_dict["mc"][mytype][myyear][2][1]
        if mypol == 'Down':
            extraOpts += path_dict["mc"][mytype][myyear][1][0]
            extraOpts += path_dict["mc"][mytype][myyear][2][1]
        extraOpts += 'DaVinci().Simulation = True;'
        extraOpts += 'DaVinci().Lumi = False;'
    else:
        extraOpts += 'DaVinci().Simulation = False;'
        extraOpts += 'DaVinci().Lumi = True;'

    print("Adding the following options " + extraOpts)

    myApplication.extraOpts = extraOpts
    theline = 'mytype = "'+mytype+'"\n'
    theline += 'mystream = "'+mystream+'"\n'
    theline += 'myst = "'+myst+'"\n'

    if os.path.isfile('steering.py'):
        os.remove('steering.py')#remove old copies left behind by accident
    f = open('steering.py', 'w')
    f.write(theline)
    f.close()
    j.inputfiles += ['steering.py']

    myApplication.options = ["/afs/cern.ch/user/h/hbernste/bddk/ganga/options_class_BDDK.py"]
    j.application = myApplication
    j.inputdata = BKQuery(path).getDataset()
    j.submit()

os.remove('steering.py')#no need to leave it behind
