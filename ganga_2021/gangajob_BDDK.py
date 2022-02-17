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
                           "######## - Specify MC file. i.e, 11198006 \n",
                    choices=['Data',
                             '11198006',
                             '11198400',
                             "11198005",
                             '11198410',
                             '11198401',
                             '13198040',
                             '13198400',
                             '13198200',
                             '13198600',
                             '12197400',
                             '12197401',
                             '12197422',
                             '12197410',
                             '12197024',
                             '12197023',
                             '11198022',
                             '11196414',
                             '11196413',
                             '11196019',
                             '12197008',
                             '11198007',
                             '12197423',
                             '12197045',
                             '11198023'])
parser.add_argument('--YEAR',default='2016',help="do you want 2016, 2017, or 2018)?",
                    choices=['2016','2017','2018'])
parser.add_argument('--POL',default='up',help="up or down?",
                    choices=['Up','Down'])
parser.add_argument('--TESTING', default='yes',help="Are you testing or running for real?",
                    choices=['yes','no'])

args = parser.parse_args()
mytype = args.TYPE
myyear = args.YEAR
mytesting = args.TESTING
mypol = args.POL
myApplication = GaudiExec()
myApplication.directory = "/afs/cern.ch/user/h/hbernste/DaVinci_v45r6"
myApplication.platform ='x86_64-centos7-gcc9-opt'

if mytype != "Data":
    if mytype in ['12197423','12197045','11198023']:
        path = new_paths_mc_09k[f"{myyear}_{mypol}"] + f"/{mytype}/ALLSTREAMS.MDST"
    else:
        path = new_paths_mc[f"{myyear}_{mypol}"] + f"/{mytype}/ALLSTREAMS.MDST"
    myst = "MDST"
    mystream = "AllStreams"
else:
    path = new_paths_data[f"{myyear}_{mypol}"]
    myst = "MDST"
    mystream = "Bhadron"

Files = []
name = mytype+"_"+myyear+"_"+mypol

if mytesting == "no":
    if mytype == "Data":
        j = Job(
                backend=Dirac(),
                name=name,
                parallel_submit=True,
                splitter=SplitByFiles(filesPerJob=10, maxFiles = -1),
                do_auto_resubmit=False,
                outputfiles=[DiracFile(namePattern='*.root')], #has to match what's in optsfile
                )
    if mytype != "Data":
        j = Job(
                backend=Dirac(),
                name=name,
                parallel_submit=True,
                splitter=SplitByFiles(filesPerJob=10, maxFiles = -1),
                do_auto_resubmit=False,
                outputfiles=[DiracFile(namePattern='*.root')], #has to match what's in optsfile
                )
if mytesting == "yes":
    j = Job(
            backend=Dirac(),
            name=name,
            parallel_submit=True,
            splitter=SplitByFiles(filesPerJob=2, maxFiles = 100),
            do_auto_resubmit=False,
            outputfiles=[DiracFile(namePattern='*.root')], #has to match what's in optsfile
            )

j.splitter.ignoremissing = True

extraOpts = 'from Configurables import DaVinci\nDaVinci().EvtMax = -1;'
extraOpts += 'DaVinci().DataType = "'+myyear+'";'
extraOpts += 'DaVinci().PrintFreq = 50000;'
extraOpts += 'DaVinci().TupleFile = "ntuple.root";'

if mytype != 'Data':
        extraOpts += all_dddb_dict['dddb']
        extraOpts += all_dddb_dict[f"{myyear}_{mypol}"]
        extraOpts += 'DaVinci().Simulation = True;'
        extraOpts += 'DaVinci().Lumi = False;'
else:
        extraOpts += 'DaVinci().Simulation = False;'
        extraOpts += 'DaVinci().Lumi = True;'

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
myApplication.options = ["/afs/cern.ch/user/h/hbernste/bddk/ganga_2021/options_class_BDDK.py"]
j.application = myApplication
j.inputdata = BKQuery(path).getDataset()
j.submit()

os.remove('steering.py')#no need to leave it behind
